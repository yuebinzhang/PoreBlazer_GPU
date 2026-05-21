!--------------------------------------------------------------------------------------
! GPU-accelerated subroutines for PoreBlazer v4.0
!
! Uses OpenACC directives with NVIDIA HPC SDK (nvfortran)
! Accelerates the two main computational hot spots:
!   1. lattice_calculations  — 3D grid × atoms distance calculations
!   2. surface_area          — Monte Carlo surface sampling
!
! Compile with: nvfortran -acc -gpu=cc89 -O3 -Minfo=accel ...
!--------------------------------------------------------------------------------------

!--------------------------------------------------------------------------------------
! GPU-accelerated lattice_calculations
!
! This is the dominant hot spot (~95% of runtime). The algorithm:
!   For each cubelet in the 3D grid (ncubesx × ncubesy × ncubesz):
!     For each atom in the structure:
!       - compute minimum-image distance
!       - check if cubelet center is inside an atom (overlap)
!       - track nearest atom distance
!       - compute LJ interaction with helium probe
!
! GPU strategy: Parallelize the 3D cubelet loop. Each GPU thread handles one
! cubelet and loops over all atoms sequentially (coalesced reads from atom data).
! A compaction pass on CPU fills the variable-length output arrays (g_cubes,
! he_cubes, n_cubes) to avoid atomic operations in the GPU kernel.
!
! For orthorhombic cells (the common case for MOFs/zeolites), the minimum-image
! calculation is inlined directly. Non-orthorhombic cells fall back to CPU.
!--------------------------------------------------------------------------------------
subroutine lattice_calculations_gpu
    use parameters
    use atoms
    use adsorbent
    use lattice
    use defaults, only: rdbl
    implicit none

    ! Cell geometry — extracted from fcell as plain types for GPU
    real*8  :: ell(3)
    logical :: ortho

    ! Loop variables
    integer :: i, j, k, l, icount

    ! Per-cubelet temporaries
    real*8  :: xc, yc, zc
    real*8  :: d1, d2, d3, rdist2, rdist2_ref
    real*8  :: rdist_surface, rdist_surface_ref, lj_local
    real*8  :: sig2_rdist2, rdist6, rdist12, lj_energy
    logical :: overlap
    integer :: imin

    ! Compaction counters
    integer :: ng, nhe, nn

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! GPU-accelerated lattice calculations (OpenACC)         !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    ! Check if cell is orthorhombic — non-ortho falls back to CPU
    ortho = fcell%orthoflag
    if(.not. ortho) then
        write(*,*) "  Non-orthorhombic cell detected — using CPU version"
        call lattice_calculations
        return
    end if

    ! Extract cell dimensions
    ell = fcell%ell

    ! Initialize output arrays
    lattice_space = 0; lattice_space_he = 0; lattice_space_n = 0
    lattice_rdist2 = 0.0d0; lattice_lj_he = 0.0d0
    g_cubes = 0; he_cubes = 0; n_cubes = 0

    !####################################################################
    ! GPU KERNEL: Parallel over all cubelets
    ! No compaction arrays — use lattice_space as marker.
    ! lattice_index computed on CPU after (saves 3 × ntot × 4 bytes VRAM).
    !####################################################################
    !$acc enter data copyin(coords_temp, asigma2, asigma2_he, aeps_he, &
    !$acc&   asigma, asigma2_n, atype, ell, hicut2, cube_size, &
    !$acc&   ncubesx, ncubesy, ncubesz, natoms) &
    !$acc&   create(lattice_space, lattice_rdist2, lattice_space_he, &
    !$acc&   lattice_space_n, lattice_lj_he)

    !$acc parallel loop collapse(3) gang vector &
    !$acc& private(xc, yc, zc, d1, d2, d3, rdist2, rdist2_ref, overlap, &
    !$acc&   rdist_surface, rdist_surface_ref, imin, lj_local, &
    !$acc&   sig2_rdist2, rdist6, rdist12, lj_energy) &
    !$acc& firstprivate(ell, cube_size, hicut2) &
    !$acc& present(coords_temp, asigma2, asigma2_he, aeps_he, asigma, &
    !$acc&   asigma2_n, atype, lattice_space, lattice_rdist2, &
    !$acc&   lattice_space_he, lattice_space_n, lattice_lj_he, &
    !$acc&   ncubesx, ncubesy, ncubesz)
    do l=1, ncubesz
        do k=1, ncubesy
            do j=1, ncubesx

                !--- Cubelet center coordinates (slanted frame) ---
                xc = dble(j-1)*cube_size + 0.5d0*cube_size
                yc = dble(k-1)*cube_size + 0.5d0*cube_size
                zc = dble(l-1)*cube_size + 0.5d0*cube_size

                !--- Initialize per-cubelet accumulators ---
                overlap = .false.
                rdist2_ref = huge(0.0d0)
                rdist_surface_ref = huge(0.0d0)
                lj_local = 0.0d0
                imin = 1

                !--- Loop over all atoms (sequential per cubelet) ---
                do i=1, natoms
                    ! Minimum-image distance (orthorhombic fast path)
                    d1 = xc - coords_temp(1, i)
                    d2 = yc - coords_temp(2, i)
                    d3 = zc - coords_temp(3, i)
                    d1 = d1 - ell(1)*anint(d1/ell(1))
                    d2 = d2 - ell(2)*anint(d2/ell(2))
                    d3 = d3 - ell(3)*anint(d3/ell(3))
                    rdist2 = d1*d1 + d2*d2 + d3*d3

                    !--- Overlap check: cube center inside an atom? ---
                    if(rdist2 < 0.25d0*asigma2(atype(i))) then
                        overlap = .true.
                        exit
                    end if

                    !--- Track nearest atom ---
                    if(rdist2 < rdist2_ref) then
                        rdist2_ref = rdist2
                        imin = i
                    end if

                    !--- Track shortest surface distance ---
                    rdist_surface = dsqrt(rdist2) - 0.5d0*asigma(atype(i))
                    if(rdist_surface < rdist_surface_ref) then
                        rdist_surface_ref = rdist_surface
                    end if

                    !--- LJ interaction with helium probe (within cutoff) ---
                    if(rdist2 < hicut2) then
                        sig2_rdist2 = asigma2_he(atype(i))/rdist2
                        rdist6 = sig2_rdist2*sig2_rdist2*sig2_rdist2
                        rdist12 = rdist6*rdist6
                        lj_energy = aeps_he(atype(i))*(rdist12 - rdist6)
                        lj_local = lj_local + lj_energy
                    end if
                end do  ! i = 1, natoms

                !--- Store LJ energy for this cubelet ---
                icount = (l-1)*ncubesx*ncubesy + (k-1)*ncubesx + j
                lattice_lj_he(icount) = lj_local

                !--- Skip if cubelet center overlaps with an atom ---
                if(overlap) cycle

                !--- Register as geometrically accessible (point probe) ---
                lattice_space(j, k, l) = 1
                lattice_rdist2(j, k, l) = rdist_surface_ref*rdist_surface_ref

                !--- Helium accessibility ---
                if(rdist2_ref > 0.25d0*asigma2_he(atype(imin))) then
                    lattice_space_he(j, k, l) = 1
                end if

                !--- Nitrogen accessibility ---
                if(rdist2_ref > asigma2_n(atype(imin))) then
                    lattice_space_n(j, k, l) = 1
                end if

            end do  ! j
        end do  ! k
    end do  ! l
    !$acc end parallel loop

    !####################################################################
    ! Data transfer: copy results back from GPU
    !####################################################################
    !$acc exit data copyout(lattice_space, lattice_rdist2, lattice_space_he, &
    !$acc&   lattice_space_n, lattice_lj_he)

    !####################################################################
    ! CPU compaction: use lattice_space markers + recompute icount
    ! This avoids storing 3 large g_valid/he_valid/n_valid arrays on GPU.
    !####################################################################
    ng = 0; nhe = 0; nn = 0
    do l=1, ncubesz
        do k=1, ncubesy
            do j=1, ncubesx
                icount = (l-1)*ncubesx*ncubesy + (k-1)*ncubesx + j
                if(lattice_space(j,k,l) == 1) then
                    ng = ng + 1
                    g_cubes(ng) = icount
                end if
                if(lattice_space_he(j,k,l) == 1) then
                    nhe = nhe + 1
                    he_cubes(nhe) = icount
                end if
                if(lattice_space_n(j,k,l) == 1) then
                    nn = nn + 1
                    n_cubes(nn) = icount
                end if
            end do
        end do
    end do

    ng_cubes = ng
    nhe_cubes = nhe
    nn_cubes = nn

    !--- Compute lattice_index on CPU (was done on GPU before, moved to save VRAM) ---
    do icount = 1, ntot
        l = (icount - 1) / (ncubesx * ncubesy) + 1
        k = modulo((icount - 1) / ncubesx, ncubesy) + 1
        j = modulo(icount - 1, ncubesx) + 1
        lattice_index(1, icount) = j
        lattice_index(2, icount) = k
        lattice_index(3, icount) = l
    end do

    !--- Allocate percolation arrays ---
    allocate(PA1(ng_cubes), PA2(ng_cubes), PA3(ng_cubes), PA4(ng_cubes))
    PA1 = 0.0d0; PA2 = 0; PA3 = 0; PA4 = 0

    write(*,'(a,i10)') "  Cubelets accessible to point probe:  ", ng_cubes
    write(*,'(a,i10)') "  Cubelets accessible to helium:       ", nhe_cubes
    write(*,'(a,i10)') "  Cubelets accessible to nitrogen:     ", nn_cubes
    write(*,*)
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! GPU lattice calculations complete                      !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return
end subroutine lattice_calculations_gpu


!--------------------------------------------------------------------------------------
! GPU-accelerated surface_area
!
! For each atom, nsample random surface points are generated and tested against
! all other atoms for overlap. The fraction of successful trials gives the
! accessible surface area.
!
! GPU strategy: Parallelize over atoms. Each GPU thread:
!   1. Generates random directions using a thread-local LCG
!   2. Creates test points and applies PBC
!   3. Checks overlap against all other atoms (inner reduction)
!   4. Accumulates successful trials
!
! Uses a thread-safe inline LCG random number generator (no external library).
! The cell is assumed orthorhombic (most common case); falls back to CPU otherwise.
!--------------------------------------------------------------------------------------
subroutine surface_area_gpu
    use parameters
    use atoms
    use adsorbent
    use lattice
    use defaults, only: rdbl, pi
    use results
    implicit none

    ! Local variables
    real*8  :: ell(3), hcut2_loc
    integer :: i, j, k, nx, ny, nz
    real*8  :: phi, costheta, theta
    real*8  :: rdist2, sfrac, sjreal, stotalreduced, vol
    real*8  :: d1, d2, d3
    real*8  :: xp, yp, zp
    logical :: deny
    logical :: ortho

    ! Thread-local RNG seed
    integer :: seed

    ! External GPU random function
    external gpu_rand
    real*8 :: gpu_rand

    ! Per-atom counters (output from GPU kernel)
    integer, allocatable :: ncount_arr(:)
    real*8  :: stotal_gpu

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! GPU-accelerated surface area (OpenACC)                 !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    ! Fall back to CPU for non-orthorhombic cells
    ortho = fcell%orthoflag
    if(.not. ortho) then
        write(*,*) "  Non-orthorhombic cell detected — using CPU version"
        call surface_area
        return
    end if

    ! Quick return if no nitrogen-accessible points
    if(nn_cubes == 0) then
        stotal = 0.0d0
        vol = fundcell_getvolume(fcell)
        stotalreduced = stotal/(vol)*1.0d4
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
            ' surface area in A^2:                 ', stotal
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
            ' surface area per volume in m^2/cm^3: ', stotalreduced
        write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
            ' surface area per mass in m^2/g:      ', &
            stotalreduced / (smass/(0.6022141d0*vol))
        write(*,*) "  No nitrogen accessible surface area"
        write(*,*)
        write(*,*) "!-------------------------------------------------------!"
        write(*,*) "! GPU surface area calculations complete                 !"
        write(*,*) "!-------------------------------------------------------!"
        write(*,*)
        write(100,*) property
        write(100,'(a,f12.2)') "S_AC_A^2 ", stotal
        write(100,'(a,f12.2)') "S_AC_m^2/cm^3 ", stotalreduced
        write(100,'(a,f12.2)') "S_AC_m^2/g ", &
            stotalreduced / (smass/(0.6022141d0*vol))
        return
    end if

    ! Extract cell info
    ell = fcell%ell
    hcut2_loc = hicut2

    ! Per-atom trial counter
    allocate(ncount_arr(natoms))
    ncount_arr = 0

    !####################################################################
    ! GPU KERNEL: Parallel over atoms
    !####################################################################
    !$acc enter data copyin(coords, asigma_n, asigma2_n, atype, ell, &
    !$acc&   coeff_surface, coeff_surface2, cube_size, lattice_space_n, &
    !$acc&   ncubesx, ncubesy, ncubesz, natoms, nsample) &
    !$acc&   create(ncount_arr)

    ! Thread-local LCG random number generator
    ! Each thread uses seed derived from atom index for reproducibility
    ! LCG: x_{n+1} = (1664525 * x_n + 1013904223) mod 2^31

    !$acc parallel loop gang vector &
    !$acc& private(j, phi, costheta, theta, xp, yp, zp, d1, d2, d3, &
    !$acc&   rdist2, deny, nx, ny, nz, seed) &
    !$acc& firstprivate(ell, coeff_surface, coeff_surface2, cube_size, &
    !$acc&   natoms, nsample, ncubesx, ncubesy, ncubesz) &
    !$acc& present(coords, asigma_n, asigma2_n, atype, ncount_arr, &
    !$acc&   lattice_space_n)
    do i=1, natoms
        ! Initialize thread-local RNG seed
        seed = 1 + i * 694847539
        ncount_arr(i) = 0

        do j=1, nsample
            !--- Generate random unit vector ---
            phi  = pi - gpu_rand(seed) * 2.0d0 * pi
            costheta = 1.0d0 - gpu_rand(seed) * 2.0d0
            theta = acos(costheta)

            !--- Test point on atom i surface ---
            xp = sin(theta) * cos(phi) * (coeff_surface * asigma_n(atype(i))) &
                 + coords(1, i)
            yp = sin(theta) * sin(phi) * (coeff_surface * asigma_n(atype(i))) &
                 + coords(2, i)
            zp = costheta * (coeff_surface * asigma_n(atype(i))) &
                 + coords(3, i)

            !--- Apply PBC (orthorhombic) ---
            d1 = xp - ell(1)*anint(xp/ell(1))
            d2 = yp - ell(2)*anint(yp/ell(2))
            d3 = zp - ell(3)*anint(zp/ell(3))

            if(d1 < 0.0d0) d1 = d1 + ell(1)
            if(d1 >= ell(1)) d1 = d1 - ell(1)
            if(d2 < 0.0d0) d2 = d2 + ell(2)
            if(d2 >= ell(2)) d2 = d2 - ell(2)
            if(d3 < 0.0d0) d3 = d3 + ell(3)
            if(d3 >= ell(3)) d3 = d3 - ell(3)
            xp = d1; yp = d2; zp = d3

            !--- Locate cubelet containing test point ---
            nx = int(xp / cube_size) + 1
            ny = int(yp / cube_size) + 1
            nz = int(zp / cube_size) + 1

            ! Periodic boundary for cubelet index
            if(nx > ncubesx) nx = nx - (nx/ncubesx)*ncubesx
            if(nx < 1) nx = nx + (1-(nx/ncubesx))*ncubesx
            if(ny > ncubesy) ny = ny - (ny/ncubesy)*ncubesy
            if(ny < 1) ny = ny + (1-(ny/ncubesy))*ncubesy
            if(nz > ncubesz) nz = nz - (nz/ncubesz)*ncubesz
            if(nz < 1) nz = nz + (1-(nz/ncubesz))*ncubesz

            !--- Reject if cubelet is not nitrogen-accessible ---
            if(lattice_space_n(nx, ny, nz) < 1) cycle

            !--- Overlap test against all other atoms ---
            deny = .false.
            do k=1, natoms
                if(k == i) cycle
                d1 = xp - coords(1, k)
                d2 = yp - coords(2, k)
                d3 = zp - coords(3, k)
                d1 = d1 - ell(1)*anint(d1/ell(1))
                d2 = d2 - ell(2)*anint(d2/ell(2))
                d3 = d3 - ell(3)*anint(d3/ell(3))
                rdist2 = d1*d1 + d2*d2 + d3*d3

                if(rdist2 < coeff_surface2 * asigma2_n(atype(k))) then
                    deny = .true.
                    exit
                end if
            end do

            if(deny) cycle
            ncount_arr(i) = ncount_arr(i) + 1
        end do  ! j = 1, nsample
    end do  ! i = 1, natoms
    !$acc end parallel loop

    !####################################################################
    ! Data transfer: copy results back
    !####################################################################
    !$acc exit data copyout(ncount_arr)

    !####################################################################
    ! CPU reduction: accumulate surface area
    !####################################################################
    stotal = 0.0d0
    do i=1, natoms
        if(ncount_arr(i) > 0) then
            sfrac = dble(ncount_arr(i)) / dble(nsample)
            sjreal = 4.0d0 * pi * coeff_surface2 * asigma2_n(atype(i)) * sfrac
            stotal = stotal + sjreal
        end if
    end do

    deallocate(ncount_arr)

    !--- Report results ---
    vol = fundcell_getvolume(fcell)
    stotalreduced = stotal / (vol) * 1.0d4

    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
        ' surface area in A^2:                 ', stotal
    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
        ' surface area per volume in m^2/cm^3: ', stotalreduced
    write(*,'(2a,f12.2)') " "//Trim(adjustl(property)), &
        ' surface area per mass in m^2/g:      ', &
        stotalreduced / (smass/(0.6022141d0*vol))
    write(*,*)
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! GPU surface area calculations complete                 !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    write(100,*) property
    write(100,'(a,f12.2)') "S_AC_A^2 ", stotal
    write(100,'(a,f12.2)') "S_AC_m^2/cm^3 ", stotalreduced
    write(100,'(a,f12.2)') "S_AC_m^2/g ", &
        stotalreduced / (smass/(0.6022141d0*vol))

    return
end subroutine surface_area_gpu


!--------------------------------------------------------------------------------------
! GPU-compatible random number generator
! Simple LCG (Linear Congruential Generator)
! Thread-safe: each thread has its own seed
! Returns uniform random in [0, 1)
!--------------------------------------------------------------------------------------
function gpu_rand(seed) result(r)
    implicit none
    !$acc routine vector
    integer, intent(inout) :: seed
    real*8  :: r
    integer :: new_seed

    ! LCG parameters (Numerical Recipes)
    ! x_{n+1} = (1664525 * x_n + 1013904223) mod 2^31
    new_seed = 1664525 * seed + 1013904223

    ! Keep in positive range
    if(new_seed < 0) new_seed = new_seed + 2147483647
    seed = mod(new_seed, 2147483647)
    if(seed < 0) seed = seed + 2147483647

    r = dble(seed) / 2147483647.0d0
    if(r < 0.0d0) r = 0.0d0
    if(r >= 1.0d0) r = 0.9999999999d0

    return
end function gpu_rand
