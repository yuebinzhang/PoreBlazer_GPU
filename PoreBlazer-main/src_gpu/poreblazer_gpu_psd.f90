!--------------------------------------------------------------------------------------
! GPU-accelerated PSD Monte Carlo — pore_distribution
!
! The PA-search algorithm:
!   1. Fill PA1-PA4 from accessible cubelets (same as CPU)
!   2. Sort PA1-PA4 by lattice_rdist2 ascending
!   3. Launch GPU kernel: each thread runs one trial
!   4. Each trial: pick random cubelet, search PA entries from LARGEST to
!      SMALLEST to find largest enclosing sphere, record histogram bin
!   5. CPU reverse cumulative sum to build PSD curve
!
! For orthorhombic cells only; falls back to CPU for non-ortho.
! Uses local array copies to avoid OpenACC allocatable issues.
!--------------------------------------------------------------------------------------
subroutine pore_distribution_gpu
    use parameters
    use atoms
    use adsorbent
    use lattice
    use distributions
    use results
    use fundcell, only: fundcell_getvolume
    use defaults, only: rdbl
    implicit none

    integer :: i, j, k, l, m, bin, nc_local, isite, ncycles, ivis, icount
    integer :: nx, ny, nz, nx1, ny1, nz1
    real*8  :: sigma2_ref, sigma_ref, rdist2, d1, d2, d3
    real*8  :: xc, yc, zc, xc1, yc1, zc1
    real*8  :: rn
    integer*8 :: seed, new_seed
    integer*8 :: t1, t2, t_rate, t_max
    real*8  :: elapsed

    ! Local copies of module arrays for GPU
    real*8, allocatable  :: lcl_PA1(:)
    integer, allocatable :: lcl_PA2(:), lcl_PA3(:), lcl_PA4(:)
    integer, allocatable :: lcl_g_cubes(:), lcl_n_cubes(:)
    real*4, allocatable  :: lcl_psd_cumul(:)

    call system_clock(t1, t_rate, t_max)

    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Starting GPU-accelerated PSD calculations (OpenACC)    !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    ncycles = ncycles_psd
    if(ncycles <= 0) then
        write(*,*) "  PSD skipped (ncycles_psd <= 0)"
        return
    end if

    !--- Fill PA1-PA4 from accessible cubelets (same fill as CPU version) ---
    ivis = 0
    if(property == "Total ") then
        do l=1, ncubesz
            do k=1, ncubesy
                do j=1, ncubesx
                    if(lattice_space(j,k,l) < 1) cycle
                    ivis = ivis + 1
                    PA1(ivis) = lattice_rdist2(j,k,l)
                    PA2(ivis) = j
                    PA3(ivis) = k
                    PA4(ivis) = l
                end do
            end do
        end do
        call sort(ng_cubes, PA1, PA2, PA3, PA4)
        nc_local = ng_cubes
    else
        do l=1, ncubesz
            do k=1, ncubesy
                do j=1, ncubesx
                    if(lattice_space_n(j,k,l) < 1) cycle
                    ivis = ivis + 1
                    PA1(ivis) = lattice_rdist2(j,k,l)
                    PA2(ivis) = j
                    PA3(ivis) = k
                    PA4(ivis) = l
                end do
            end do
        end do
        call sort(nn_cubes, PA1, PA2, PA3, PA4)
        nc_local = nn_cubes
    end if

    if(nc_local == 0) then
        write(*,*) "  No accessible cubelets — PSD is all zeros"
        psd_cumul = 0.0; psd = 0.0
        return
    end if

    write(*,'(a,i10,a,i6)') "  PSD: nc_local=", nc_local, "  ncycles=", ncycles

    !--- Copy to local arrays for GPU ---
    allocate(lcl_PA1(nc_local), source=PA1(1:nc_local))
    allocate(lcl_PA2(nc_local), source=PA2(1:nc_local))
    allocate(lcl_PA3(nc_local), source=PA3(1:nc_local))
    allocate(lcl_PA4(nc_local), source=PA4(1:nc_local))
    if(property == "Total ") then
        allocate(lcl_g_cubes(ng_cubes), source=g_cubes(1:ng_cubes))
    else
        allocate(lcl_n_cubes(nn_cubes), source=n_cubes(1:nn_cubes))
    end if
    allocate(lcl_psd_cumul(nbins+100))
    lcl_psd_cumul = 0.0

    !######################################################################
    ! GPU KERNEL: Parallel Monte Carlo over trials
    !######################################################################

    if(property == "Total ") then
        !$acc data copyin(lcl_PA1, lcl_PA2, lcl_PA3, lcl_PA4, lcl_g_cubes) &
        !$acc&      create(lcl_psd_cumul)

        !$acc parallel loop gang vector_length(32) &
        !$acc& private(isite, nx, ny, nz, nx1, ny1, nz1, icount, &
        !$acc&   sigma2_ref, sigma_ref, rdist2, d1, d2, d3, &
        !$acc&   xc, yc, zc, xc1, yc1, zc1, m, bin, rn, seed, new_seed) &
        !$acc& firstprivate(ncycles, nc_local, ng_cubes, &
        !$acc&   cube_size, nbins, ncubesx, ncubesy, ncubesz) &
        !$acc& present(lcl_PA1, lcl_PA2, lcl_PA3, lcl_PA4, &
        !$acc&   lcl_psd_cumul, lcl_g_cubes)
        do i = 1, ncycles
            sigma2_ref = 0.0d0
            seed = 1_8 + i * 694847539_8
            new_seed = mod(1664525_8 * seed + 1013904223_8, 2147483647_8)
            seed = new_seed
            rn = dble(seed) / 2147483647.0d0
            isite = int(rn * dble(ng_cubes)) + 1
            if(isite > ng_cubes) isite = ng_cubes
            icount = lcl_g_cubes(isite)
            nz = (icount - 1) / (ncubesx * ncubesy) + 1
            ny = modulo((icount - 1) / ncubesx, ncubesy) + 1
            nx = modulo(icount - 1, ncubesx) + 1
            xc = dble(nx-1)*cube_size + 0.5d0*cube_size
            yc = dble(ny-1)*cube_size + 0.5d0*cube_size
            zc = dble(nz-1)*cube_size + 0.5d0*cube_size
            do m = nc_local, 1, -1
                nx1 = lcl_PA2(m); ny1 = lcl_PA3(m); nz1 = lcl_PA4(m)
                xc1 = dble(nx1-1)*cube_size + 0.5d0*cube_size
                yc1 = dble(ny1-1)*cube_size + 0.5d0*cube_size
                zc1 = dble(nz1-1)*cube_size + 0.5d0*cube_size
                d1 = xc - xc1; d2 = yc - yc1; d3 = zc - zc1
                rdist2 = d1*d1 + d2*d2 + d3*d3
                if(rdist2 > lcl_PA1(m)) then
                    cycle
                else
                    sigma2_ref = lcl_PA1(m)
                    sigma_ref = sqrt(sigma2_ref)
                    bin = int(2.0d0*sigma_ref / binsize) + 1
                    if(bin > nbins+100) bin = nbins+100
                    if(bin < 1) bin = 1
                    !$acc atomic update
                    lcl_psd_cumul(bin) = lcl_psd_cumul(bin) + 1.0
                    exit
                end if
            end do
        end do
        !$acc end parallel loop

        !$acc update host(lcl_psd_cumul)
        !$acc end data
    else
        !$acc data copyin(lcl_PA1, lcl_PA2, lcl_PA3, lcl_PA4, lcl_n_cubes) &
        !$acc&      create(lcl_psd_cumul)

        !$acc parallel loop gang vector_length(32) &
        !$acc& private(isite, nx, ny, nz, nx1, ny1, nz1, icount, &
        !$acc&   sigma2_ref, sigma_ref, rdist2, d1, d2, d3, &
        !$acc&   xc, yc, zc, xc1, yc1, zc1, m, bin, rn, seed, new_seed) &
        !$acc& firstprivate(ncycles, nc_local, nn_cubes, &
        !$acc&   cube_size, nbins, ncubesx, ncubesy, ncubesz) &
        !$acc& present(lcl_PA1, lcl_PA2, lcl_PA3, lcl_PA4, &
        !$acc&   lcl_psd_cumul, lcl_n_cubes)
        do i = 1, ncycles
            if(nn_cubes == 0) cycle
            sigma2_ref = 0.0d0
            seed = 1_8 + i * 694847539_8
            new_seed = mod(1664525_8 * seed + 1013904223_8, 2147483647_8)
            seed = new_seed
            rn = dble(seed) / 2147483647.0d0
            isite = int(rn * dble(nn_cubes)) + 1
            if(isite > nn_cubes) isite = nn_cubes
            icount = lcl_n_cubes(isite)
            nz = (icount - 1) / (ncubesx * ncubesy) + 1
            ny = modulo((icount - 1) / ncubesx, ncubesy) + 1
            nx = modulo(icount - 1, ncubesx) + 1
            xc = dble(nx-1)*cube_size + 0.5d0*cube_size
            yc = dble(ny-1)*cube_size + 0.5d0*cube_size
            zc = dble(nz-1)*cube_size + 0.5d0*cube_size
            do m = nc_local, 1, -1
                nx1 = lcl_PA2(m); ny1 = lcl_PA3(m); nz1 = lcl_PA4(m)
                xc1 = dble(nx1-1)*cube_size + 0.5d0*cube_size
                yc1 = dble(ny1-1)*cube_size + 0.5d0*cube_size
                zc1 = dble(nz1-1)*cube_size + 0.5d0*cube_size
                d1 = xc - xc1; d2 = yc - yc1; d3 = zc - zc1
                rdist2 = d1*d1 + d2*d2 + d3*d3
                if(rdist2 > lcl_PA1(m)) then
                    cycle
                else
                    sigma2_ref = lcl_PA1(m)
                    sigma_ref = sqrt(sigma2_ref)
                    bin = int(2.0d0*sigma_ref / binsize) + 1
                    if(bin > nbins+100) bin = nbins+100
                    if(bin < 1) bin = 1
                    !$acc atomic update
                    lcl_psd_cumul(bin) = lcl_psd_cumul(bin) + 1.0
                    exit
                end if
            end do
        end do
        !$acc end parallel loop

        !$acc update host(lcl_psd_cumul)
        !$acc end data
    end if

    !######################################################################
    ! Copy results back and build cumulative histogram
    !######################################################################

    ! Build cumulative histogram (reverse cumulative sum)
    do bin = nbins+99, 1, -1
        lcl_psd_cumul(bin) = lcl_psd_cumul(bin) + lcl_psd_cumul(bin+1)
    end do

    psd_cumul = lcl_psd_cumul
    ffv = psd_cumul(1) / real(ncycles, 4)
    if(property == "Total ") then
        ffv = ffv * dble(ng_cubes) / dble(ntot)
    else
        ffv = ffv * dble(nn_cubes) / dble(ntot)
    end if
    free_volume = ffv * fundcell_getvolume(fcell)

    ! Write output files
    call write_psd_output

    call system_clock(t2)
    elapsed = dble(t2-t1)/dble(t_rate)
    write(*,'(a,f8.3,a)') "  [TIMING] PSD GPU took ", elapsed, " seconds"

    deallocate(lcl_PA1, lcl_PA2, lcl_PA3, lcl_PA4, lcl_psd_cumul)
    if(allocated(lcl_g_cubes)) deallocate(lcl_g_cubes)
    if(allocated(lcl_n_cubes)) deallocate(lcl_n_cubes)

    return
end subroutine pore_distribution_gpu


!--------------------------------------------------------------------------------------
! Write PSD output files (shared between CPU and GPU PSD)
!--------------------------------------------------------------------------------------
subroutine write_psd_output
    use parameters
    use lattice
    use distributions
    use results
    use fundcell, only: fundcell_getvolume, fundcell_unslant
    use vector, only: vectype
    use defaults, only: rdbl
    implicit none
    type(vectype) :: atvec1
    integer :: i, isite, nx, ny, nz
    real*8  :: x, y, z, deldis1, deldis2, deldis

    open(13, file=Trim(adjustl(property))//'_psd_cumulative.txt', status='unknown')
    open(14, file=Trim(adjustl(property))//'_psd.txt', status='unknown')

    write(13, *) '# Cumulative accessible volume distribution as a function of probe diameter'
    write(13, *) '# '
    write(13, *) '# d(probe)                Volume Fraction'

    if(psd_cumul(1) > 0.0) then
        psd_cumul = psd_cumul / psd_cumul(1)
    end if

    do i=1, nbins
        write(13,*) binsize*real(i-1)-binsize/2.0, psd_cumul(i)
    end do

    psd(1)=0.0
    psd(2)=0.0
    psd(nbins)=0.0

    do i=2, nbins-1
        deldis1=psd_cumul(i+1)
        deldis2=psd_cumul(i-1)
        deldis=deldis1-deldis2
        psd(i)=-1.0*(deldis)/(binsize*2.0)
    end do

    write(14,*) '# Derivative distribution function -dV(r)/dr (or -dV(d)/dd) vs d'
    do i=2, nbins-1
        write(14,*) binsize*real(i-1)-binsize/2.0, psd(i)
    end do

    close(13); close(14)

    open(20, file="probe_occupiable_volume.xyz", status="unknown")
    write(20,*) 0
    write(20,*)
    close(20)

    write(*,*) Trim(adjustl(property)), " cumulative PSD and differential"
    write(*,*) "PSD have been stored in files ", &
        Trim(adjustl(property))//"_psd_cumulative.txt"
    write(*,*) "and ", Trim(adjustl(property))//"_ psd.txt"
    write(*,*)
    write(*,*) "!-------------------------------------------------------!"
    write(*,*) "! Pore size distribution calculations complete          !"
    write(*,*) "!-------------------------------------------------------!"
    write(*,*)

    return
end subroutine write_psd_output
