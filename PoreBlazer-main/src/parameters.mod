V35 :0x24 parameters
14 poreblazer.f90 S624 0
05/21/2026  09:25:42
use matrix private
use vector private
use fundcell private
enduse
D 388 26 2051 432 2024 7
D 477 20 500
S 624 24 0 0 0 9 1 0 4986 5 8000 A 0 0 0 0 B 0 10 0 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 10 0 0 0 0 0 0 parameters
S 626 23 0 0 0 9 2024 624 5006 4 0 A 0 0 0 0 B 400000 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 fundamental_cell
S 627 23 0 0 0 9 2125 624 5023 4 0 A 0 0 0 0 B 400000 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 fundcell_getvolume
R 1336 26 10 vector =
R 1337 26 11 vector *
R 1338 26 12 vector +
R 1339 26 13 vector -
R 1340 26 14 vector /
R 1769 26 10 matrix =
R 1770 26 11 matrix +
R 1771 26 12 matrix -
R 1772 26 13 matrix *
R 1773 26 14 matrix /
R 2008 26 13 fundcell =
R 2009 26 14 fundcell +
R 2010 26 15 fundcell -
R 2011 26 16 fundcell *
R 2012 26 17 fundcell /
R 2024 25 29 fundcell fundamental_cell
R 2051 5 56 fundcell orthoflag fundamental_cell
R 2052 5 57 fundcell ell fundamental_cell
R 2053 5 58 fundcell half_ell fundamental_cell
R 2054 5 59 fundcell ell_inv fundamental_cell
R 2055 5 60 fundcell anglebc fundamental_cell
R 2056 5 61 fundcell angleac fundamental_cell
R 2057 5 62 fundcell angleab fundamental_cell
R 2058 5 63 fundcell origin fundamental_cell
R 2059 5 64 fundcell eff fundamental_cell
R 2060 5 65 fundcell volume fundamental_cell
R 2061 5 66 fundcell slantmatrix fundamental_cell
R 2062 5 67 fundcell unslantmatrix fundamental_cell
R 2063 5 68 fundcell lvec fundamental_cell
R 2064 5 69 fundcell width fundamental_cell
R 2065 5 70 fundcell minwidth fundamental_cell
R 2066 5 71 fundcell mx fundamental_cell
R 2067 5 72 fundcell my fundamental_cell
R 2068 5 73 fundcell mz fundamental_cell
R 2125 14 130 fundcell fundcell_getvolume
S 2285 6 4 0 0 388 2291 624 10636 4 0 A 0 0 0 0 B 0 12 0 0 0 0 0 0 0 0 0 0 2296 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 fcell
S 2286 6 4 0 0 10 2287 624 11032 4 0 A 0 0 0 0 B 0 13 0 0 0 0 0 0 0 0 0 0 2297 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 hicut
S 2287 6 4 0 0 10 2288 624 11038 4 0 A 0 0 0 0 B 0 13 0 0 0 8 0 0 0 0 0 0 2297 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 hicut2
S 2288 6 4 0 0 10 2289 624 11045 4 0 A 0 0 0 0 B 0 14 0 0 0 16 0 0 0 0 0 0 2297 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 coeff_surface
S 2289 6 4 0 0 10 2290 624 11059 4 0 A 0 0 0 0 B 0 14 0 0 0 24 0 0 0 0 0 0 2297 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 coeff_surface2
S 2290 6 4 0 0 10 1 624 11074 4 0 A 0 0 0 0 B 0 15 0 0 0 32 0 0 0 0 0 0 2297 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 temp
S 2291 6 4 0 0 6 2292 624 11079 4 0 A 0 0 0 0 B 0 16 0 0 0 432 0 0 0 0 0 0 2296 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 nsample
S 2292 6 4 0 0 6 2293 624 8997 4 0 A 0 0 0 0 B 0 0 0 0 0 436 0 0 0 0 0 0 2296 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 iseed
S 2293 6 4 0 0 6 1 624 11087 4 0 A 0 0 0 0 B 0 18 0 0 0 440 0 0 0 0 0 0 2296 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 vis_option
S 2294 3 0 0 0 6 0 1 0 0 0 A 0 0 0 0 B 0 0 0 0 0 0 0 0 0 20 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 6
S 2295 6 4 0 0 477 1 624 11098 4 0 A 0 0 0 0 B 0 19 0 0 0 0 0 0 0 0 0 0 2298 0 0 0 0 0 0 0 0 0 0 624 0 0 0 0 property
S 2296 11 0 0 0 9 2072 624 11107 40800000 805000 A 0 0 0 0 B 0 20 0 0 0 444 0 0 2285 2293 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _parameters$0
S 2297 11 0 0 0 9 2296 624 11121 40800000 805000 A 0 0 0 0 B 0 20 0 0 0 40 0 0 2286 2290 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _parameters$2
S 2298 11 0 0 0 9 2297 624 11135 40800000 805000 A 0 0 0 0 B 0 20 0 0 0 20 0 0 2295 2295 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 _parameters$1
A 500 2 0 0 0 6 2294 0 0 0 500 0 0 0 0 0 0 0 0 0 0 0
Z
Z
