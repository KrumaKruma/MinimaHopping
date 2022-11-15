SciCORE:
ml Python/intel2019
f2py -c --fcompiler=intelem --f90flags="-qopenmp" -lmkl symmetry_penalty_v3.f90 -m sympenlib

ml purge
ml PyTorch
f2py3 -c -lscalapack -lopenblas --f90flags="-fopenmp" -lgomp symmetry_penalty_v3.f90 -m sympenlib
f2py3 -c -lscalapack -lopenblas --f90flags="-fopenmp -Ofast -ffast-math" -lgomp symmetry_penalty_v3.f90 -m sympenlib
f2py3 -c -lscalapack -lopenblas --f90flags="-fopenmp -Ofast -march=znver1 -mtune=znver1 -mfma -mavx2 -fomit-frame-pointer -ffast-math" -lgomp symmetry_penalty_v3.f90 -m sympenlib

BLISS:::
ml foss/2021a
f2py3 -c --f90flags="-fopenmp -fPIC" -lgomp -lpthread symmetry_penalty_v3.f90 /kernph/gubmor00/lapack-testing/flame/amd-libflame/lib/LP64/libflame.a /kernph/gubmor00/lapack-testing/blis/amd-blis/lib/LP64/libblis.a -m sympenlib



gfortran:
gfortran symmetry_penalty_v3.f90  -llapack -lblas -fopenmp

f2py3 -c -llapack -lblas symmetry_penalty_v3.f90 -m sympenlib --fcompiler=gfortran --f90flags="-fopenmp -O0" -lgomp




Timings reiner Fortran Code:
gfortran with -lopenblas -Ofast -march=native 
same as ifort -mkl -Ofast

Andris FLAGS:
Hi everyone,

I was doing some tests on the new AMD nodes with my code using openmp.
I could confirm what Andris had already discovered.
Using the following openmp settings improved the performance by 30% or so.
export OMP_PROC_BIND=TRUE
export OMP_PLACES=threads

But I also found that gfortran generated much faster machine code than ifort.
Using the latest version available on scicore (10.3.0) brought an additional speedup but also older versions were faster than ifort.

In my tests with 20 openmp threads I got the following timings for my code.
ifort: 3:40 min
gfortran (v9): 2:41 min
gfortran (v10): 2:11 min

My flags were the following:
gfortran: -Ofast -march=znver1 -mtune=znver1 -mfma -mavx2 -fomit-frame-pointer -ffast-math
ifort: -march=core-avx2 -align array64byte -fma -fomit-frame-pointer -mkl -Ofast -fast

Has anyone of you had similar experiences with gfortran being faster or maybe found some flags that work better with ifort?

Best wishes
Jonas




Wichtig!
https://github.com/numpy/numpy/issues/11908
Fopenmp:
I usually do this for gfortran so fopenmp in f90 flag AND -lgomp for GCC compiler
f2py -c file.f90 -m file_module --f90flags='-fopenmp' -lgomp
and for ifort
f2py -c file.f90 -m file_module --fcompiler=intelem --f90flags='-qopenmp' -liomp5

f2py3 -c -L/usr/lib/x86_64-linux-gnu/ -llapack -lblas symmetry_penalty_v3.f90 -m sympenlib --fcompiler=gfortran --f90flags=" -O0" -lgomp
