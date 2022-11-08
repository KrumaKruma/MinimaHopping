SciCORE:
ml purge
ml PyTorch
f2py3 -c -lm -lopenblas -lpthread symmetry_penalty_v3.f90 -m sympenlib --fcompiler=gfortran --f90flags="-fopenmp"


gfortran:
gfortran symmetry_penalty_v3.f90  -llapack -lblas -fopenmp

f2py3 -c -llapack -lblas symmetry_penalty_v3.f90 -m sympenlib --fcompiler=gfortran --f90flags="-fopenmp -O0" -lgomp



Wichtig!
https://github.com/numpy/numpy/issues/11908
Fopenmp:
I usually do this for gfortran so fopenmp in f90 flag AND -lgomp for GCC compiler
f2py -c file.f90 -m file_module --f90flags='-fopenmp' -lgomp
and for ifort
f2py -c file.f90 -m file_module --fcompiler=intelem --f90flags='-qopenmp' -liomp5

f2py3 -c -L/usr/lib/x86_64-linux-gnu/ -llapack -lblas symmetry_penalty_v3.f90 -m sympenlib --fcompiler=gfortran --f90flags=" -O0" -lgomp
