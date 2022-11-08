SciCORE:
ml purge
ml PyTorch
f2py3 -c -lm -lopenblas -lpthread bazant_lib.f90 -m bazant --fcompiler=gfortran --f90flags="-fopenmp"

gfortran:
f2py3 -c -llapack -lblas bazant_lib.f90 -m bazant --fcompiler=gfortran --f90flags="-fopenmp -O0" -lgomp
