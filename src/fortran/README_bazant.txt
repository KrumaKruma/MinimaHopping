f2py3 -c -llapack -lblas bazant_lib.f90 -m bazant --fcompiler=gfortran --f90flags="-fopenmp -O0" -lgomp
