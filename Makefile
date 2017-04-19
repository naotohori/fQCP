.SUFFIXES: .F90 .o

OPT = -O3
FC = gfortran
PYF = gnu95
PYOPT = -O3

#OPT = -O3 -xHost
#FC = ifort
#PYF = intelem
#PYOPT = -O3 -xHost

FLAG = $(OPT)

OBJS = move_coords.o qcprot.o qcprmsd.o calcrmsd.o calcrotation.o superimpose.o main.o

all:	testQCP

testQCP:	$(OBJS)
	$(FC) $(FLAG) -o testQCP $(OBJS)

clean:
	/bin/rm -f *.o testQCP
	/bin/rm -rf *.so *.so.dSYM *.pyf

##  To see the list of fortran compiler available in f2py 
#   $ f2py -c --help-fcompiler
f2py:
	f2py -m CalcRMSD -h CalcRMSD.pyf calcrmsd.F90
	f2py --fcompiler=$(PYF)  --opt=$(PYOPT) -c CalcRMSD.pyf move_coords.F90 qcprmsd.F90 calcrmsd.F90
	f2py -m CalcROT -h CalcROT.pyf calcrotation.F90
	f2py --fcompiler=$(PYF)  --opt=$(PYOPT) -c CalcROT.pyf move_coords.F90 qcprot.F90 calcrotation.F90
 
# suffix rule
.F90.o:
	$(FC) $(FLAG) -c $< 
