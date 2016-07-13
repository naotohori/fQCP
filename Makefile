.SUFFIXES: .F90 .o

OPT = -O3
FC = gfortran

FLAG = $(OPT)

OBJS = qcprot.o main.o

all:	testQCP

testQCP:	$(OBJS)
	$(FC) $(FLAG) -o testQCP $(OBJS)

clean:
	/bin/rm -f *.o testQCP

# suffix rule
.F90.o:
	$(FC) $(FLAG) -c $< 
