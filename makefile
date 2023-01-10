FORTRAN=gfortran
#COPTS = -c  -O -pedantic -fbounds-check -Wall -g -fbacktrace 
#COPTS = -c  -O -pedantic -fbounds-check 
COPTS = -c  -O -pedantic -g
PROG = sabor
SUFFIXES =.f90.o
.SUFFIXES : .f90 .o
$(SUFFIXES):
	$(FORTRAN) $(COPTS) $*.f90

OBJS	      = longr.o\
		parmmage.o \
		imprime.o \
		init.o \
		intmatvec.o\
		intgradc.o\
		ajout.o\
		maillage.o\
		meshtools.o\
		fsource.o\
		kscalaire.o\
		conditioninitiale.o\
		matrixinitDDFV.o\
		scmem.o\
		ubord.o\
		tenseur.o\
		tenseurv.o\
		coefftransmv.o\
		coefftransm.o\
		newtoncDDFV.o\
		newtoneDDFV.o\
		newtonuDDFV.o\
		newtonvDDFV.o\
		plotvtkmod.o\
		gliomaddfv.o
$(PROG): $(OBJS)
	$(FORTRAN) $(OBJS) -o $(PROG)

clean:;	rm -f *.o *.mod objets/*



