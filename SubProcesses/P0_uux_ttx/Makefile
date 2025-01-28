### Compilers r###
CC=g++## for C++
FF=gfortran## for Fortran

###  Directories ###
SRCDIR=${PWD}/src
INCDIR=${PWD}/inc
LHAPDFDIR=$//home/matteo/LHAPDF/include
LHAPDFDIR2=$//home/matteo/LHAPDF/share/LHAPDF
LIBDIR=${PWD}/lib
BLDDIR=${PWD}/bld
GSLDIR=$//home/matteo/gsl-2.8
BOOSTDIR=$//usr/include/boost
EIGENDIR=$//home/matteo/eigen-3.4.0
MADLOOPDIR=${PWD}
MADLOOPRES=${PWD}/MadLoop5_resources
MPLIBDIR=${PWD}/../../lib
CUTOOLSDIR=${PWD}/../../Source/CutTools/includects
REGIDIR=${PWD}/../../Source/IREGI/src
PSDIR=${PWD}/../PhaseSpaceGenerator

### Flags ###
CFLAG=-O3 -Wall -I${INCDIR} -I${LHAPDFDIR} -I${LHAPDFDIR2} -I${MPLIBDIR} -I${MADLOOPDIR} -I${GSLDIR} -I${BOOSTDIR} -I${EIGENDIR} -I${MADLOOPRES}
FFLAG=-O3 -I${INCDIR} -I${PWD}/RonbergFort -I${PWD}/HSVirtualH1gg -I${PSDIR}

### Paths ###
VPATH=${SRCDIR}

### Source files ###
CFILES=$(wildcard ${SRCDIR}/*.cpp)
FFILES=$(wildcard ${SRCDIR}/*.f)

### Object files ###
COBJS=$(subst .cpp,.o,$(subst ${SRCDIR},${BLDDIR},${CFILES}))
FOBJS=$(subst .f,.o,$(subst ${SRCDIR},${BLDDIR},${FFILES}))

### Libraries ###
LIB=${LIBDIR}/libresum.a
GSLIB=-L/usr/lib/x86_64-linux-gnu -lgsl -lgslcblas
STLIB= -L/usr/lib/gcc/x86_64-linux-gnu/13 -lm -lgfortran 
LHAPDFLIB= -L/home/matteo/LHAPDF/lib -lLHAPDF
MADLOOPLIB=-L$(MADLOOPDIR) -lMadLoop
MPLIB= -L${MPLIBDIR} -ldhelas -lmodel
LCTLIB= -L${CUTOOLSDIR} -lcts
REGILIB= -L${REGIDIR} -liregi
INTLIB= -L${PWD}/RonbergFort -lnin
PSLIB= -L${PSDIR} -lgenps
#INTLIB2= -L${PWD}/HSVirtualH1gg -lHFgg2ttx
### Commands ###
all: RUN
lib: ${LIB}

RUN: main.cpp ${LIB}
	${CC} ${CFLAG} -o $@ main.cpp ${LIB} ${STLIB} ${GSLIB} ${LHAPDFLIB} ${INTLIB} ${PSLIB} $(MADLOOPLIB) -ldl ${MPLIB} ${LCTLIB} ${REGILIB} -lquadmath -lgfortran 


${LIB}:	${COBJS} ${FOBJS}
	ar -ruc $@ ${BLDDIR}/*.o

${BLDDIR}%.o: ${SRCDIR}%.f
	cd ${BLDDIR}; ${FF} -c ${FFLAG} $?;

${BLDDIR}%.o: ${SRCDIR}%.cpp
	cd ${BLDDIR}; ${CC} -c ${CFLAG} $<;
	
### Cleaning ###
clean:
	rm -f RUN ${LIBDIR}/*.a ${BLDDIR}/*.o *.log output/*; clear;

