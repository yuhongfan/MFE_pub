PROGDIR = .

################### mac with gfortran and openmpi ######################
FORT = mpif90
CPP = cpp
CPPFLAG = -P -traditional -DJCROSSB
LDOPTS =
FFLAGS = -O3 -fallow-argument-mismatch -fdefault-real-8 -ftree-vectorize
###############################################################

################### cheyenne ################################
#FORT = mpif90
#CPP = cpp
#CPPFLAG = -P -traditional -DJCROSSB
#LDOPTS =
#FFLAGS = -O3 -r8 -fp-model strict -march=corei7 -axAVX
###############################################################

POSTCPP = $(PROGDIR)/postcpp

.SUFFIXES: .o .f90 .F

MODULES = ModPar.o ModPhysics.o ModGrid.o ModSundry.o ModControl.o \
ModBoundary.o ModField.o ModDel.o ModWork.o ModConductivity.o \
ModIO.o ModInterp.o ModGetQtys.o \
ModUserSetup.o ModBval.o ModInitialization.o \
ModRHS.o ModNonIdealRHS.o ModIteration_sts.o

OBJ = main.o

$(OBJ): $(MODULES)

.F.o:
	$(CPP) $(CPPFLAG) $*.F $(POSTCPP)/$*.f90
	$(FORT) $(FFLAGS) $(POSTCPP)/$*.f90 -c

xmfe: $(OBJ)
	$(FORT) $(FFLAGS) -o xmfe $(OBJ) $(MODULES) $(LDOPTS)

clean:
	rm *.o *.mod xmfe postcpp/*.f90
