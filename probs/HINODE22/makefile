PROGDIR = .

################### yellowstone################################
#FORT = mpif90
#CPP = cpp
#CPPFLAG = -P -traditional -DRTPGR -DJCROSSB -DFQTANH -DTHCONDUCT -DHEAT_RADLOSS
#LDOPTS =
#FFLAGS = -O3 -r8 -fp-model strict
###############################################################

################### cheyenne ################################
FORT = mpif90
CPP = cpp
CPPFLAG = -P -traditional -DRTPGR -DJCROSSB -DFQTANH -DTHCONDUCT -DHEAT_RADLOSS
LDOPTS =
FFLAGS = -O3 -r8 -fp-model strict -march=corei7 -axAVX
###############################################################

################### topaz #####################################
#FORT = ifort
#CPP = cpp
#CPPFLAG = -P -traditional -DRTPGR -DJCROSSB -DFQTANH -DTHCONDUCT -DHEAT_RADLOSS
#LDOPTS = -lmpi
#FFLAGS = -O3 -r8 -fp-model strict
###############################################################

################### pleiades ##################################
#FORT = ifort
#CPP = cpp
#CPPFLAG = -P -traditional -DRTPGR -DJCROSSB -DFQTANH -DTHCONDUCT -DHEAT_RADLOSS
#LDOPTS = -lmpi
###FFLAGS = -O3 -r8 -axAVX -xSSE4.1 -fp-model strict -vec-report1
###FFLAGS = -O3 -r8 -axAVX -xSSE4.1 -fp-model strict
###FFLAGS = -O3 -r8 -axCORE-AVX2 -xSSE4.2 -fp-model strict
###FFLAGS = -O3 -r8 -axCORE-AVX512 -xSSE4.2 -fp-model strict
#FFLAGS = -O3 -r8 -axCORE-AVX512,CORE-AVX2,AVX -xSSE4.2 -fp-model strict
##############################################################

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
