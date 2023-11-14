FORT = ifort
LDOPTS = -L /glade/u/home/yfan/lib -lfishpack
FLAGS = -mcmodel=medium -O3

.SUFFIXES: .o .f90 .f

OBJ = problemsize.o blk.o grid.o input_namelists.o blkdata.o \
initialization.o input.o main_windinit_pdfi_chromfix.o lint.o

.f90.o:
	$(FORT) $(FLAGS) $*.f90 -c
.f.o:
	$(FORT) $(FLAGS) $*.f -c

xinit: $(OBJ)
	$(FORT) $(FLAGS) -o xinit -shared-intel $(OBJ) $(LDOPTS)

clean:
	rm *.o *.mod xinit
