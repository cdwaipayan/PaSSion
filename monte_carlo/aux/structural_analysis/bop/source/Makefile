
GPROG =	../bin/BOP

OBJS1 = commons.o
OBJS2 =	orderparam.o cell_list.o stack.o rotations_module.o main.o

OBJS = ${OBJS1} ${OBJS2}
LDFLAGS = -L.
DEFS =
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P

###################################### COMPILERS AND COMPILER FLAGS ######################################
FC = gfortran
FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -fimplicit-none -fno-automatic
SWITCH=gfortran
LIBS  =

###################################### RULES AND TARGETS ######################################
.SUFFIXES:
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<

default: $(GPROG)

$(GPROG): $(CHDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(EXTRAS)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS) $(OBJS) $(CHDUM) $(AMHDUM) $(AMB9DUM) $(LDFLAGS) $(LIBS)

clean:
	rm -f ${GPROG} $(OBJS) *.mod $(EXTRAS) *.lst *.o nag.f90

rebuild:
	make clean
	make

###################################### DEPENDENCIES ######################################
${OBJS2}: ${OBJS1}

main:				commons.o cell_list.o orderparam.o rotations_module.o
orderparam.o:     	commons.o cell_list.o stack.o
cell_list.o:		commons.o
rotations_module.o: commons.o
stack.o:
