BIN = ../bin
GPROG =	$(BIN)/FULL_MONTE

OBJS1 = commons.o
OBJS2 =	accumul.o initialise.o keywords.o main.o move.o move_mc.o umbrella.o readinput.o \
        output.o potential.o orderparam.o free_energy_mc.o viewconfig.o \
        cluster_move.o cell_list.o rotations_module.o stack.o unitcells.o glj.o \
		kihara.o pglj.o kf.o cpp.o etp.o dmbl_glj.o kf_rec.o hdmbl.o

OBJS = ${OBJS1} ${OBJS2}
LDFLAGS = -L.
DEFS =
# CPP = /usr/bin/cpp
CPP = /lib/cpp
CPFLAGS = -traditional -P

###################################### COMPILERS AND COMPILER FLAGS ######################################
FC = gfortran
FFLAGS= -ffixed-line-length-132 -g -fbounds-check -Wuninitialized -O -fimplicit-none -fno-automatic 
# -fallow-argument-mismatch
SWITCH=gfortran
LIBS  = drand48.c

###################################### RULES AND TARGETS ######################################
.SUFFIXES:
.SUFFIXES: .o .f .f90

.f90.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<
.f.o:
	$(FC) $(FFLAGS) ${SEARCH_PATH} -c $<

default: $(GPROG)

$(GPROG): $(CHDUM) $(AMHDUM) $(AMB9DUM) $(OBJS) $(EXTRAS) | $(BIN)
	$(FC) $(FFLAGS) ${SEARCH_PATH} -o $@ $(EXTRAS) $(OBJS) $(CHDUM) $(AMHDUM) $(AMB9DUM) $(LDFLAGS) $(LIBS)

$(BIN):
	mkdir $@

clean:
	rm -f ${GPROG} $(OBJS) *.mod $(EXTRAS) *.lst *.o nag.f90

rebuild:
	make clean
	make

###################################### DEPENDENCIES ######################################
${OBJS2}: ${OBJS1}
accumul.o:        	commons.o orderparam.o
initialise.o:     	commons.o cell_list.o rotations_module.o orderparam.o cluster_move.o
unitcells.o:		commons.o
keywords.o:       	commons.o
main.o:           	commons.o
move.o:           	commons.o free_energy_mc.o
move_mc.o:        	commons.o orderparam.o cell_list.o rotations_module.o cluster_move.o
cluster_move.o:		commons.o cell_list.o rotations_module.o
umbrella.o:			commons.o orderparam.o cell_list.o
output.o:         	commons.o 
potential.o:      	commons.o
free_energy_mc.o: 	commons.o cell_list.o rotations_module.o cluster_move.o
readinput.o:      	commons.o
glj.o:            	commons.o
kihara.o:			commons.o rotations_module.o
kf.o:				commons.o rotations_module.o
pglj.o:				commons.o rotations_module.o
cpp.o:				commons.o rotations_module.o
etp.o:				commons.o rotations_module.o
dmbl_glj.o:	  		commons.o rotations_module.o
kf_rec.o:	  		commons.o rotations_module.o
hdmbl.o:           	commons.o
orderparam.o:     	commons.o cell_list.o stack.o cluster_move.o
viewconfig.o:     	commons.o
cell_list.o:
rotations_module.o:
stack.o:
