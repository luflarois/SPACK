################################################
# Choose compilation option RUN / DEBUG
OPTION=DEBUG
COMPILE = ifort 
#
#FLAG_DEBUG= -g  #-Mbounds -Bstatic
#FLAG_RUN = -malign-double
FLAG_DEBUG= -static -O0 -CB -traceback# -Mr8 #-g 
FLAG_RUN = -static 
FLAG=$(FLAG_$(OPTION))

################################################
#
HEADER := \
	parametre \
	ficcom \
	nficfort \
#
SRC := \
        jxspack.f90 \
	commande.f90 \
	cinet.f90 \
	gestion.f90 \
	generator.f90 \
	initcinet.f90 \
	initialisation.f90 \
	lectci.f90 \
	lectcinet.f90 \
        write.f90 \
        lectciaq.f90
####################################################
#
OBJ := $(SRC:%f=%o)
#

ALL : SPACK
	@echo

SPACK : $(OBJ) Makefile
	@$(COMPILE) $(FLAG) -o $@ $(OBJ)
	mv SPACK ../Application/SPACK
	@echo xxxxxxx SPACK has been generated in ../Application.

.f.o:
	$(COMPILE) $(FLAG) -c $(SRC)






