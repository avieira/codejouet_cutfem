default : ALL

INCLUDE_DIR = -I./GraphBLAS/Include -Iinclude
LIB_DIR = -L./GraphBLAS/build
DEBUG = -g
EXEC = main.exe
OBJ = obj
CC = cc
CF = gfortran
OPTS_C = -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-function -std=c11 -lm -lgraphblas 
OPTS_F = -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-function -std=f95 -lm #-fno-underscoring 
#CC = icx
#CF = ifx
#OPTS_C = -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-command-line-argument -std=c11 -lm 
##OPTS_F = -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-function -Wno-unused-command-line-argument -stand f95 -lm 
LINK_OPTS = -Wl,-rpath,/home/castor/Documents/travaux/bridge_fortran_sparsesuite/GraphBLAS/build ./GraphBLAS/build/libgraphblas.so.10.0.2

SOURCES = $(wildcard *.c) 
SOURCES_F = $(wildcard *.F90)
OBJECTS = $(patsubst %.c, $(OBJ)/%.o, $(SOURCES))
OBJECTS_F = $(patsubst %.F90, $(OBJ)/%.o, $(SOURCES_F))

ALL : $(EXEC)

$(OBJ)/%.o: %.c
	$(CC) $(DEBUG) $(INCLUDE_DIR) $(LIB_DIR) $(OPTS_C) -c $< -o $@

$(OBJ)/%.o: %.F90
	$(CF) $(DEBUG) $(OPTS_F) -c $< -o $@

$(OBJ)/ale_connectivity_mod.o: ale_connectivity_mod.F
	$(CF) $(DEBUG) -Wall -Wextra -Wno-unused-but-set-variable -Wno-unused-function -lm -c $< -o $@
                                              
$(OBJ)/grid_struct_multicutcell_mod.o : $(OBJ)/polygon_mod.o
$(OBJ)/multicutcell_solver_mod.o : $(OBJ)/polygon_mod.o $(OBJ)/ale_connectivity_mod.o $(OBJ)/grid_struct_multicutcell_mod.o $(OBJ)/integer_linked_list.o $(OBJ)/riemann_solver_mod.o
$(OBJ)/bridge_fortran_c.o : $(OBJ)/multicutcell_solver_mod.o


$(EXEC) : $(OBJECTS) $(OBJECTS_F)
	$(CF) $(DEBUG) $(INCLUDE_DIR) $(LIB_DIR) $^ $(OPTS_F) $(LINK_OPTS) -o $(EXEC)

clean:
	rm *.mod $(OBJ)/* $(EXEC)