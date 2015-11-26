all: simple2 simple unsteady

CFLAGS=-O 

include ${PETSC_DIR}/conf/variables
include $(PETSC_DIR)/conf/rules

INCLUDES+=-I$(PETSC_INCLUDE)  
LIBS+=$(PETSC_LIB) -lX11 -llapack -lblas

simple: simple.o chkopts 
	mpic++  -o $@ $<  $(LIBS)

simple2: simple2.o chkopts
	mpic++  -o $@ $<  $(LIBS)

unsteady: unsteady.o chkopts
	mpic++ -o $@ $< $(LIBS)

