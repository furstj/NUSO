CC  = mpic++
CXX = mpic++
CXXFLAGS = -std=c++11 -O3

CFLAGS  += `pkg-config --cflags PETSc`
CXXFLAGS  += `pkg-config --cflags PETSc`
LDFLAGS += `pkg-config --libs PETSc`

all: ex 0x_dmplex

ex: ex1 ex11

0x_dmplex: 01_dmplex_basic 02_dmplex_gmsh 03_dmplex_gmsh


