CC  = mpic++
CXX = mpic++
CXXFLAGS = -std=c++11 -O3

CFLAGS  += `pkg-config --cflags PETSc`
CXXFLAGS  += `pkg-config --cflags PETSc`
LDFLAGS += `pkg-config --libs PETSc`

all: 01_scalar_upwind 02_scalar_upwind 03_scalar_upwind_ts

01_scalar_upwind: 01_scalar_upwind.o CartesianGrid.o Grid.o 

02_scalar_upwind: 02_scalar_upwind.o CartesianGrid.o Grid.o

03_scalar_upwind_ts: 03_scalar_upwind_ts.o CartesianGrid.o Grid.o

CartesianGrid.o: CartesianGrid.cpp CartesianGrid.hpp Grid.hpp
