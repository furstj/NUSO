all: simple2 simple unsteady

CC  = mpicc
CXX = mpic++
CXXFLAGS = -std=c++11

CFLAGS  += `pkg-config --cflags PETSc` 
CXXFLAGS  += `pkg-config --cflags PETSc` 
LDFLAGS += `pkg-config --libs PETSc`

