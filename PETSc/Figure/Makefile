CC  = mpicc
CXX = mpic++
CXXFLAGS = -std=c++11 -g

CFLAGS  += `pkg-config --cflags PETSc`
CXXFLAGS  += `pkg-config --cflags PETSc`
LDFLAGS += `pkg-config --libs PETSc`

all: blur blur-ts

clean:
	rm -f blur blur-ts

