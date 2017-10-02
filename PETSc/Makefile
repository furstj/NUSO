CC  = mpicc
CXX = mpic++
CXXFLAGS = -std=c++11

CFLAGS  += `pkg-config --cflags PETSc`
CXXFLAGS  += `pkg-config --cflags PETSc`
LDFLAGS += `pkg-config --libs PETSc`

all: 0x_vector 1x_matrix 2x_dm 3x_snes 4x_ts

clean:
	rm -f 01_vector 02_vector 03_vector 04_vector 05_vector \
		10_matrix 11_matrix 12_matrix_metis  \
		20_dm_laplace2d \
		30_snes_trivial 31_snes_mf \
		40_ts_trivial 41_ts_implicit

0x_vector: 01_vector 02_vector 03_vector 04_vector 05_vector

1x_matrix: 10_matrix 11_matrix 12_matrix_metis 

2x_dm: 20_dm_laplace2d

3x_snes: 30_snes_trivial 31_snes_mf

4x_ts: 40_ts_trivial 41_ts_implicit

