/*
  Program demonstrujici precislovani a distribuci matice

  Matice o rozmerech n x n (n sude) ma takovou strukturu, ze liche prvky tvori kliku, 
  sude prvky tvori druhou kliku a navice je vrchol 0 spojen s vrcholem 1.
  
  Program je urcen pro 2 procesory, mel by ale fungovat pro libovolny pocet.


  Vhodne rozdeleni na 2 procesory bude takove, ze sude prvky budou na 
  procesoru 0 a liche prvky na procesoru 1.

  tj. vektor 
  0 1 2 3 4 5 6 7 8 9

  bude prerovnan a rozdelen takto:
  0 2 4 6 8 | 1 3 5 7 9

  Chceme-li tedy zapsat prvek na pozici 3 v puvodnim (aplikacnim) cislovani, 
  musime ho ulozit do PETSc vektoru na pozici 6.
 */

#include <iostream>
#include <utility>

#include <petscvec.h>
#include <petscmat.h>
#include <petscao.h>
#include <petscksp.h>


// Vytvoreni matice s danou permutaci
// na diagonale je i+n, mimo diagonalu je -1
Mat CreateMatrix(int n, AO ao, int localSize) {
  Mat A;
  MatCreate(PETSC_COMM_WORLD, &A);
  MatSetSizes(A, localSize, localSize, n, n);
    
  MatSetFromOptions(A);
  MatSetUp(A);

  PetscInt istart, iend;
  MatGetOwnershipRange(A, &istart, &iend);

  int ia, ja;
  for (int i=istart; i<iend; i++) {
    ia = i;
    if (ao != PETSC_NULL) AOPetscToApplication(ao, 1, &ia);

    for (int j=i%2; j<n; j+=2) {
      ja = j;
      if (ao != PETSC_NULL) AOPetscToApplication(ao, 1, &ja);
      
      if (i==j) 
	MatSetValue(A, ia, ja, n+i, INSERT_VALUES);
      else
	MatSetValue(A, ia, ja, -1, INSERT_VALUES);
    }
  }
  ia = 0;
  ja = 1;
  if (ao != PETSC_NULL) AOPetscToApplication(ao, 1, &ia);
  if (ao != PETSC_NULL) AOPetscToApplication(ao, 1, &ja);

  MatSetValue(A, ia, ja, -1, INSERT_VALUES);
  MatSetValue(A, ja, ia, -1, INSERT_VALUES);
  
  MatAssemblyBegin(A, MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(A, MAT_FINAL_ASSEMBLY);

  return A;
}

/* Vytvoreni nenulove struktury matice. 
   PETSc podporuje prime vytvoreni pomoci MatCreateMPIAdj. 
   Zde je ale pouzito vytvoreni pomoci MatConvert
*/
Mat CreateAdjacencyMatrix(int n) {
  Mat A = CreateMatrix(n, PETSC_NULL, PETSC_DECIDE);

  Mat Adj;
  MatConvert(A, MATMPIADJ, MAT_INITIAL_MATRIX, &Adj);
  
  MatDestroy(&A);

  return Adj;
}


int main(int argc,char **args)
{
  const int n = 10; // Velikost matice
  int myCpu;
  int mySize;

  // Inicializace 
  PetscInitialize( &argc , &args , (char *)0 , 0 );
  MPI_Comm_rank(PETSC_COMM_WORLD, &myCpu);
  MPI_Comm_size(PETSC_COMM_WORLD, &mySize);

  // Vytvoreni nenulove struktury matice 
  Mat Adj = CreateAdjacencyMatrix(n);
  PetscPrintf(PETSC_COMM_WORLD, "Adjacency matrix;\n");
  MatView(Adj, PETSC_VIEWER_STDOUT_WORLD);

  // Vypocet rozdeleni matice na procesory 
  MatPartitioning part;
  MatPartitioningCreate(PETSC_COMM_WORLD, &part);
  MatPartitioningSetAdjacency(part, Adj);
  MatPartitioningSetFromOptions(part);
  
  IS partCpu;
  MatPartitioningApply(part, &partCpu);
  PetscPrintf(PETSC_COMM_WORLD, "\nPartitioning:\n");
  ISView(partCpu, PETSC_VIEWER_STDOUT_WORLD);
  
  // Precislovani vrcholu grafu 
  IS permutation;
  ISPartitioningToNumbering(partCpu, &permutation);
  AO ao;
  AOCreateBasicIS(permutation, PETSC_NULL, &ao);
  PetscPrintf(PETSC_COMM_WORLD, "\nRenumbering:\n");
  AOView(ao, PETSC_VIEWER_STDOUT_WORLD);

  // Kolik prvku by melo padnout na dany procesor
  int localSize[mySize];
  ISPartitioningCount(partCpu, mySize, localSize);
  if (myCpu==0) {
    std::cout << "Local sizes:" << std::endl;
    for (int i=0; i<mySize; i++) 
      std::cout << localSize[i] << std::endl;
  }

  // Vytvoreni prerovnane matice a prave strany
  Mat A = CreateMatrix(n, ao, localSize[myCpu]);
  PetscPrintf(PETSC_COMM_WORLD, "\nMatrix:\n");
  MatView(A, PETSC_VIEWER_STDOUT_WORLD);

  Vec b;
  VecCreateMPI(PETSC_COMM_WORLD, PETSC_DECIDE, n, &b);
  if (myCpu==0) {
    int ib[n];
    double bb[n];
    for (int i=0; i<n; i++) { ib[i]=i; bb[i]=i; } 
    AOPetscToApplication(ao, n, ib);
    VecSetValues(b, n, ib, bb, INSERT_VALUES);
  };
  VecAssemblyBegin(b); VecAssemblyEnd(b);
  PetscPrintf(PETSC_COMM_WORLD, "\nrhs:\n");
  VecView(b, PETSC_VIEWER_STDOUT_WORLD);

  // Vyreseni soustavy
  Vec x;  
  VecDuplicate(b, &x);

  KSP solver;
  KSPCreate(PETSC_COMM_WORLD, &solver);
  KSPSetOperators(solver, A, A);
  KSPSetFromOptions(solver);
  KSPSetUp(solver);
  KSPSolve(solver, b, x);

  PetscPrintf(PETSC_COMM_WORLD, "\nsolution (permuted):\n");
  VecView(x, PETSC_VIEWER_STDOUT_WORLD);

  // Zpetna permutace vysledku vysledku
  Vec xp; VecDuplicate(x,&xp);
  VecScatter scat;
  VecScatterCreate(x, permutation, xp, PETSC_NULL, &scat);
  VecScatterBegin(scat, x, xp, INSERT_VALUES, SCATTER_FORWARD);
  VecScatterEnd(scat, x, xp, INSERT_VALUES, SCATTER_FORWARD);
  
  PetscPrintf(PETSC_COMM_WORLD, "\nsolution :\n");
  VecView(xp, PETSC_VIEWER_STDOUT_WORLD);

  VecScatterDestroy(&scat);
  KSPDestroy(&solver);
  VecDestroy(&xp);
  VecDestroy(&x);
  VecDestroy(&b);
  MatDestroy(&A);
  AODestroy(&ao);
  ISDestroy(&permutation);
  ISDestroy(&partCpu);
  MatPartitioningDestroy(&part);
  MatDestroy(&Adj);
  PetscFinalize();
  
  return 0;
}
