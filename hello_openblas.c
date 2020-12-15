#include <cblas.h>
#include <stdio.h>

void main()
{
  int i = 0;
  double A[6] = {1.0, 2.0, 1.0, -3.0, 4.0, -1.0}; 
  double B[6] = {1.0, 2.0, 1.0, -3.0, 4.0, -1.0};
  double C[9] = {.5, .5, .5, .5, .5, .5, .5, .5, .5};

  // C = a*A*B(T) + b*C

  cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, 3, 3, 2, 1, A, 3, B, 3, 0, C, 3);

  double D[9] = {.5, .5, .5, .5, .5, .5, .5, .5, .5};
  cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 3, 3, 2, 1, A, 2, B, 2, 0, D, 3);
  
  printf("\nC = ");
  for (i = 0; i < 9; i++)
    printf("%lf ", C[i]);
  printf("\n");

  printf("\nD = ");
  for (i = 0; i < 9; i++)
    printf("%lf ", D[i]);
  printf("\n");
}