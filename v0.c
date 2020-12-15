#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <cblas.h>

// // Definition of the kNN result struct
// typedef struct knnresult
// {
//     int *nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
//     double *ndist; //!< Distance of nearest neighbors          [m-by-k]
//     int m;         //!< Number of query points                 [scalar]
//     int k;         //!< Number of nearest neighbors            [scalar]
// } knnresult;

// //! Compute k nearest neighbors of each point in X [n-by-d]
// /*!

//   \param  X      Corpus data points              [n-by-d]
//   \param  Y      Query data points               [m-by-d]
//   \param  n      Number of corpus points         [scalar]
//   \param  m      Number of query points          [scalar]
//   \param  d      Number of dimensions            [scalar]
//   \param  k      Number of neighbors             [scalar]

//   \return  The kNN result
// */
// knnresult kNN(double *X, double *Y, int n, int m, int d, int k)
// {
// }

/**
 * Function that calculates an m×n Euclidean distance matrix D
 * D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).');
 * X and Y lists with CblasColMajor format
 * X[n*d] and Y[m*d]
 */
void distanceMatrix(double *X, double *Y, int n, int m, int d)
{
    /* Allocate memory for sumX and sumY lists*/
    double *sumX = (double *)malloc(n * sizeof(double));
    double *sumY = (double *)malloc(m * sizeof(double));

    /* Check if the memory has been successfully allocated */
    if (sumX == NULL || sumY == NULL)
    {
        printf("Not enough memory to allocate for sumX and sumY ");
        free(sumX);
        free(sumY);
        return;
    }

    /* Calculate sumX list*/
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < d; j++)
        {
            sumX[i] += pow(X[i + j * n], 2);
        }
    }

    /* Calculate sumY list*/
    for (int i = 0; i < m; i++)
    {
        for (int j = 0; j < d; j++)
        {
            sumY[i] += pow(Y[i + j * m], 2);
        }
    }

    printf("\nSumX = ");
    for (int i = 0; i < n; i++)
    {
        printf(" %lf ", sumX[i]);
    }

    printf("\nSumY = ");
    for (int i = 0; i < m; i++)
    {
        printf(" %lf ", sumY[i]);
    }
    
    /* Allocate memory for C and D lists*/
    double alpha = 1.0, beta = 0;
    double *C = (double *)malloc(n * m * sizeof(double));
    double *D = (double *)malloc(n * m * sizeof(double));

    /* Check if the memory has been successfully allocated */
    if (C == NULL || D == NULL)
    {
        printf("Not enough memory to allocate for C and D ");
        free(C);
        free(D);
        return;
    }

    /* Calculate list C = X*Y */
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, n, m, d, alpha, X, n, Y, m, beta, C, n);

    printf("\nC = ");
    for (int i = 0; i < n * m; i++)
    {
        printf("%lf ", C[i]);   
    }

    /* Calculate list D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).'); */
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            D[i + j * n] = sqrt(sumX[i] - 2 * C[i + j * n] + sumY[j]);
        }
    }   

    printf("\n\nD = ");
    for (int i = 0; i < n * m; i++)
    {
        printf("%lf ", D[i]);
    }
    printf("\n");

    /* Deallocate used memory */
    free(sumX);
    free(sumY);
    free(C);
    free(D);
}

int main()
{   
    int n = 3, m = 2, d = 2, k = 2;

    /* Allocate memory for X and Y lists*/
    double *X = (double *)malloc(n * d * sizeof(double));
    double *Y = (double *)malloc(m * d * sizeof(double));

    /* Check if the memory has been successfully allocated */
    if (X == NULL || Y == NULL)
    {
        printf("Memory not allocated.");
        free(X);
        free(Y);
        return -1;
    }

    //double X[6] = {1.0, 2.0, 1.0, -3.0, 4.0, -1.0}; // 3x2
    X[0] = 1.0;
    X[1] = 2.0;
    X[2] = 1.0;
    X[3] = -3.0;
    X[4] = 4.0;
    X[5] = -1.0;

    Y[0] = 1.0;
    Y[1] = 2.0;
    Y[2] = -3.0;
    Y[3] = 4.0;

    distanceMatrix(X, Y, n, m, d);

    /* Deallocate used memory */
    free(X);
    free(Y);

    return 0;
}