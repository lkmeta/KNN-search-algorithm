#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <cblas.h>
#include <omp.h>
#include <time.h>
#include <limits.h>

#define SWAPD(x,y) { double temp = x; x = y; y = temp; }
#define SWAPI(x,y) { int temp = x; x = y; y =temp; }


 // Definition of the kNN result struct
typedef struct knnresult
{
    int *nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double *ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;         //!< Number of query points                 [scalar]
    int k;         //!< Number of nearest neighbors            [scalar]
} knnresult;

// //! Compute k nearest neighbors of each point in X [n-by-d]
// /*!

//   \param  X      Corpus data points              [n-by-d]
//   \param  Y      Query data points               [m-by-d]
//   \param  n      Number of corpus points         [scalar]
//   \param  m      Number of query points          [scalar]
//   \param  d      Number of dimensions            [scalar]
//   \param  k      Number of neighbors             [scalar]

//function printing a matrix of doubles
void printMatrix(double* A, int size){
    for(int i=0;i<size;++i){
        printf("%lf ", A[i]);
    }
    printf("\n");
}

//function creating a matrix of doubles with random values
double createRandomMatrix(double* A, int size){
    for(int i=0;i<size;++i){
        A[i] = (rand() % 21);
    }
}

/**
 * Function that calculates an m×n Euclidean distance matrix D
 * D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).');
 * X and Y lists with CblasColMajor format
 * X[n*d] and Y[m*d]
 */
double* distanceMatrix(double *X, double *Y, int n, int m, int d)
{
    /* Allocate memory for sumX and sumY lists*/
    double *sumX = (double *)malloc(n * sizeof(double));
    double *sumY = (double *)malloc(m * sizeof(double));

    if(sumX==NULL){
        printf("Error in distanceMatrix: Couldn't allocate memory for sumX");
        exit(-1);
    }

    if(sumY==NULL){
        printf("Error in distanceMatrix: Couldn't allocate memory for sumX");
        exit(-1);
    }

    /* Calculate sumX list*/
    //mporei parallhla
    for (int i = 0; i < n; i++)
    {
        sumX[i]=0;
        for (int j = 0; j < d; j++)
        {
            sumX[i] += pow(X[i*d+j], 2);
        }
    }

    /* Calculate sumY list*/
    //mporei parallhla
    for (int i = 0; i < m; i++)
    {
        sumY[i]=0;
        for (int j = 0; j < d; j++)
        {
            sumY[i] += pow(Y[i*d+j], 2);
        }
    }

//    printf("\nSumX = ");
//    printMatrix(sumX, n);

//    printf("\nSumY = ");
//    printMatrix(sumY, m);
    
    /* Allocate memory for C and D lists*/
    double alpha = 1.0;
    double beta = 0;
    double *C = (double *)malloc(n * m * sizeof(double));

    if(C==NULL){
        printf("Error in distanceMatrix: Couldn't allocate memory for C");
        exit(-1);
    }

    /* Calculate list C = X*Y */
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, n, m, d, alpha, X, d, Y, d, beta, C, m);

//    printf("\nC = ");
//    printMatrix(C, n*m);

    /* Calculate list D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).'); */
    #pragma omp parallel for schedule(dynamic)
    for(int i=0; i<n; i++)
    {
        for (int j = 0; j < m; j++)
        {
            C[i*m+j] = sqrt(sumX[i] - 2 * C[i*m+j] + sumY[j]);
        }
    }   

    printf("\n\nD = ");
    printMatrix(C, n*m);

    /* Deallocate used memory */
    free(sumX);
    free(sumY);

    return C;
}

// Partition using Lomuto partition scheme
int partition(double* A, int* B, int left, int right, int pivotIndex)
{
    // Pick pivotIndex as pivot from the array
    double pivot = A[pivotIndex];
 
    // Move pivot to end
    SWAPD(A[pivotIndex], A[right]);
    SWAPI(B[pivotIndex], B[right]);
 
    // elements less than pivot will be pushed to the left of pIndex
    // elements more than pivot will be pushed to the right of pIndex
    // equal elements can go either way
    int pIndex = left;
 
    // each time we finds an element less than or equal to pivot, pIndex
    // is incremented and that element would be placed before the pivot.
    for (int i = left; i < right; i++)
    {
        if (A[i]-pivot<0.001f)
        {
            SWAPD(A[i], A[pIndex]);
            SWAPI(B[i], B[pIndex]);
            pIndex++;
        }
    }
 
    // Move pivot to its final place
    SWAPD(A[pIndex], A[right]);
    SWAPI(B[pIndex], B[right]);
 
    // return pIndex (index of pivot element)
    return pIndex;
}

// Returns the k-th smallest element of list within left..right
// (i.e. left <= k <= right). The search space within the array is
// changing for each round - but the list is still the same size.
// Thus, k does not need to be updated with each round.
double quickSelect(double* A, int* B, int left, int right, int k) //change to void later, keep it now for testing
{
    // If the array contains only one element, return that element
    if (left == right)
        return A[left];
 
    // select a pivotIndex between left and right
    int pivotIndex = left + rand() % (right - left + 1);
 
    pivotIndex = partition(A, B, left, right, pivotIndex);
 
    // The pivot is in its final sorted position
    if (k == pivotIndex)
        return A[k];
 
    // if k is less than the pivot index
    else if (k < pivotIndex)
        return quickSelect(A, B, left, pivotIndex - 1, k);
 
    // if k is more than the pivot index
    else
        return quickSelect(A, B, pivotIndex + 1, right, k);
}

/*
//for a specific column of D find the k smallest elements and sore their values in ndist and their indices in nidx;
void checkElems(double* D, int n, int m, int k, double kElem, int pointNum, int* nidx, double* ndist){

    int elemCount=0;

    int left = pointNum*n;
    int right = (pointNum+1)*n - 1;

    for(int i=0;i<n && elemCount<k;++i){ //no need to keep checking the elements if we have already found the k smallest ones
        if(D[i*m+pointNum]-kElem<0.001f){
            if(D[i*m+pointNum]<0.001f){
                if(i==n-1)break;
                else continue;
            }
            nidx[pointNum*k+elemCount] = i; 
            ndist[pointNum*k+elemCount] = D[i*m+pointNum];
            elemCount++;
        }
    }
}
*/

//Function that finds the knn of a specific point 'pointNum' of Y. Called if k>=n
void addElems(double* D, int n, int m, int k, int pointNum, int* nidx, double* ndist){
    int elemCount=0;
    for(int i=0;i<n;++i){
        if(D[i*m+pointNum]<0.001f){
            if(i==n-1)break;
            else continue;
        }
        ndist[pointNum*k+elemCount] = D[i*m+pointNum];
        nidx[pointNum*k+elemCount] = i;
        elemCount++;
    }
    for(int i=elemCount;i<k;++i){
        ndist[pointNum*k+i] = INFINITY;
        nidx[pointNum*k+i] = INT_MAX;
    }
}

void addElems2(double* A, int* B,int n, int k, int pointNum, double* ndist, int* nidx){
    int elemCount=0;
    for(int i=0;i<=k && elemCount<k;++i){
        if(A[i]<0.001f){
            if(i==n-1) break;
            else continue;
        }
        nidx[pointNum*k+elemCount] = B[i];
        ndist[pointNum*k+elemCount] = A[i];
        elemCount++;
    }
}

/*
Function that finds the knn of a specific point 'pointNum' of Y. Called if k<n
*/
void kSelect(double* D,int pointNum, int n, int m, int k, int* nidx, double* ndist, int flag){

    if(k>=n){
        //O(n)
        addElems(D,n,m,k,pointNum,nidx,ndist);
    }
    else{
        //we create another subarray because if we apply quickselect in the original the elements might change positions
        double* mPointDistanceMatrix = (double*)malloc(n*sizeof(double));
        int* mPointIndexMatrix = (int*)malloc(n*sizeof(int));

        if(mPointDistanceMatrix==NULL){
            printf("Error in kSelect: Couldn't allocate memory for mPointDistanceMatrix");
            exit(-1);
        }

        if(mPointIndexMatrix==NULL){
            printf("Error in kSelect: Couldn't allocate memory for mPointIndexMatrix");
            exit(-1);
        }

        //now mPointDistanceMatrix contains the distances the points of X have from point pointNum of Y
        //now mPointIndexMatrix contains the indexes of the points of X for a specific query point OF y
        //O(n)
        for(int i=0;i<n;++i){
            mPointDistanceMatrix[i] = D[i*m+pointNum];
            mPointIndexMatrix[i] = i;
        }

//        printf("mPointDistanceMatrix=\n");
//        printMatrix(mPointDistanceMatrix,n);

        double kElem;

        //O(n)
        if(flag==1){
            kElem = quickSelect(mPointDistanceMatrix,mPointIndexMatrix,0,n-1,k);
        }
        else{
            kElem = quickSelect(mPointDistanceMatrix,mPointIndexMatrix,0,n-1,k-1);
        }
//        printf("kElem=%lf\n", kElem);

        //O(n)
        //checkElems(D,n,m,k,kElem,pointNum,nidx,ndist);

        //O(k)
        addElems2(mPointDistanceMatrix,mPointIndexMatrix,n,k,pointNum,ndist,nidx);

        free(mPointDistanceMatrix);
        free(mPointIndexMatrix);
    }
}

//function that implements the knn algorithm for a given query set Y and corpus set X
knnresult kNN(double* X, double* Y, int n, int m, int d, int k){
    
    int YisX;

    if(Y==X) YisX = 1;
    else YisX=0;

    double *D = distanceMatrix(X,Y,n,m,d);

    int* nidx = (int *)malloc(m * k * sizeof(int));

    if(nidx==NULL){
        printf("Error in kNN: Couldn't allocate memory for nidx");
        exit(-1);
    }

    double* ndist = (double *)malloc(m * k * sizeof(double));

    if(ndist==NULL){
        printf("Error in kNN: Couldn't allocate memory for ndist");
        exit(-1);
    }

    //for each point of the query set Y, find its knn
    for(int i=0;i<m;++i){
//        printf("\nm=%d\n", i);
        kSelect(D,i,n,m,k,nidx,ndist,YisX);
    }

    free(D);

    //struct to be returned
    knnresult retVal;

    retVal.k = k;
    retVal.m = m;
    retVal.ndist = ndist;
    retVal.nidx = nidx;

    return retVal;
}


/*
int main(int argc, char* argv[])
{   
    //currently change them by hand, will fix later
    int n = 4, m = 3, d = 2, k = 3;

    srand(time(NULL));

    int threadNum=atoi(argv[1]);    //number of threads
    printf("\nYou have chosen %d threads \n",threadNum);

    // Allocate memory for X and Y lists
    double *X = (double *)malloc(n * d * sizeof(double));
    double *Y = (double *)malloc(m * d * sizeof(double));

    if(X==NULL){
        printf("Error in main: Couldn't allocate memory for X");
        exit(-1);
    }

    if(Y==NULL){
        printf("Error in main: Couldn't allocate memory for Y");
        exit(-1);
    }

    //create set X and y
    createRandomMatrix(X,n*d);
    createRandomMatrix(Y,m*d);

    printf("X= \n");
    printMatrix(X, n*d);

    printf("Y= \n");
    printMatrix(Y, m*d);

    knnresult result = kNN(X,X,n,n,d,k);

    printf("\nndist = \n");
    printMatrix(result.ndist, m*k);

    printf("nidx = \n");
    for(int i=0;i<m*k;++i){
        printf("%d ", result.nidx[i]);
    }
    printf("\n");

    // Deallocate used memory 
    free(X);
    free(Y);

    free(result.nidx);
    free(result.ndist);

    return 0;
}
*/