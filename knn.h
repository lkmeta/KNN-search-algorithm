#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <cblas.h>
#include <time.h>
#include <limits.h>

#define SWAPD(x, y)      \
    {                    \
        double temp = x; \
        x = y;           \
        y = temp;        \
    } //swap two integers
#define SWAPI(x, y)   \
    {                 \
        int temp = x; \
        x = y;        \
        y = temp;     \
    } //swap two doubles

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

/**
 * Function printing matrix of doubles. Although it is stored as one dimensional array, actually
 * it is a 2-dimensional matrix stored in row major format.
 * We change line after we print each row of the matrix for clearer output.
 * Input:
 *      double* A: matrix to be printed
 *      int size: number of elements of the matrix
 *      int lineSize: number of elements in each row of the matrix
 * Output:
 *      None
**/
void printMatrix(double *A, int size, int lineSize)
{
    for (int i = 0; i < size; ++i)
    {
        if(i%lineSize==0 && i!=0){
            printf("\n");
        }
        printf("%lf ", A[i]);
    }
    printf("\n");
}

/**
 * Function creating random doubles within a specific range.
 * Input:
 *      double min: minimum value the random number can take
 *      double max: maximum value the random number can take
 * Output:
 *      the random number generated
**/
double randfrom(double min, double max) 
{
    double range = (max - min); 
    double div = RAND_MAX / range;
    return min + (rand() / div);
}

/**
 * Function creating a matrix of doubles with random values in [-100,100].
 * Input:
 *      double* A: the matrix we want to give the random values
 *      int size: number of elements of the matrix
 * Output:
 *      None
**/
void createRandomMatrix(double *A, int size)
{
    for (int i = 0; i < size; ++i)
    {
        A[i] = randfrom(-100.0, 100.0);
    }
}

/**
 * Function that for a given crpuse set X and query set Y calculates the Euclidean distance matrix D of
 * these points. The distance matrix D is calculated using the formula
 * D = sqrt(sum(X.^2,2) - 2 * X*Y.' + sum(Y.^2,2).')
 * X and Y lists are in treated in a row major format
 * In order to save some space we calculate the distance matrix row by row instead of all at once. This way
 * we don't have to compute sumX for all elements at once but one element each time.
 * Input:
 *      double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
 *      double* Y: mxd matrix containing the query data points (m data points with d coordinates each)
 *      int n: number of corpus points
 *      int m: number of query points
 *      int d: number of dimensions
 * Output:
 *      double* C: the Euclidean distance matrix (actually the D distance matrix)
 **/
double *distanceMatrix(double *X, double *Y, int n, int m, int d)
{   
    double dsquared;    //the square of the distance between two points

    double sumX;
    double *sumY = (double *)malloc(m * sizeof(double));

    //for the C=aA*B+b*C formula used in the cblas_dgemm function
    double alpha = 1.0;
    double beta = 0;

    //no need to create separate C and D matrixes, we calculate C = X*Y' first and then add to it sumX and sumY
    //in an appropriate way to get the D distance matrix (we save space this way)
    double *C = (double *)malloc(n * m * sizeof(double));

    if (C == NULL)
    {
        printf("Error in distanceMatrix: Couldn't allocate memory for C");
        exit(-1);
    }

    //Calculate sumY list
    for (int i = 0; i < m; i++)
    {
        sumY[i] = 0;
        for (int j = 0; j < d; j++)
        {
            sumY[i] += pow(Y[i * d + j], 2);
        }
    }

    //Calculate the corresponding row of the D matrix for each point in the corpus set X
    for (int i = 0; i < n; ++i)
    {

        //calculate sumX for the i-th corpus point each time
        sumX = 0;
        for (int k = 0; k < d; ++k)
        {
            sumX += pow(X[i * d + k], 2);
        }

        //calculate the i-th row of the matrix C=X*Y'
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, 1, m, d, alpha, &X[i * d], d, Y, d, beta, &C[i * m], m);

        //calculate the i-th row of the matrix D=(sumX)ee'−2C+ee'(sumY)'
        for (int j = 0; j < m; j++)
        {
            dsquared = sumX - 2 * C[i * m + j] + sumY[j];
            //due to rounding errors because of the way doubles are stored, distances slightly bigger or smaller
            //than 0 are rounded to 0
            if(dsquared<1e-5){
                dsquared =0;
            }
            C[i * m + j] = sqrt(dsquared);
        }
    }

    //deallocate used memory
    free(sumY);

    return C;
}

/** 
* Partition using Lomuto partition scheme
* We do the partition by applying comparisons to the elements of A but we move the elements of B in the same way
* so that they are rearranged in correspodence those of A.
* O(n) complexity
* Input:
*       double* A: first matrix in which we apply partition. Its elements are used in the comparisons with the pivot
*       int* B: second matrix in which we apply partition (rearranged the same way as A)
*       int left: left element of the array we are partitioning
*       int right: right element of the array we are partitioning
*       int pivotIndex: index of the chosen pivot element
* Output:
*       int pIndex: index of pivot element on the rearranged matrix
**/
int partition(double *A, int *B, int left, int right, int pivotIndex)
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
        //we allow for our doubles a tolerance error of 1e-6 (not every number can be stored exactly as a binary)
        if (A[i] - pivot < 1e-6)
        {
            SWAPD(A[i], A[pIndex]);
            SWAPI(B[i], B[pIndex]);
            pIndex++;
        }
    }

    // move pivot to its final place
    SWAPD(A[pIndex], A[right]);
    SWAPI(B[pIndex], B[right]);

    // return pIndex (index of pivot element)
    return pIndex;
}

/**
 * Function that implements the QuickSort algorithm. The algorithm is applied in the ndist array but the elements
 * of nidx are rearranged in an according way to match those of ndist.
 * Input:
 *      double* ndist: first matrix in which we apply quicksort.
 *      int* nidx: second matrix in which we apply quicksort (rearranged the same way as ndist) 
 *      int left: starting index, 
 *      int right: ending index
 * Output:
 *      None
**/
void quickSort(double* ndist, int* nidx, int left, int right) 
{ 
    if (left < right) 
    {
        // select a pivotIndex between left and right
        int pivotIndex = left + rand() % (right - left + 1);

        pivotIndex = partition(ndist, nidx, left, right, pivotIndex); 
  
        // Separately sort elements before 
        // partition and after partition 
        quickSort(ndist, nidx, left, pivotIndex - 1); 
        quickSort(ndist, nidx, pivotIndex + 1, right); 
    } 
}

/**
* Returns the k-th smallest element of list within left..right
* (i.e. left <= k <= right). The search space within the array is
* changing for each round - but the list is still the same size.
* Thus, k does not need to be updated with each round.
* O(logn) complexity
* Input:
*       double* A: first matrix in which we apply quickselect
*       int* B: A's corresponding matrix in which we apply quickselect
*       int left: left element of the array we are partitioning
*       int right: right element of the array we are partitioning
*       int k: the index of the element we are looking for (we are searching the k-th smallest element)
* Output:
*       double A[k]: the value of the k-th smallest element of matrix A
**/
void quickSelect(double *A, int *B, int left, int right, int k)
{
    // If the array contains only one element, return that element
    if (left == right)
        return;

    // select a pivotIndex between left and right
    int pivotIndex = left + rand() % (right - left + 1);

    pivotIndex = partition(A, B, left, right, pivotIndex);

    // The pivot is in its final sorted position
    if (k == pivotIndex)
        return;

    // if k is less than the pivot index
    else if (k < pivotIndex)
        return quickSelect(A, B, left, pivotIndex - 1, k);

    // if k is more than the pivot index
    else
        return quickSelect(A, B, pivotIndex + 1, right, k);
}

/**
 * Function that finds the knn of a specific point 'pointNum' of Y. Called if k>=n
 * Matrix ndist is arranged in such a way that each of the first k-1 elements of the matrix is smaller than the k-th
 * element. So, in order to get k nearest neighbors we just have to choose the first k elements of the particular column
 * in D. nidx contains the indexes of these elements in the original corpus set X and is arranged in such a way
 * that nidx[i] has the index of the point that its distance from the given query point is ndist[i]. This function
 * is called when the nearest neighbors we are looking for are more than the corpus set points. In this case,
 * we add all of the points in X and fill the remaining of the ndist with INFINITY and and nidx with INT_MAX
 * This way, if there is a new corpuse set X in the future, the infinity-valued elements will alway be bigger and
 * eventually removed when the new points come in (the distances from the query point will be finite).
 * O(n) complexity
 * Input:
 *      double* D: the all-to-all distance matrix
 *      int n: number of corpus set X points
 *      int m: number of query set Y points
 *      int k: number of nearest neighbors we are looking for
 *      int pointNum: index of the query point we are examining
 *      int* nidx: array containing the indexes of the points of set X (i-th element of nidx in correspondence
 *      to the i-th element of ndist) nearest to the query point
 *      double* ndist: array containing the distances of the k nearest points of X to the query point
 *      with index pointnum in Y
 * Output:
 *      None
**/
void addElems(double *D, int n, int m, int k, int pointNum, int *nidx, double *ndist)
{
    int elemCount = 0;
    for (int i = 0; i < n; ++i)
    {

        //ignore the element if it is same as the query point (distance between the two points is 0)
/*
        if (D[i * m + pointNum] < 0.001f)
        {
            if (i == n - 1)
                break;
            else
                continue;
        }
*/
        ndist[pointNum * k + elemCount] = D[i * m + pointNum];
        nidx[pointNum * k + elemCount] = i;
        elemCount++;
    }

    //in the remaining positions put these values so that they can be removed first when new corpus set points arrive
    for (int i = elemCount; i < k; ++i)
    {
        ndist[pointNum * k + i] = INFINITY;
        nidx[pointNum * k + i] = INT_MAX;
    }
}

/**
 * Function that finds the knn of a specific point 'pointNum' of Y. Called if k<n.
 * Matrix ndist is arranged in such a way that each of the first k-1 elements of the matrix is smaller than the k-th
 * element. So, in order to get k nearest neighbors we just have to choose the first k elements of the particular column
 * in D. nidx contains the indexes of these elements in the original corpus set X and is arranged in such a way
 * that nidx[i] has the index of the point that its distance from the given query point is ndist[i]. This function
 * is called when the nearest neighbors we are looking for are less than the corpus set points. In this case, we choose
 * the first k non-zero elements of A (and their corresponding indexes in X stored in the k first elements of B).
 * O(k) complexity
 * Input:
 *      double* A: matrix containing the distances of the corpus set points from the query point with index pointnum
 *      int* B: matrix containing the indexes of the points of X whose distances form the query point are in the
 *      corresponding position in A (i.e the index of the corpus set point, whose distance from the query point
 *      is A[i], is B[i])
 *      int n: number of corpus set X points
 *      int k: number of nearest neighbors we are looking for
 *      int pointNum: index of the query point we are examining
 *      int* nidx: array containing the indexes of the points of set X (i-th element of nidx in correspondence
 *      to the i-th element of ndist) nearest to the query point
 *      double* ndist: array containing the distances of the k nearest points of X to the query point
 *      with index pointNum in Y
 * Output:
 *      None
**/
void addElems2(double *A, int *B, int n, int k, int pointNum, double *ndist, int *nidx)
{
    int elemCount = 0;
    for (int i = 0; i <= k && elemCount < k; ++i)
    {

        //ignore the element if its value is zero (if the query point is compared to itself)
/*        if (A[i] < 0.001f)
        {
            if (i == n - 1)
                break;
            else
                continue;
        }
*/
        nidx[pointNum * k + elemCount] = B[i];
        ndist[pointNum * k + elemCount] = A[i];
        elemCount++;
    }
}

/**
* Function that finds the knn of a specific point 'pointNum' of Y.
* If n<k then just consider all elements of X as neighbors of the query point
* Else, create mPointDistanceMatrix which contains the distances of the query point's knn and mPointIndexMatrix which
* the indexes of those points in X and set ndist nidx equal to them respectively.
* O(n) complexity
* Input:
*       double* D, the all-to-all distance matrix
*       int pointNum: index of the query point we are examining
*       int n: number of corpus set X points
*       int m: number of query set Y points
*       int k: number of nearest neighbors we are looking for
*       int* nidx: array containing the indexes of the points of set X (i-th element of nidx in correspondence
*       to the i-th element of ndist) nearest to the query point
*       double* ndist: array containing the distances of the k nearest points of X to the query point
*       with index pointNum in Y
*       int flag: 0 if X and Y are different matrices, 0 if they are the same matrix
* Output:
*       None
**/
void kSelect(double *D, int pointNum, int n, int m, int k, int *nidx, double *ndist, int flag)
{

    if (k >= n)
    {
        quickSort(ndist, nidx, 0, n-1);
        addElems(D, n, m, k, pointNum, nidx, ndist);
    }
    else
    {
        //we create another subarray because if we apply quickselect in the original the elements change positions
        double *mPointDistanceMatrix = (double *)malloc(n * sizeof(double));
        int *mPointIndexMatrix = (int *)malloc(n * sizeof(int));

        if (mPointDistanceMatrix == NULL)
        {
            printf("Error in kSelect: Couldn't allocate memory for mPointDistanceMatrix");
            exit(-1);
        }

        if (mPointIndexMatrix == NULL)
        {
            printf("Error in kSelect: Couldn't allocate memory for mPointIndexMatrix");
            exit(-1);
        }

        //now mPointDistanceMatrix contains the distances the points of X have from point pointNum of Y
        //now mPointIndexMatrix contains the indexes of the points of X for a specific query point OF y
        for (int i = 0; i < n; ++i)
        {
            mPointDistanceMatrix[i] = D[i * m + pointNum];
            mPointIndexMatrix[i] = i;
        }

        //if X==Y then find k+1 smallest elements because we the same point cannot be considered as a neighbor of itself
        if (flag == 1)
        {
            quickSelect(mPointDistanceMatrix, mPointIndexMatrix, 0, n - 1, k);
        }
        //else find the k smallest elements
        else
        {
            quickSelect(mPointDistanceMatrix, mPointIndexMatrix, 0, n - 1, k - 1);
        }

        //sort the elements in nondecreasing order
        quickSort(mPointDistanceMatrix, mPointIndexMatrix,0,k-1);

        //create the final ndist and nidx for the query point with index pointNum
        addElems2(mPointDistanceMatrix, mPointIndexMatrix, n, k, pointNum, ndist, nidx);

        free(mPointDistanceMatrix);
        free(mPointIndexMatrix);
    }
}

/**
 * Function that implements the knn algorithm for a given query set Y and corpus set X.
 * At first we check if Y is the same as X, then we calculate the all-to-all distance matrix D
 * and then for each point in Y we find its k nearest neighbors
 * Input:
 *      double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
 *      double* Y: mxd matrix containing the query data points (m data points with d coordinates each)
 *      int n: number of corpus points
 *      int m: number of query points
 *      int d: number of dimensions
 *      int k: number of nearest neighbors we are looking for
 * Output:
 *      knnresult retVal: structure containing the info about the knn of each point of Y
**/
knnresult kNN(double *X, double *Y, int n, int m, int d, int k)
{

    int YisX;

    //check if X and Y are the same
    if (Y == X)
        YisX = 1;
    else
        YisX = 0;

    //calculate distance matrix
    double *D = distanceMatrix(X, Y, n, m, d);

    int *nidx = (int *)malloc(m * k * sizeof(int));

    if (nidx == NULL)
    {
        printf("Error in kNN: Couldn't allocate memory for nidx");
        exit(-1);
    }

    double *ndist = (double *)malloc(m * k * sizeof(double));

    if (ndist == NULL)
    {
        printf("Error in kNN: Couldn't allocate memory for ndist");
        exit(-1);
    }

    //for each point of the query set Y, find its knn
    for (int i = 0; i < m; ++i)
    {
        kSelect(D, i, n, m, k, nidx, ndist, YisX);
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

/**
 * Function implementing the insertion sort algorithm.
 * The algorithm is applied in the arr array but the elements of the indexes array
 * are rearranged in an according way to match those of arr.
 * Input:
 *      double* arr: first matrix in which we apply insertion sort.
 *      int* nidx: second matrix in which we apply insertion sort (rearranged the same way as arr) 
 *      int n: number of elements
 * Output:
 *      None
**/
void insertionSort(double *arr, int *indexes, int n)
{
    int i, j;
    double key;
    int keyIndex;

    for (i = 1; i < n; i++)
    {
        key = arr[i];
        keyIndex = indexes[i];
        j = i - 1;

        /* Move elements of arr[0..i-1], that are 
          greater than key, to one position ahead 
          of their current position */
        while (j >= 0 && arr[j] - key > -1e-6)
        {
            arr[j + 1] = arr[j];
            indexes[j + 1] = indexes[j];
            j--;
        }
        arr[j + 1] = key;
        indexes[j + 1] = keyIndex;
    }
}

/**
 * Function that finds the k smallest elements of the two lists. Combines the two lists into one and then applies
 * k-select to it. Then, we only keep the first k elements of this list which are (due to k-select)
 * the k smallest ones. The lists actually contain the distances of two different subsets of corpus points from
 * the query points. By doing this for every subset of corpus points, we eventually examine the whole corpus set
 * and get the final knn of each query point.
 * O(m) time complexity.
 * Input:
 *      knnresult old: struct containing the first list to be merged (updated list is stored to this one)
 *      knnresult new:struct containing the second list to be merged
 *      int m: number of query points
 *      int k: int k: number of nearest neighbors we are looking for
 *      int offset: offset of the corpus sub-set points in respect to the original corpus set
 * Output:
 *      None
**/
void mergeLists(knnresult old, knnresult new, int m, int k, int offset)
{

    //array containing the elements of both old.ndist and new.ndist arrays combined for a particular query point each time
    double *ndistComb = (double *)malloc(2 * k * sizeof(double));
    if (ndistComb == NULL)
    {
        printf("Error in mergeLists: Couldn't allocate memory for ndistComb");
        exit(-1);
    }

    //array containing the elements of both old.nidx and new.nidx arrays combined for a particular query point each time
    int *nidxComb = (int *)malloc(2 * k * sizeof(int));
    if (nidxComb == NULL)
    {
        printf("Error in mergeLists: Couldn't allocate memory for nidxComb");
        exit(-1);
    }

    //add the elemements to ndistComb and nidxComb
    for (int i = 0; i < m; ++i)
    {
        for (int j = 0; j < k; ++j)
        {
            ndistComb[j] = old.ndist[i * k + j];
            ndistComb[k + j] = new.ndist[i * k + j];
            nidxComb[j] = old.nidx[i * k + j];
            //we have to add offset to the indexes because they are in respect
            //to the corpus subset of the process an not the whole corpus set
            nidxComb[k + j] = new.nidx[i * k + j] + offset;
        }

        //find the k-th smallest element of ndistComb and rearrange ndistComb and nidxComb accordingly
        quickSelect(ndistComb, nidxComb, 0, 2 * k - 1, k-1);

        for (int j = 0; j < k; ++j)
        {
            old.ndist[i * k + j] = ndistComb[j];
            old.nidx[i * k + j] = nidxComb[j];
        }
    }

    //deallocate memory
    free(ndistComb);
    free(nidxComb);
}
