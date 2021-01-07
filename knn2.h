#ifndef KNN2_H
#define KNN2_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define SWAPD(x,y) { double temp = x; x = y; y = temp; }    //swap two doubles
#define SWAPI(x,y) { int temp = x; x = y; y =temp; }        //swap two integers

typedef struct queryPoint queryPoint;
typedef struct Node Node;
typedef struct knnresult knnresult;


// Definition of the query point struct
struct queryPoint
{
    double *coord;    //!< d coords for query point               [1-by-d]
    int d;            //!< Number of coords
    int *nidx;        //!< Indices (0-based) of nearest neighbors [m-by-k]
    double *ndist;    //!< Distance of nearest neighbors          [m-by-k]
    int k;            //!< Number of nearest neighbors            [scalar]
    int numOfIndexes; //!< Counter for Indexes in a queryPoint
    int flag;         //!< Counter for the visited nodes from the VPT
    double tau;       //!< Radius for the searchVPT process
};


// Definition of the node struct
struct Node
{
    Node *left;         //left child of the node
    Node *right;        //right child of the node
    int p;              //index of vantage point
    double mu;          //median distance
    double *dists;      //matrix containing the distances of the points of the node from its vantage point
    int *indx;          //matrix containing the indexes of the points of the node
    int numOfIndexes;   //number of points contained in the node
};


// Definition of the kNN result struct
struct knnresult
{
    int *nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double *ndist; //!< Distance of nearest neighbors          [m-by-k]
    int m;         //!< Number of query points                 [scalar]
    int k;         //!< Number of nearest neighbors            [scalar]
};


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
         if (A[i]-pivot<1e-6)
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
*       None
**/
void quickSelect(double* A, int* B, int left, int right, int k)
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
 * Function implementing the insertion sort algorithm.
 * The algorithm is applied in the arr array but the elements of the indexes array
 * are rearranged in an according way to match those of arr.
 * O(n^2) time complexity but it is more efficient than other sorting algorithms if we deal with a small
 * number of elements or the elements are almost sorted already
 * Input:
 *      double* arr: first matrix in which we apply insertion sort.
 *      int* indexes: second matrix in which we apply insertion sort (rearranged the same way as arr) 
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
 * Function checking whether a specific point has been sampled already by checking whether the index
 * of this point belongs already in the array containing the indexes of the sampled points
 * Input:
 *      int* sampleIndex: array containing the indexes of the sampled points (-1 means that this position has not
 *      been filled with a sampled point yet)
 *      int sampleSize: total number of points to be sampled
 *      int index: index of the point that we want to check if it has been sampled
 * Output:
 *      int result: 0 if not sampled already, 1 otherwise
**/
int sampledAlready(int *sampleIndex, int sampleSize, int index)
{
    int result = 0;
    
    for (int i = 0; i < sampleSize; ++i)
    {
        //no need to check after that because these positions have not been filled yet
        if (sampleIndex[i] == -1)
        {
            break;
        }

        //if already sampled, stop searching the rest of the sampled points
        if (index == sampleIndex[i])
        {
            result = 1;
            break;
        }
    }

    return result;
}


/**
 * Function that creates a sample of points by choosing random points from the whole set of points.
 * Structured in such a way that all of the randomly selected points are unique.
 * Input:
 *      double* X: array containing the coordinates of all the points we can choose from
 *      double* sample: array containing the coordinates of the sampled points
 *      int* indexes: array containing the indexes of the points contained in X (because X itself might be a subset
 *      of the set containing the whole points)
 *      int sampleSize: number of sample points
 *      int n: number of points in X
 *      ind d: number of dimensions
 * Output:
 *      int* sampleIndex: indexes of the sampled points in the matrix containing all the points
**/
int *sampleSet(double *X, double *sample, int *indexes, int sampleSize, int n, int d)
{
    int count = 0;  //counts how mny elements have been sampled already
    int index;      //randomly chosen index 

    int *sampleIndex = (int *)malloc(sampleSize * sizeof(int));
    if (sampleIndex == NULL)
    {
        printf("Error in sampleSet: Couldn't allocate memory for sampleIndex");
        exit(-1);
    }

    //initialize the index array with -1 to declare that these positions are to be filled
    for (int i = 0; i < sampleSize; ++i)
    {
        sampleIndex[i] = -1;
    }

    while (count < sampleSize)
    {
        index = rand() % n;

        //if it has not been sampled already, add it to the sampled points
        if (sampledAlready(sampleIndex, sampleSize, indexes[index]) == 0)
        {
            sampleIndex[count] = indexes[index];
            for (int j = 0; j < d; ++j)
            {
                sample[count * d + j] = X[index * d + j];
            }
            count++;
        }
    }

    return sampleIndex;
}


/**
 * Function that calculates median of a list.
 * O(nlogn) time complexity.
 * Input:
 *      double *sampleDistances: array containing the values out of which we will choose the median distance
 *      int *indexes: array containing the indexes of the points who produced the distances stored in the
 *      sampleDistances array
 *      int sampleSize: number of sampled points
 * Output:
 *      double mu: median value
**/
double findMedian(double *sampleDistances, int *indexes, int sampleSize)
{
    double mu;

    insertionSort(sampleDistances, indexes, sampleSize);

    int middle = (sampleSize + 1) / 2 - 1;

    //if number of elements is not divisible by two choose the middle point
    if (sampleSize % 2 == 1)
    {
        mu = sampleDistances[middle];
    }

    //otherwise calculate the mean value of the two middle points
    else
    {
        mu = (sampleDistances[middle] + sampleDistances[middle + 1]) / 2;
    }

    return mu;
}


/**
 * Function that selects the vantage point. First we sample a number of points. Then for each one of these
 * sampled points, we create a second sample set and calculate the 2nd moment of this point based on the sample
 * points of the second set. Finally we keep the point of the first sample set that has the biggest 2nd moment.
 * This point is the most appropriate to become a vantage point since bigger spread (2nd moment) means that it is
 * more likely of it to be to the edge of the region defined by the points. This method does not ensure that
 * we always choose the best vantage point but produces better results than picking a random point as the vantage
 * point without increasing the time complexity too much.
 * Input:
 *      double* X: array containing the coordinates of all the points we can choose from
 *      int* indexes: array containing the indexes of the points contained in X (because X itself might be a subset
 *      of the set containing the whole points)
 *      int n: number of points in X
 *      int d: number of dimensions
 * Output:
 *      int bestP: the index of the chosen vantage point in the array containing the whole points
**/
int selectVP(double *X, int *indexes, int n, int d)
{

    //number of sample points is set to 10
    int sampleSize = 10;

    //unless the total number of points is less
    while (sampleSize >= n)
    {
        sampleSize /= 2;
    }

    //first set of sampled points
    double *P = (double *)malloc(sampleSize * d * sizeof(double));
    if (P == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for P");
        exit(-1);
    }

    //second set of sampled points
    double *D = (double *)malloc(sampleSize * d * sizeof(double));
    if (D == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for D");
        exit(-1);
    }

    //distances between a specific point of the first sample set and the points of the second sample set
    double *sampleDistances = (double *)malloc(sampleSize * sizeof(double));
    if (sampleDistances == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for sampleDistances");
        exit(-1);
    }

    //indexes of the of the elements sampled on the first sample set
    int *sampleIndex = sampleSet(X, P, indexes, sampleSize, n, d);

    int middle;
    double mu;
    double bestSpread = 0.0;
    double spread;
    int bestP;      //index of the point which is more appropriate to be a vantage point

    //for each point of the first sample set
    for (int i = 0; i < sampleSize; ++i)
    {
        //indexes of the of the elements sampled on the second sample set
        int *sampleIndex2 = sampleSet(X, D, indexes, sampleSize, n, d);

        //for each point of the second sample set calculate its distance from the specific point of the first
        //sample set
        for (int j = 0; j < sampleSize; ++j)
        {
            sampleDistances[j] = 0;
            for (int k = 0; k < d; ++k)
            {
                sampleDistances[j] += pow(P[i * d + k] - D[j * d + k], 2);
            }
            if(sampleDistances[j]<1e-5){
                sampleDistances[j] =0;
            }
            sampleDistances[j] = sqrt(sampleDistances[j]);
        }
        
        mu = findMedian(sampleDistances, sampleIndex2, sampleSize);
        free(sampleIndex2);

        spread = 0;

        //calculate the 2nd moment for this point of the first sample set
        for (int j = 0; j < sampleSize; ++j)
        {
            spread += pow(sampleDistances[j] - mu, 2);
        }
        spread /= sampleSize;
        
        //keep the point with the biggest spread value so far
        if (spread > bestSpread)
        {
            bestSpread = spread;
            bestP = sampleIndex[i];
        }
    }

    free(P);
    free(D);
    free(sampleDistances);
    free(sampleIndex);

    return bestP;
}


/**
 * Function that finds the index of the vantage points on the subset of points used to create a specific subtree
 * of the whole VPT. The set S of points used as an argument in the function making the VPT is not always the same
 * as the set X that contains all the points given. That is why it is important to distinguish the difference
 * between the index of a point in S and the index of a point in X. The value returned from the function that 
 * chooses the vantage point is the index of this point in X. This function finds the index of this vantage point
 * in S. This function is useful when we need to calculate distances from the vantage point, where we only have
 * the S matrix and not the whole x matrix.
 * Input:
 *      int vp: the index of the vantage point in the matrix containing all the points of the problem
 *      int* indexes: the indexes of the points contained in the given sub set of points to create the subtree
 *      int n: number of points in the subset
 * Output:
 *      int VPindex: index of the vantage point in the subset
**/
int findVPIndx(int vp, int *indexes, int n)
{
    int VPindex;
    for (int i = 0; i < n; ++i)
    {
        if (indexes[i] == vp)
        {
            VPindex = i;
            break;
        }
    }
    return VPindex;
}


/**
 * Function that creates the vantage point tree. First, the function checks whether it has to create a leaf or a
 * non-leaf node. If the node created is a leaf, we simply add all the given points to the leaf and set the index
 * of the vantage point to be -1 (works as a flag to separate the leafs from non-leaf nodes). If the node created
 * is not a leaf then we choose an appropriate vantage point, find the median of the points and separate the points
 * into two categories: those closer to the vantage point than the median that go to the left child and those farther
 * from the vantage point than the median that go to the right child. Finally, we call makeVPT to create the left
 * and then the right subtree.
 * Input:
 *      double* S: matrix containing the coordinates of the points out of which we are going to create the particular
 *      node and all of its sub-nodes (if there are any)
 *      int n: number of corpus set X points
 *      int d: number of dimensions
 *      int* indexes: the indexes of the given points on the initial point set containing all the points
 *      int B: the least number of points each leaf node may contain
 * Output:
 *      Node* nd: root of the tree created
**/
Node *makeVPT(double *S, int n, int d, int *indexes, int B)
{
    Node *nd;

    if (n==0)
    {
        nd = NULL;
        return nd;
    }

    nd = (Node *)malloc(sizeof(Node));
    if (nd == NULL)
    {
        printf("Error in makeVPT: Couldn't allocate memory for nd");
        exit(-1);
    }

    //if the points are not enough to be separated into a parent node and two children add them all in one node
    if (n < 2 * B + 1)
    {
        nd->p = -1;
        nd->mu = -1.0;
        nd->left = NULL;
        nd->right = NULL;
        nd->dists = (double *)malloc(n * sizeof(double));
        if (nd->dists == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->dists");
            exit(-1);
        }
        
        nd->indx = (int *)malloc(n * sizeof(int));
        if (nd->indx == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->indx");
            exit(-1);
        }

        memcpy(nd->indx, indexes, n * sizeof(int));

        // add length of indx array
        nd->numOfIndexes = n;
    }

    //otherwise create a arent node and its two children
    else
    {
        //find an appropriate vantage point
        nd->p = selectVP(S, indexes, n, d);

        //the vp's index on the S matrix containing the points to create the tree
        int pInd = findVPIndx(nd->p, indexes, n);

        nd->dists = (double *)malloc(n * sizeof(double));
        if (nd->dists == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->dists");
            exit(-1);
        }

        //calculate the distances of the rest of the points from the vantage point
        for (int i = 0; i < n; ++i)
        {
            nd->dists[i] = 0;
            for (int j = 0; j < d; ++j)
            {
                nd->dists[i] += pow(S[i * d + j] - S[pInd * d + j], 2);
            }
            if(nd->dists[i]<1e-5){
                nd->dists[i] = 0;
            }
            nd->dists[i] = sqrt(nd->dists[i]);
        }

        nd->indx = (int *)malloc(n * sizeof(int));
        if (nd->indx == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->indx");
            exit(-1);
        }

        memcpy(nd->indx, indexes, n * sizeof(int));

        //since it is not a leaf, the node will contain only one point, the vantage point
        nd->numOfIndexes = 1;

        //we use a temporary array for the distances because calculating the median rearranges the elements in
        //an unwanted way
        double *distsTemp = (double *)malloc(n * sizeof(double));
        if (distsTemp == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for distsTemp");
            exit(-1);
        }

        memcpy(distsTemp, nd->dists, n * sizeof(double));

        nd->mu = findMedian(distsTemp, indexes, n);

        free(distsTemp);

        double *right;      //coordinates of the points that will go to the right child
        int *rightIndexes;  //indexes of the points that will go to the right child
        
        double *left;       //coordinates of the points that will go to the left child
        int *leftIndexes;   //indexes of the points that will go to the left child

        int lSize = 0;      //number of points going to the left child
        int rSize = 0;      //number of points going to the right child

        //check how many points are closer to the vp than the median or farther from it
        for (int i = 0; i < n; ++i)
        {
            //exclude comparison of vantage point with itself
            if (nd->dists[i] < 1e-7)
            {
                if (i == n - 1)
                    break;
                else
                    continue;
            }

            if (nd->dists[i] - nd->mu < 1e-7)
            {
                lSize++;
            }
            else
            {
                rSize++;
            }
        }

        right = (double *)malloc(rSize * d * sizeof(double));
        if (right == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for right");
            exit(-1);
        }

        rightIndexes = (int *)malloc(rSize * sizeof(int));
        if (rightIndexes == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for rightIndexes");
            exit(-1);
        }

        left = (double *)malloc(lSize * d * sizeof(double));
        if (left == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for left");
            exit(-1);
        }

        leftIndexes = (int *)malloc(lSize * sizeof(int));
        if (leftIndexes == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for leftIndexes");
            exit(-1);
        }

        int lCounter = 0;   //counts how many elements have been added to the left child
        int rCounter = 0;   //counts how many elements have been added to the right child

        //separate elements going to the left or the right child
        for (int i = 0; i < n; ++i)
        {
            //exclude vantage point from being added to one of the children
            if (nd->dists[i] < 1e-7)
            {
                if (i == n - 1)
                    break;
                else
                    continue;
            }

            if (nd->dists[i] - nd->mu < 1e-7)
            {
                if (lCounter >= lSize)
                    break;
                else
                {
                    leftIndexes[lCounter] = nd->indx[i];
                    for (int j = 0; j < d; ++j)
                    {
                        left[lCounter * d + j] = S[i * d + j];
                    }
                    lCounter++;
                }
            }
            else
            {
                if (rCounter >= rSize)
                    break;
                else
                {
                    rightIndexes[rCounter] = nd->indx[i];
                    for (int j = 0; j < d; ++j)
                    {
                        right[rCounter * d + j] = S[i * d + j];
                    }
                    rCounter++;
                }
            }
        }

        //create the left subtree
        nd->left = makeVPT(left, lSize, d, leftIndexes, B);

        //create the right subtree
        nd->right = makeVPT(right, rSize, d, rightIndexes, B);

        free(right);
        free(rightIndexes);
        free(left);
        free(leftIndexes);
    }

    return nd;
}


/**
 * Function that adds the neighbors from the node into the query point.
 * First, we add all the neighbors from current node into the query point.
 * Then we sort the list with the distances. Finally, we reallocate the list 
 * and we keep the nearest neighbors.
 * 
 * Input:
 *      Node *currentNode    : current node in the tree
 *      queryPoint *currentP : current query point we are looking for its k nearest neighbors 
 *      double *S            : nxd matrix containing the corpus data points (n points with d coordinates each)
 * Output:
 *      None
**/
void addElements(Node *currentNode, queryPoint *currentP, double *S, int offset)
{

    //check if Node is empty
    if (currentNode->numOfIndexes == 0)
    {
        printf("Error in addElements: No elements.");
        return;
    }

    int length = currentNode->numOfIndexes;          // current length
    int newLength = length + currentP->numOfIndexes; // possible next length

    //reallocate memory because of new possible length
    currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * newLength);
    currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * newLength);

    //add new indexes and distances from curent node into the query point
    for (int i = 0; i < length; i++)
    {
        //if the node is not a leaf add the vantage point
        if (currentNode->p != -1)
        {
            //temporary idx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->p;

            //calculate temporary dist to add
            currentP->ndist[i + currentP->numOfIndexes] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                currentP->ndist[i + currentP->numOfIndexes] += pow(S[currentP->d * currentNode->p -offset + j] - currentP->coord[j], 2);
            }
            if(currentP->ndist[i + currentP->numOfIndexes]<1e-5){
                currentP->ndist[i + currentP->numOfIndexes] = 0;
            }
            currentP->ndist[i + currentP->numOfIndexes] = sqrt(currentP->ndist[i + currentP->numOfIndexes]);
        }

        //if the node is a leaf, iclude all the corpus set points of the leaf a neighbors
        else
        {
            //temporary idx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->indx[i];

            //calculate temporary dist to add
            currentP->ndist[i + currentP->numOfIndexes] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                currentP->ndist[i + currentP->numOfIndexes] += pow(S[currentP->d * currentNode->indx[i] - offset + j] - currentP->coord[j], 2);
            }
            if(currentP->ndist[i + currentP->numOfIndexes]<1e-5){
                currentP->ndist[i + currentP->numOfIndexes] = 0;
            }
            currentP->ndist[i + currentP->numOfIndexes] = sqrt(currentP->ndist[i + currentP->numOfIndexes]);
        }
    }

    //sort the new ndist and nidx lists created
    if (newLength != 1)
    {
        insertionSort(currentP->ndist, currentP->nidx, newLength);
    }

    //find the real number of indexes that we will keep
    for (int i = 0; i < length; i++)
    {
        if (currentP->k - currentP->numOfIndexes > 0)
        {
            currentP->numOfIndexes++;
        }
    }

    //reallocate the length of distances and indexes lists according to the number of neighbors we have found
    if (currentP->numOfIndexes != newLength)
    {
        currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * (currentP->numOfIndexes));
        currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * (currentP->numOfIndexes));
    }

    //update tau to be the distance from the farthest neighbor so far
    if (currentP->k == currentP->numOfIndexes)
    {
        currentP->tau = currentP->ndist[currentP->numOfIndexes - 1];
    }

}

/**
* Function that searches for the k nearest distances from the query point in the VPT.
* First, we check if we are searching in a leaf node or not. 
* When we are in a leaf node we push the nearest neighbors in query point and we return it.
* When we are not in a leaf node, we check for the following possible cases:
*   1. if tempdist is lower than the current radius we use addElements
*   2. if tempdist is lower than the current median we search in left child
*   3. if tempdist is higher than the current median we search in right child
* Finally, we check for intersection case of the radius tau and the radius defined by the median
* of a node and then we return the currentP with all the k nearest neighbors. 
* Input:
*       Node *root         : root node of the VPT
*       queryPoint *queryP : query point for which we are looking the k nearest neighbors 
*       double *S          : nxd matrix containing the corpus data points (n points with d coordinates each)
* Output:
*       None 
**/
void searchVPT(Node *root, queryPoint *queryP, double *S, int offset)
{

    if (root == NULL)
    {
        printf("Error in searchVPT: Root is NULL.");
        exit(-1);
    }

    if (queryP == NULL)
    {
        printf("Error in searchVPT: QueryPoint is NULL.");
        exit(-1);
    }

    //increase the counter for the visited nodes
    queryP->flag++;

    //if the node is a leaf
    if (root->right == NULL && root->left == NULL)
    {
        //push nearest neighbors in query point
        addElements(root, queryP, S, offset);
        return;
    }
    //if the node is not a leaf
    else
    {

        double tempDist = 0;

        //calculating distance between queryPoint and vp
        for (int i = 0; i < queryP->d; ++i)
        {
            tempDist += pow(S[queryP->d * root->p - offset + i] - queryP->coord[i], 2);
        }
        if(tempDist<1e-5){
                tempDist = 0;
            }
        tempDist = sqrt(tempDist);

        //if the vantage point is inside the tau radius consider it as a neighbor
        if (tempDist < queryP->tau)
        {
            addElements(root, queryP, S, offset);
        }

        //if the query point is closer to the points of the left child
        if (tempDist < root->mu)
        {
            searchVPT(root->left, queryP, S, offset);

            //if we haven't found all neighbors or there is an intersection between the median radius and tau radius
            //check the other child too
            if (queryP->k > queryP->numOfIndexes || tempDist > root->mu - queryP->tau)
            {
                searchVPT(root->right, queryP, S, offset);
            }
        }

        //if the query point is closer to the points of the right child
        if (tempDist >= root->mu)
        {
            searchVPT(root->right, queryP, S, offset);

            //if we haven't found all neighbors or there is an intersection between the median radius and tau radius
            //check the other child too
            if (queryP->k > queryP->numOfIndexes || tempDist < root->mu + queryP->tau)
            {
                searchVPT(root->left, queryP, S, offset);
            }
        }
    }

    return;
}


/**
 * Function that implements the knn algorithm for a given query set Y and corpus set X.
 * For each point in Y we create its queryPoint struct and then we search in the VPT for its k nearest neighbors
 * Input:
 *      double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
 *      double* Y: mxd matrix containing the query data points (m data points with d coordinates each)
 *      Node* root: pointer of type Node to the root of the VPTree
 *      int n: number of corpus points
 *      int m: number of query points
 *      int d: number of dimensions
 *      int k: number of nearest neighbors we are looking for
 *      int offset: offset of the beginning of the local corpus set X in respect to the total corpus set X
 *      containing all the corpus points
 * Output:
 *      knnresult retVal: structure containing the info about the knn of each point of Y
**/
knnresult kNN(double *X, double *Y, Node* root, int n, int m, int d, int k, int offset)
{
    
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

    queryPoint *p = (queryPoint*)malloc(sizeof(queryPoint));
    if (p == NULL)
    {
        printf("failed to allocate memory for query point\n");
        exit(-1);
    }
    
    double *tempCoord = (double *)malloc(d * sizeof(double));
    if (tempCoord == NULL)
    {
        printf("failed to allocate memory for tempCoord\n");
        exit(-1);
    }

    //for each point of the query set Y, find its knn using searchVPT
    for (int j = 0; j < m; j++)
    {
        //initialize the parameters of the queryPoint p
        p->d = d;
        p->k = k;
        p->flag = 0;
        p->tau = INFINITY;

        for (int i = 0; i < d; i++)
        {
            tempCoord[i] = Y[j*d + i];
        }

        p->coord = tempCoord;
        p->nidx = NULL;
        p->ndist = NULL;
        p->numOfIndexes = 0;

        //search for its kNN in the local tree
        searchVPT(root, p, X, offset);

        for (int i = 0; i < p->numOfIndexes; i++)
        {
            nidx[j*k+i] = p->nidx[i];
        }

        for (int i = 0; i < p->numOfIndexes; i++)
        {
            ndist[j*k+i] = p->ndist[i];
        }

    }

    free(p->nidx);
    free(p->ndist);
    free(tempCoord);
    free(p);

    //struct to be returned
    knnresult retVal;

    retVal.k = k;
    retVal.m = m;
    retVal.ndist = ndist;
    retVal.nidx = nidx;

    return retVal;
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
void mergeLists(double* ndistOld,int* nidxOld, double* ndistNew, int* nidxNew, int m, int k, int offset)
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
            ndistComb[j] = ndistOld[i * k + j];
            ndistComb[k + j] = ndistNew[i * k + j];
            nidxComb[j] = nidxOld[i * k + j];
            //we have to add offset to the indexes because they are in respect
            //to the corpus subset of the process an not the whole corpus set
            nidxComb[k + j] = nidxNew[i * k + j] + offset;
        }

        //find the k-th smallest element of ndistComb and rearrange ndistComb and nidxComb accordingly
        quickSelect(ndistComb, nidxComb, 0, 2 * k - 1, k - 1);

        for (int j = 0; j < k; ++j)
        {
            ndistOld[i * k + j] = ndistComb[j];
            nidxOld[i * k + j] = nidxComb[j];
        }
    }

    //deallocate memory
    free(ndistComb);
    free(nidxComb);
}

/**
 * Function freeing all the memory allocated to create the VPTree. Works in a recursive function, freeing
 * the allocated memory first from the leafs and then going back up to the root.
 * Input:
 *      Node* node: the node which we want to deallocate its memory
 * Output:
 *      None
**/
void freeNode(Node* node){
    
    //if it is not a leaf
    if(node->numOfIndexes==1){

        //free its left child
        if(node->left!=NULL){
            freeNode(node->left);
            free(node->left);
        }

        //free its right child
        if(node->right != NULL){
            freeNode(node->right);
            free(node->right);
        }

        free(node->dists);
        free(node->indx);
    }

    //if it is not a leaf
    if((node->left == NULL) && (node->right == NULL)){
        free(node->dists);
        free(node->indx);
    }

    return;
}

#endif