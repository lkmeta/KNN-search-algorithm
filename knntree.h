#ifndef KNNTREE_H
#define KNNTREE_H

/**
 * TODO: ~mergelist -> yes or no about the instertionSort? EGW LEW NA BALOUME
 *       ~queryPoint -> yes or no about the flag variable? EGW LEW NA FYGEI
 * 
 * add more COMMENTS: sampledAlready, sampleSet
 * 
*/

#include <stdio.h>
#include <stdlib.h>

#include <math.h>
#include <time.h>
#include <string.h>

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

typedef struct knnresult knnresult;
typedef struct queryPoint queryPoint;
typedef struct Node Node;

// Definition of the kNN result struct
struct knnresult
{
    int *nidx;     //!< Indices (0-based) of nearest neighbors      [m-by-k]
    double *ndist; //!< Distance of nearest neighbors               [m-by-k]
    int m;         //!< Number of query points                      [scalar]
    int k;         //!< Number of nearest neighbors                 [scalar]
};

// Definition of the query point struct
struct queryPoint
{
    double *coord;    //!< d coords for query point                    [1-by-d]
    int d;            //!< Number of coords                            [scalar]
    int *nidx;        //!< Indices (0-based) of nearest neighbors      [m-by-k]
    double *ndist;    //!< Distance of nearest neighbors               [m-by-k]
    int k;            //!< Number of nearest neighbors                 [scalar]
    int numOfIndexes; //!< Counter for Indeces in a queryPoint         [scalar]
    int flag;         //!< Counter for the visited nodes from the VPT  [scalar]
    double tau;       //!< Radius for the searchVPT process            [scalar]
};

// Definition of the Node struct
struct Node
{
    Node *left;       //!< Left child Node                             [ Node ]
    Node *right;      //!< Right child Node                            [ Node ]
    int p;            //!< Vantage Point for the Node                  [scalar]
    double mu;        //!< Median value for the Node                   [scalar]
    double *dists;    //!< Distance of nearest neighbors               [m-by-k]
    int *indx;        //!< Indices (0-based) of nearest neighbors      [m-by-k]
    int numOfIndexes; //!< Counter for Indeces in a Node               [scalar]
};

//function printing a matrix of doubles
void printMatrix(double *A, int size)
{
    for (int i = 0; i < size; ++i)
    {
        printf("%lf ", A[i]);
    }
    printf("\n");
}

//function creating a matrix of doubles with random values in [-100,100]
double createRandomMatrix(double *A, int size)
{
    for (int i = 0; i < size; ++i)
    {
        A[i] = (double)rand() / RAND_MAX * 2.0 - 100.0;
    }
}

/** 
* Partition using Lomuto partition scheme
* We do the partition applying comparisons to the elements of A but we move the elements of B in the same way
* so that they are rearranged in correspodence those of A.
* O(n) complexity
* Input:
*       double* A: first matrix in which we apply partition. Its elements are used inthe comparisons with the pivot
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
        //we allow for our doubles a tolerance error of 0.001 (not every number can be stored exactly as a binary)
        if (A[i] - pivot < 1e-6)
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

/* Function to sort an array using insertion sort */
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
        while (j >= 0 && arr[j] - key > 1e-7)
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
 * Function that calculates if the sampled already.
 * If already sampled then returns 1 else 0.
 * Input:
 *      double* X      : nxd matrix containing the corpus data points (n points with d coordinates each)
 *      double *sample : list with sample
 *      int *indexes   : list with indexes for sample 
 *      int sampleSize : num of samples in the list
 *      int n          : number of corpus points
 *      int d          : number of dimensions
 * Output:
 *      int *sampleIndex        : list with sample indexes
**/
int sampledAlready(int *sampleIndex, int sampleSize, int index)
{
    int result = 0;

    for (int i = 0; i < sampleSize; ++i)
    {
        if (sampleIndex[i] == -1)
        {
            break;
        }
        if (index == sampleIndex[i])
        {
            result = 1;
            break;
        }
    }

    return result;
}

/**
 * Function that calculates samplesSet indexes list.
 * Input:
 *      double* X        : nxd matrix containing the corpus data points (n points with d coordinates each)
 *      double *sample   : list with sample
 *      int *indexes     : list with indexes for sample 
 *      int sampleSize   : num of samples in the list
 *      int n            : number of corpus points
 *      int d            : number of dimensions
 * Output:
 *      int *sampleIndex : list with sample indexes
**/
int *sampleSet(double *X, double *sample, int *indexes, int sampleSize, int n, int d)
{
    int count = 0;
    int index;

    int *sampleIndex = (int *)malloc(sampleSize * sizeof(int));
    if (sampleIndex == NULL)
    {
        printf("Error in sampleSet: Couldn't allocate memory for sampleIndex");
        exit(-1);
    }

    for (int i = 0; i < sampleSize; ++i)
    {
        sampleIndex[i] = -1;
    }

    while (count < sampleSize)
    {
        index = rand() % n;

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
 * Function that calculates median from a list.
 * O(nlogn) time complexity.
 * Input:
 *      double *sampleDistances : list with distances
 *      int *indexes            : list with indexes for sampleDistances
 *      int sampleSize          : num of samples in sampleDistances
 * Output:
 *      double mu               : median value
**/
double findMedian(double *sampleDistances, int *indexes, int sampleSize)
{
    double mu;

    //first we sort the array
    insertionSort(sampleDistances, indexes, sampleSize);

    int middle = (sampleSize + 1) / 2 - 1;

    //check for even or odd case
    if (sampleSize % 2 == 1)
    {
        mu = sampleDistances[middle];
    }
    else
    {
        mu = (sampleDistances[middle] + sampleDistances[middle + 1]) / 2;
    }

    return mu;
}

/**
 * Function that selects the vantage point.
 * Input: 
 *      double* X    : nxd matrix containing the corpus data points (n points with d coordinates each)
 *      int *indexes : list with indexes for sampleDistances
 *      int n        : number of corpus points
 *      int d        : number of dimensions
 * Output:
 *      int bestP    : best vantage point value
**/
int selectVP(double *X, int *indexes, int n, int d)
{
    //pick a random sampleSize
    int sampleSize = 5;

    while (sampleSize >= n)
    {
        sampleSize /= 2;
    }

    double *P = (double *)malloc(sampleSize * d * sizeof(double));
    if (P == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for P");
        exit(-1);
    }

    double *D = (double *)malloc(sampleSize * d * sizeof(double));
    if (D == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for D");
        exit(-1);
    }

    double *sampleDistances = (double *)malloc(sampleSize * sizeof(double));
    if (sampleDistances == NULL)
    {
        printf("Error in selectVP: Couldn't allocate memory for sampleDistances");
        exit(-1);
    }

    int *sampleIndex = sampleSet(X, P, indexes, sampleSize, n, d);

    int middle;
    double mu;
    double bestSpread = 0.0;
    double spread;
    int bestP;

    for (int i = 0; i < sampleSize; ++i)
    {
        int *sampleIndex2 = sampleSet(X, D, indexes, sampleSize, n, d);

        for (int j = 0; j < sampleSize; ++j)
        {
            sampleDistances[j] = 0;
            for (int k = 0; k < d; ++k)
            {
                sampleDistances[j] += pow(P[i * d + k] - D[j * d + k], 2);
            }
            sampleDistances[j] = sqrt(sampleDistances[j]);
        }

        mu = findMedian(sampleDistances, sampleIndex2, sampleSize);
        free(sampleIndex2);

        spread = 0;
        for (int j = 0; j < sampleSize; ++j)
        {
            spread += pow(sampleDistances[j] - mu, 2);
        }
        spread /= sampleSize;

        if (spread > bestSpread)
        {
            bestSpread = spread;
            bestP = sampleIndex[i];
        }
    }

    //deallocate memory
    free(P);
    free(D);
    free(sampleDistances);
    free(sampleIndex);

    return bestP;
}

/**
 * Function that finds the index from the vantage point.
 * O(n) complexity
 * Input:
 *      int vp       : value of the vantage point
 *      int *indexes : list with n indexes
 * Output:
 *      int VPindex  : index of the vp in the list
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
* Function that creates the vantage point tree. 
* Firstly, we create the root node with the vantage point, median, indexes and distances array
* and then we continue similarly into left and right nodes.
* The process ends when "n < 2 * B + 1" which means that we are in a leaf node.
* Input:
*       double* S    : nxd matrix containing the corpus data points (n points with d coordinates each)
*       int n        : number of corpus set X points
*       int d        : number of dimensions
*       int *indexes : list with n indexes
*       int B        : the least number of points leaves nodes contains
* Output:
*       Node *nd     : root node from the Vantage Point Tree
**/
Node *makeVPT(double *S, int n, int d, int *indexes, int B)
{
    //allocate memory for root node
    Node *nd;

    //check if array has 0 indexes
    if (n == 0)
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

    if (n < 2 * B + 1)
    {
        //leaf node case

        nd->p = -1;
        nd->mu = -1.0;
        nd->left = NULL;
        nd->right = NULL;
        nd->dists = (double *)malloc(n * sizeof(double));
        nd->indx = (int *)malloc(n * sizeof(int));
        nd->numOfIndexes = n; //add length of indx array
        if (nd->indx == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->indx");
            exit(-1);
        }

        memcpy(nd->indx, indexes, n * sizeof(int));
    }
    else
    {
        //node case

        //choose vantage point for the node
        nd->p = selectVP(S, indexes, n, d);

        //find index for the choosen vp
        int pInd = findVPIndx(nd->p, indexes, n);

        //allocate memory for the distances and calculate them
        nd->dists = (double *)malloc(n * sizeof(double));
        if (nd->dists == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for nd->dists");
            exit(-1);
        }

        for (int i = 0; i < n; ++i)
        {
            nd->dists[i] = 0;
            for (int j = 0; j < d; ++j)
            {
                nd->dists[i] += pow(S[i * d + j] - S[pInd * d + j], 2);
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

        //add length of indx array for current node
        nd->numOfIndexes = 1;

        //use temporary dists list to find the median value
        double *distsTemp = (double *)malloc(n * sizeof(double));
        if (distsTemp == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for distsTemp");
            exit(-1);
        }

        memcpy(distsTemp, nd->dists, n * sizeof(double));

        nd->mu = findMedian(distsTemp, indexes, n);

        free(distsTemp);

        //child nodes process
        double *right;
        int *rightIndexes;
        double *left;
        int *leftIndexes;

        int lSize = 0;
        int rSize = 0;

        for (int i = 0; i < n; ++i)
        {
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

        int lCounter = 0;
        int rCounter = 0;

        for (int i = 0; i < n; ++i)
        {
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

        nd->left = makeVPT(left, lSize, d, leftIndexes, B);
        nd->right = makeVPT(right, rSize, d, rightIndexes, B);
    }

    return nd;
}

/**
 * Function that adds the neighbors from the node into the query point.
 * Firstly, we add all the neighbors from current node into the query point.
 * Then we sort the list with the distances. Finally, we reallocate the list 
 * and we keep the nearest neighbors.
 * 
 * Input:
 *      Node *currentNode    : current node in the tree
 *      queryPoint *currentP : current query point we are looking for the k nearest neighbors 
 *      double *S            : nxd matrix containing the corpus data points (n points with d coordinates each)
 * Output:
 *      None
**/
void addElements(Node *currentNode, queryPoint *currentP, double *S)
{

    //used lists in this funtion with indexes and dists from cuppent point
    // currentP->nidx;
    // currentP->ndist;

    //check if Node is empty
    if (currentNode->numOfIndexes == 0)
    {
        printf("Error in addElements: No indexes.");
        return;
    }

    int length = currentNode->numOfIndexes;          // current length
    int newLength = length + currentP->numOfIndexes; // possible next length

    //reallocate memory cause of new possible length
    currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * newLength);
    currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * newLength);

    //add new indexes and distances from cureent node into the query point
    for (int i = 0; i < length; i++)
    {
        //node case
        if (currentNode->p != -1)
        {
            //temporary indx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->p;

            //calculate temporary dist to add
            currentP->ndist[i + currentP->numOfIndexes] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                currentP->ndist[i + currentP->numOfIndexes] += pow(S[currentP->d * currentNode->p + j] - currentP->coord[j], 2);
            }
            currentP->ndist[i + currentP->numOfIndexes] = sqrt(currentP->ndist[i + currentP->numOfIndexes]);
        }
        //leaf case
        else
        {
            //temporary indx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->indx[i];

            //calculate temporary dist to add
            currentP->ndist[i + currentP->numOfIndexes] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                currentP->ndist[i + currentP->numOfIndexes] += pow(S[currentP->d * currentNode->indx[i] + j] - currentP->coord[j], 2);
            }
            currentP->ndist[i + currentP->numOfIndexes] = sqrt(currentP->ndist[i + currentP->numOfIndexes]);
        }
    }

    //sort the distances and and indexes
    if (newLength != 1)
    {
        insertionSort(currentP->ndist, currentP->nidx, newLength);
    }

    //find the real number of indexes that will keep
    for (int i = 0; i < length; i++)
    {
        if (currentP->k - currentP->numOfIndexes > 0)
        {
            currentP->numOfIndexes++;
        }
    }

    //reallocate the length from distances and indexes lists
    if (currentP->numOfIndexes != newLength)
    {
        currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * (currentP->numOfIndexes));
        currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * (currentP->numOfIndexes));
    }

    //change tau
    if (currentP->k == currentP->numOfIndexes)
    {
        currentP->tau = currentP->ndist[currentP->numOfIndexes - 1];
    }
}

/**
* Function that searches the k nearest distances for the query point in the VPT.
* Firstly, we check if we are searching in a leaf node or not. 
* When we are in a leaf node we push the nearest neighbors in query point and we return it.
* When we are not in a leaf node, we check for the following possible cases:
*   1. if tempdist is lower than the current radius we use addElements
*   2. if tempdist is lower than the current median we search in left child
*   3. if tempdist is higher than the current median we search in right child
* Finally, we check for intersection case and then we return the currentP with all the k nearest neighbors. 
* Input:
*       Node *root         : root node from the VPT
*       queryPoint *queryP : query point that we are looking for the k nearest neighbors 
*       double *S          : nxd matrix containing the corpus data points (n points with d coordinates each)
* Output:
*       queryPoint *queryP : query point with the k nearest neighbors 
**/
queryPoint *searchVPT(Node *root, queryPoint *queryP, double *S)
{
    //check if root or query point are null
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

    // Increase the counter for the visited nodes
    queryP->flag++;

    if (root->right == NULL && root->left == NULL)
    {
        //leaf node case

        //push nearest neighbors in query point
        addElements(root, queryP, S);

        return queryP;
    }
    else
    {
        //node case

        double tempDist = 0;

        //calculating distance between queryPoint and vp
        for (int i = 0; i < queryP->d; ++i)
        {
            tempDist += pow(S[queryP->d * root->p + i] - queryP->coord[i], 2);
        }
        tempDist = sqrt(tempDist);

        if (tempDist < queryP->tau)
        {
            //push nearest neighbors in query point
            addElements(root, queryP, S);
        }

        if (tempDist < root->mu)
        {
            //search in left child process
            queryP = searchVPT(root->left, queryP, S);

            // Check intersection
            if (queryP->k > queryP->numOfIndexes || tempDist > root->mu - queryP->tau)
            {
                //search in right child process
                queryP = searchVPT(root->right, queryP, S);
            }
        }

        if (tempDist >= root->mu)
        {
            //search in right child process
            queryP = searchVPT(root->right, queryP, S);

            // Check intersection
            if (queryP->k > queryP->numOfIndexes || tempDist < root->mu + queryP->tau)
            {
                //search in left child process
                queryP = searchVPT(root->left, queryP, S);
            }
        }
    }

    return queryP;
}

/**
 * Function that implements the knn algorithm for a given query set Y and corpus set X.
 * At first we check if Y is the same as X, after we create the Vantage Point Tree
 * and then for each point in Y we search in the VPT for the k nearest neighbors
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

    //allocate memory for nidx and ndist lists
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

    //allocate memory indexes number
    int *indexes = (int *)malloc(n * sizeof(int));
    if (indexes == NULL)
    {
        printf("Error in main: Couldn't allocate memory for indexes");
        exit(-1);
    }

    for (int i = 0; i < n; ++i)
    {
        indexes[i] = i;
    }

    //choose B for the least number of points leaves nodes contains
    int B = 1;

    //create VPT Process
    Node *root = makeVPT(X, n, d, indexes, B);
    if (root == NULL)
    {
        printf("Error in main: Couldn't allocate memory for root");
        exit(-1);
    }

    //search VPT Process
    //allocate memory for one query point
    queryPoint *p;

    //for each point of the query set Y, find its knn using searchVPT
    for (int j = 0; j < m; j++)
    {
        p = malloc(sizeof(queryPoint));
        if (p == NULL)
        {
            printf("failed to allocate memory for query point\n");
            exit(-1);
        }

        p->d = d;
        p->k = k;
        p->flag = 0;
        p->tau = INFINITY;

        double *tempCoord = (double *)malloc(d * sizeof(double));
        if (tempCoord == NULL)
        {
            printf("failed to allocate memory for tempCoord\n");
            exit(-1);
        }

        for (int i = 0; i < d; i++)
        {
            tempCoord[i] = Y[i + d * j];
        }

        p->coord = tempCoord;
        p->nidx = NULL;
        p->ndist = NULL;
        p->numOfIndexes = 0;

        p = searchVPT(root, p, X);

        for (int i = 0; i < p->numOfIndexes; i++)
        {
            nidx[j * k + i] = p->nidx[i];
            ndist[j * k + i] = p->ndist[i];
        }

        free(tempCoord);
    }

    //struct to be returned
    knnresult retVal;

    retVal.k = k;
    retVal.m = m;
    retVal.ndist = ndist;
    retVal.nidx = nidx;

    free(p);
    free(indexes);
    free(root);

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
        quickSelect(ndistComb, nidxComb, 0, 2 * k - 1, k - 1);

        insertionSort(ndistComb, nidxComb, k);

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

#endif
