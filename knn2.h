#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#define SWAPD(x,y) { double temp = x; x = y; y = temp; }
#define SWAPI(x,y) { int temp = x; x = y; y =temp; }

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

struct Node
{
    Node *left;
    Node *right;
    int p;
    double mu;
    double *dists;
    int *indx;
    int numOfIndexes;
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

// // Returns the k-th smallest element of list within left..right
// // (i.e. left <= k <= right). The search space within the array is
// // changing for each round - but the list is still the same size.
// // Thus, k does not need to be updated with each round.
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

/* Function to sort an array using insertion sort*/
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

//0 an den to vrei, 1 an to vrei
int sampledAlready(int *sampleIndex, int sampleSize, int index)
{
    int result = 0;
    /*
    printf("sampleIndex=\n");
    for(int i=0;i<sampleSize;++i){
        printf("%d ",sampleIndex[i]);
    }
*/
    for (int i = 0; i < sampleSize; ++i)
    {
        if (sampleIndex[i] == -1)
        {
            break;
        }
        if (index == sampleIndex[i])
        {
            //printf("Found %d it in position %d\n",index,i);
            result = 1;
            break;
        }
    }

    return result;
}

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
    int i = 0;
    while (count < sampleSize)
    {
        //printf("\ni=%d, count=%d\n", i,count);
        index = rand() % n;
        //printf("index=%d\n", index);

        if (sampledAlready(sampleIndex, sampleSize, indexes[index]) == 0)
        {
            //printf("Didn't find it\n");
            sampleIndex[count] = indexes[index];
            for (int j = 0; j < d; ++j)
            {
                sample[count * d + j] = X[index * d + j];
            }
            count++;
        }

        // printf("sampleIndex=\n");
        // for(int k=0;k<sampleSize;++k){
        //     printf("%d ",sampleIndex[k]);
        // }
        // printf("\n");

        i++;
    }

    return sampleIndex;
}

double findMedian(double *sampleDistances, int *indexes, int sampleSize)
{
    double mu;

    // printf("\nbefore insertionSort nd->indx=\n");
    // for(int i=0;i<sampleSize;++i){
    //     printf("%d ", indexes[i]);
    // }

    insertionSort(sampleDistances, indexes, sampleSize);
    /*
    printf("\nSorted sampleDistances=\n");
    for(int j=0;j<sampleSize;++j){
        printf("%lf ",sampleDistances[j]);
    }
*/

    // printf("\nafter insertionSort nd->indx=\n");
    // for(int i=0;i<sampleSize;++i){
    //     printf("%d ", indexes[i]);
    // }

    int middle = (sampleSize + 1) / 2 - 1;

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

int selectVP(double *X, int *indexes, int n, int d)
{

    // printf("S=\n");
    // for(int i=0;i<n*d;++i){
    //     printf("%lf ", X[i]);
    // }

    int sampleSize = 5; //na dw an yparxei isws allh kalyterh epilogh

    while (sampleSize >= n)
    {
        sampleSize /= 2;
    }

    //    printf("samplesize=%d\n", sampleSize);

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

    // printf("\nP=\n");
    // for(int i=0;i<sampleSize*d;++i){
    //     printf("%lf ",P[i]);
    // }

    int middle;
    double mu;
    double bestSpread = 0.0;
    double spread;
    int bestP;

    for (int i = 0; i < sampleSize; ++i)
    {
        //        printf("\ni=%d, p=%lf %lf\n",i,P[i*d], P[i*d+1]);
        int *sampleIndex2 = sampleSet(X, D, indexes, sampleSize, n, d);
        /*        printf("\nD=\n");
        for(int j=0;j<sampleSize*d;++j){
            printf("%lf ", D[j]);
        }
*/
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
        /*        printf("\nsampleDistances=\n");
        for(int j=0;j<sampleSize;++j){
            printf("%lf ",sampleDistances[j]);
        }
*/
        mu = findMedian(sampleDistances, sampleIndex2, sampleSize);
        free(sampleIndex2);
        //        printf("mu=%lf\n", mu);
        spread = 0;
        for (int j = 0; j < sampleSize; ++j)
        {
            spread += pow(sampleDistances[j] - mu, 2);
        }
        spread /= sampleSize;
        //        printf("spread=%lf\n",spread);
        if (spread > bestSpread)
        {
            bestSpread = spread;
            bestP = sampleIndex[i];
        }
    }

    //printf("bestspread=%lf, bestP=%d\n",bestSpread,bestP);

    free(P);
    free(D);
    free(sampleDistances);
    free(sampleIndex);

    return bestP;
}

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

Node *makeVPT(double *S, int n, int d, int *indexes, int B)
{
    Node *nd;

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
        nd->p = -1;
        nd->mu = -1.0;
        nd->left = NULL;
        nd->right = NULL;
        nd->dists = (double *)malloc(n * sizeof(double));
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

    else
    {
        nd->p = selectVP(S, indexes, n, d); //to index tou kalyterou stoixeio ston arxiko pinaka X

        // printf("\np = %d", nd->p);
        // printf("\nIndexes= ");
        // for (int i = 0; i < n; ++i)
        // {
        //     printf("%d ", indexes[i]);
        // }
        // printf("\n");

        int pInd = findVPIndx(nd->p, indexes, n); //to index tou vp ston pinaka S ths sugkekrimenhs klhshs

        // printf("Index of nd->p=%d\n",pInd);
        // printf("nd->p=\n");
        // for(int i=0;i<d;++i){
        //     printf("%lf ",S[pInd*d+i]);
        // }
        // printf("\n");

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

        // add length of indx array
        nd->numOfIndexes = 1;

        // printf("before median nd->dists=\n");
        // for(int i=0;i<n;++i){
        //     printf("%lf ", nd->dists[i]);
        // }

        // printf("\nbefore median nd->indx=\n");
        // for(int i=0;i<n;++i){
        //     printf("%d ", nd->indx[i]);
        // }

        double *distsTemp = (double *)malloc(n * sizeof(double)); //giati h findMedian xalaei thn seira twn stoixeiwn
        if (distsTemp == NULL)
        {
            printf("Error in makeVPT: Couldn't allocate memory for distsTemp");
            exit(-1);
        }

        memcpy(distsTemp, nd->dists, n * sizeof(double));

        nd->mu = findMedian(distsTemp, indexes, n);

        free(distsTemp);
        //        free(indexes); //isws na kanei to free auth pou thn pernaei san orisma

        // printf("\nmu=%lf\n",nd->mu);

        // printf("after median nd->dists=\n");
        // for(int i=0;i<n;++i){
        //     printf("%lf ", nd->dists[i]);
        // }

        // printf("\nafter median nd->indx=\n");
        // for(int i=0;i<n;++i){
        //     printf("%d ", nd->indx[i]);
        // }

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

        //printf("lSize=%d, rSize=%d",lSize,rSize);

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
            { //<=
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

        // printf("\nleftIndexes=\n");
        // for(int i=0;i<lSize;++i){
        //     printf("%d ", leftIndexes[i]);
        // }

        // printf("\nleft=\n");
        // for(int i=0;i<lSize*d;++i){
        //     printf("%lf ", left[i]);
        // }

        // printf("\nrightIndexes=\n");
        // for(int i=0;i<rSize;++i){
        //     printf("%d ", rightIndexes[i]);
        // }

        // printf("\nright=\n");
        // for(int i=0;i<rSize*d;++i){
        //     printf("%lf ", right[i]);
        // }

        nd->left = makeVPT(left, lSize, d, leftIndexes, B);
        nd->right = makeVPT(right, rSize, d, rightIndexes, B);
    }

    return nd;
}

void addElements(Node *currentNode, queryPoint *currentP, double *S, int offset)
{

    // Lists with indexes and dists from cuppent point
    // currentP->nidx;
    // currentP->ndist;

    // Check if Node is empty
    if (currentNode->numOfIndexes == 0)
    {
        printf("Error in addElements: No idx.");
        return;
    }

    int length = currentNode->numOfIndexes;          // current length
    int newLength = length + currentP->numOfIndexes; // possible next length

    /* PRINT OUT THE INPUTS AND THEIR VALUES WHILE STARTING THE addElements Process */
    // printf("\nStart addElement Process. ");
    // printf("\ncurrentNode: p = %d , mu = %lf ", currentNode->p, currentNode->mu);
    // printf("\nlen = %d ", length);
    // printf(", numOfIndexes = %d ", (currentP->numOfIndexes));
    // printf(", NEW LEN = %d ", newLength);
    // printf("\nBEFORE the addElements Process\nInd = ");
    // for (int i = 0; i < (currentP->numOfIndexes); i++)
    // {
    //     printf(" %d ", currentP->nidx[i]);
    // }
    // printf("\nDIST = ");
    // for (int i = 0; i < (currentP->numOfIndexes); i++)
    // {
    //     printf(" %lf ", currentP->ndist[i]);
    // }

    currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * newLength);
    currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * newLength);

    for (int i = 0; i < length; i++)
    {
        // eimai se node
        if (currentNode->p != -1)
        {
            // temporary idx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->p;

            // calculate temporary dist to add
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
        // eimai se fyllo
        else
        {
            // temporary idx to add
            currentP->nidx[i + currentP->numOfIndexes] = currentNode->indx[i];

            // calculate temporary dist to add
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

    if (newLength != 1)
    {
        insertionSort(currentP->ndist, currentP->nidx, newLength);
    }

    /* PRINT OUT RESULTS AFTER insertionSort */
    // printf("\nAFTER insertionSort");
    // printf("\nIDX = ");
    // for (int i = 0; i < newLength; i++)
    // {
    //     printf(" %d ", currentP->nidx[i]);
    // }
    // printf("\nDIST = ");
    // for (int i = 0; i < newLength; i++)
    // {
    //     printf(" %lf ", currentP->ndist[i]);
    // }

    for (int i = 0; i < length; i++)
    {
        if (currentP->k - currentP->numOfIndexes > 0)
        {
            currentP->numOfIndexes++;
        }
    }

    if (currentP->numOfIndexes != newLength)
    {
        currentP->nidx = (int *)realloc(currentP->nidx, sizeof(int) * (currentP->numOfIndexes));
        currentP->ndist = (double *)realloc(currentP->ndist, sizeof(double) * (currentP->numOfIndexes));
    }

    if (currentP->k == currentP->numOfIndexes)
    {
        currentP->tau = currentP->ndist[currentP->numOfIndexes - 1];
    }

    /* PRINT OUT THE RESULTS FROM EVERY ADDING ELEMENT PROCESS */
    // printf("\nNEW NUMOFINDEXES = %d ", currentP->numOfIndexes);
    // printf("\nfinal Ind = ");
    // for (int i = 0; i < currentP->numOfIndexes; i++)
    // {
    //     printf(" %d ", currentP->nidx[i]);
    // }
    // printf("\nfinal Dists = ");
    // for (int i = 0; i < currentP->numOfIndexes; i++)
    // {
    //     printf(" %lf ", currentP->ndist[i]);
    // }
    // printf("\nnumOfInd = %d ", currentP->numOfIndexes);
    // printf("\ntau = %lf ", currentP->tau);
}

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

    /* PRINT OUT THE INPUTS AND THEIR VALUES WHILE STARTING THE searchVPT Process */
    // printf("Start searchVPT process. \n");
    // printf("root: p = %d , mu = %lf ", root->p, root->mu);
    // printf("\nnumOfInd: %d ", root->numOfIndexes);
    // printf("\nIND: ");
    // for (int i = 0; i < root->numOfIndexes; i++)
    // {
    //     printf(" %d ", root->indx[i]);
    // }

    // Increase the counter for the visited nodes
    queryP->flag++;

    // eimai se fyllo
    if (root->right == NULL && root->left == NULL)
    {
        // Push nearest neighbors in query point
        addElements(root, queryP, S, offset);
        return;
    }
    // eimai se node
    else
    {

        double tempDist = 0;

        // Calculating distance between queryPoint and vp
        for (int i = 0; i < queryP->d; ++i)
        {
            tempDist += pow(S[queryP->d * root->p - offset + i] - queryP->coord[i], 2);
        }
        if(tempDist<1e-5){
                tempDist = 0;
            }
        tempDist = sqrt(tempDist);

        // printf("\ntempDist = %lf ", tempDist);

        if (tempDist < queryP->tau)
        {
            // Push nearest neighbors in query point
            addElements(root, queryP, S, offset);
        }

        if (tempDist < root->mu)
        {
            // printf("\n\nSearching on left child.\n");
            searchVPT(root->left, queryP, S, offset);

            // printf("\ntempDist = %lf, mu = %lf, tau = %lf\n", tempDist, root->mu, queryP->tau);
            // Check intersection
            if (queryP->k > queryP->numOfIndexes || tempDist > root->mu - queryP->tau)
            {
                // printf("\n\nintersection: Searching on right child.\n");
                searchVPT(root->right, queryP, S, offset);
            }
        }

        if (tempDist >= root->mu)
        {
            // printf("\n\nSearching on right child.\n");
            searchVPT(root->right, queryP, S, offset);

            // printf("\ntempDist = %lf, mu = %lf, tau = %lf\n", tempDist, root->mu, queryP->tau);
            // Check intersection
            if (queryP->k > queryP->numOfIndexes || tempDist < root->mu + queryP->tau)
            {
                // printf("\n\nintersection: Searching on left child.\n");
                searchVPT(root->left, queryP, S, offset);
            }
        }
        // Done with searching VPT
    }

    return;
}

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

    // search VPT Process

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

    for (int j = 0; j < m; j++)
    {
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

        // printf("\n\nNEW SEARCH\nSearching point: ");
        // printf("x = %lf ", p->coord[0]);
        // printf(", y = %lf \n", p->coord[1]);
        searchVPT(root, p, X, offset);

        //printf("searchVPT Process done.\n");
        //printf("%d nearest neighbors \nIndices: ", p->k);
        for (int i = 0; i < p->numOfIndexes; i++)
        {
            //printf(" %d ", p->nidx[i]);
            nidx[j*k+i] = p->nidx[i];
        }

        //printf("\nDistances: ");
        for (int i = 0; i < p->numOfIndexes; i++)
        {
            //printf(" %lf ", p->ndist[i]);
            ndist[j*k+i] = p->ndist[i];
        }

        //printf("\nFLAG = %d ", p->flag);

        // nidx[j] = p->nidx;
        // ndist[j] = p->ndist;

    }

    free(p->nidx);
    free(p->ndist);
    free(tempCoord);
    free(p);

    // printf("\nNIDX: ");
    // printMatrix((double)nidx, m*k);
    //printf("\nNDIST: ");
    //printMatrix(ndist,m*k);

    //struct to be returned
    knnresult retVal;

    retVal.k = k;
    retVal.m = m;
    retVal.ndist = ndist;
    retVal.nidx = nidx;

    return retVal;
}

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
