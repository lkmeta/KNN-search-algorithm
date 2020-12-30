#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

typedef struct queryPoint queryPoint;
typedef struct Node Node;

// Definition of the query point struct
struct queryPoint
{
    double *coord; //!< d coords for query point               [1-by-d]
    int d;
    int *nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double *ndist; //!< Distance of nearest neighbors          [m-by-k]
    int k;         //!< Number of nearest neighbors            [scalar]
    int numOfIndexes;
    int flag; //!< 0 if we haven't found yet the nearest neighbors else 1
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
        while (j >= 0 && arr[j] - key > 0.001f)
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

    srand(time(NULL));

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

        // int pInd = findVPIndx(nd->p, indexes, n);
        // // printf("Index of nd->p=%d\n",pInd);
        // // printf("nd->p=\n");
        // // for(int i=0;i<d;++i){
        // //     printf("%lf ",S[pInd*d+i]);
        // // }
        // // printf("\n");

        // for (int i = 0; i < n; ++i)
        // {
        //     nd->dists[i] = 0;
        //     for (int j = 0; j < d; ++j)
        //     {
        //         nd->dists[i] += pow(S[i * d + j] - S[pInd * d + j], 2);
        //     }
        //     nd->dists[i] = sqrt(nd->dists[i]);
        // }

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
            if (nd->dists[i] <= 0.0001f)
            {
                if (i == n - 1)
                    break;
                else
                    continue;
            }
            if (nd->dists[i] - nd->mu <= 0.0001f)
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
            if (nd->dists[i] <= 0.0001f)
            {
                if (i == n - 1)
                    break;
                else
                    continue;
            }
            if (nd->dists[i] - nd->mu <= 0.0001f)
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

void addElements(Node *currentNode, queryPoint *currentP, double *S)
{
    // Check if Node is empty
    if (currentNode->numOfIndexes == 0)
    {
        printf("Error in addElements: No idx.");
        return;
    }

    int length = currentNode->numOfIndexes;
    printf("\nWORKING");
    printf("\nlength = %d ", length);

    int *tempIdx = (int *)malloc(length * sizeof(int));
    if (tempIdx == NULL)
    {
        printf("Error in addElements: Couldn't allocate memory for tempIdx");
        exit(-1);
    }

    double *tempDist = (double *)malloc(length * sizeof(double));
    if (tempDist == NULL)
    {
        printf("Error in addElements: Couldn't allocate memory for tempDist");
        exit(-1);
    }

    printf("\ncurrentNode: p = %d , mu = %lf ", currentNode->p, currentNode->mu);
    // printf("\nnumOfInd: %d ", currentNode->numOfIndexes);
    // printf("\nIND: ");
    // for (int i = 0; i < currentNode->numOfIndexes; i++)
    // {
    //     printf(" %d ", currentNode->indx[i]);
    // }

    for (int i = 0; i < length; i++)
    {   
        // temporary idx to add
        tempIdx[i] = currentNode->indx[i];

        // eimai se node
        if (currentNode->p != -1)
        {
            // calculate temporary dist to add
            tempDist[i] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                tempDist[i] += pow(S[2 * currentNode->p + j] - currentP->coord[j], 2);
            }
            tempDist[i] = sqrt(tempDist[i]);
        }
        // eimai se fyllo
        else
        {
            // calculate temporary dist to add
            tempDist[i] = 0;
            for (int j = 0; j < currentP->d; ++j)
            {
                tempDist[i] += pow(S[2 * currentNode->indx[i] + j] - currentP->coord[j], 2);
            }
            tempDist[i] = sqrt(tempDist[i]);
        }
    }

    // printf("\ntempIdx = ");
    // for (int i = 0; i < currentNode->numOfIndexes; i++)
    // {
    //     printf(" %d ", tempIdx[i]);
    // }
    // printf("\ntempDist = ");
    // for (int i = 0; i < currentNode->numOfIndexes; i++)
    // {
    //     printf(" %lf ", tempDist[i]);
    // }

    insertionSort(tempDist, tempIdx, length);

    printf("\nAFTER insertionSort");
    printf("\ntempIdx = ");
    for (int i = 0; i < currentNode->numOfIndexes; i++)
    {
        printf(" %d ", tempIdx[i]);
    }
    printf("\ntempDist = ");
    for (int i = 0; i < currentNode->numOfIndexes; i++)
    {
        printf(" %lf ", tempDist[i]);
    }
    printf("\nnumOfInd = %d ", currentP->numOfIndexes);

    // Final lists with indexes and dists for this point
    int *addedIdx = (int *)malloc(length * sizeof(int));
    if (addedIdx == NULL)
    {
        printf("Error in addElements: Couldn't allocate memory for addedIdx");
        exit(-1);
    }

    double *addedDists = (double *)malloc(length * sizeof(double));
    if (addedDists == NULL)
    {
        printf("Error in addElements: Couldn't allocate memory for addedDists");
        exit(-1);
    }

    for (int i = 0; i < length; i++)
    {
        // elegxos an phra hdh ta k pio kontina
        if (currentP->k == currentP->numOfIndexes || currentP->flag == 1)
        {
            printf("\nk nearest neighbors are found");
            currentP->flag = 1;
            break;
        }
        else if (currentP->k - currentP->numOfIndexes > 0)
        {
            addedIdx[i] = tempIdx[i];
            addedDists[i] = tempDist[i];
            currentP->numOfIndexes++;
        }
    }

    printf("\nFinally");
    printf("\ntempIdx = ");
    for (int i = 0; i < currentNode->numOfIndexes; i++)
    {
        printf(" %d ", addedIdx[i]);
    }
    printf("\ntempDist = ");
    for (int i = 0; i < currentNode->numOfIndexes; i++)
    {
        printf(" %lf ", addedDists[i]);
    }
    printf("\nnumOfInd = %d ", currentP->numOfIndexes);

    /**
     * TODO: 
     *      -add tis final lists sto currentP kai na krathsw ta k pio kontina
     * 
    */
}

queryPoint *searchVPT(Node *root, queryPoint *p, double *S)
{
    printf("... inside searchVPT ...\n");

    Node *currentNode = root;
    if (currentNode == NULL)
    {
        printf("Error in searchVPT: Couldn't allocate memory for currentNode");
        exit(-1);
    }

    // printf("currentNode: p = %d , mu = %lf ", currentNode->p, currentNode->mu);
    // printf("\nnumOfInd: %d ", currentNode->numOfIndexes);
    // printf("\nIND: ");
    // for (int i = 0; i < currentNode->numOfIndexes; i++)
    // {
    //     printf(" %d ", currentNode->indx[i]);
    // }

    queryPoint *currentP = p;
    if (currentP == NULL)
    {
        printf("Error in searchVPT: Couldn't allocate memory for currentP");
        exit(-1);
    }

    // Push nearest neighbors in query point
    addElements(currentNode, currentP, S);

    if (currentNode->right == NULL && currentNode->left == NULL)
    {
        return currentP;
    }
    else
    {

        double tempDist = 0;

        // Calculating distance between queryPoint and vp
        for (int i = 0; i < p->d; ++i)
        {
            tempDist += pow(S[2 * currentNode->p + i] - p->coord[i], 2);
        }
        tempDist = sqrt(tempDist);

        printf("\ntempDist = %lf ", tempDist);

        // Check for searching left or right here
        if (tempDist < currentNode->mu)
        {
            printf("\n\nSearching on left child.\n");
            currentP = searchVPT(currentNode->left, currentP, S);

            // Check intersection
            if (currentP->flag == 0)
            {
                printf("\n\nintersection: Searching on right child.\n");
                currentP = searchVPT(currentNode->right, currentP, S);
            }
        }
        else
        {
            printf("\n\nSearching on right child.\n");
            currentP = searchVPT(currentNode->right, currentP, S);

            // Check intersection
            if (currentP->flag == 0)
            {
                printf("\n\nintersection: Searching on left child.\n");
                currentP = searchVPT(currentNode->left, currentP, S);
            }
        }
        // Done with searching VPT
    }

    return currentP;
}

int main()
{

    int n = 10;
    int d = 2;

    double *X = (double *)malloc(n * d * sizeof(double));
    if (X == NULL)
    {
        printf("Error in main: Couldn't allocate memory for X");
        exit(-1);
    }
    //createRandomMatrix(X,n*d);
    X[0] = 1.0;
    X[1] = 3.0;
    X[2] = 4.0;
    X[3] = 2.0;
    X[4] = -2.0;
    X[5] = 0.0;
    X[6] = 5.0;
    X[7] = -1.0;
    X[8] = 7.0;
    X[9] = -3.0;
    X[10] = -4.0;
    X[11] = 6.0;
    X[12] = -8.0;
    X[13] = 1.0;
    X[14] = 3.0;
    X[15] = 4.0;
    X[16] = -7.0;
    X[17] = 5.0;
    X[18] = -1.0;
    X[19] = -2.0;
    // X[20]=-6.0;
    // X[21]=-4.0;

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

    /*
    int sampleSize = 5;
    double* D = (double *)malloc(sampleSize*d*sizeof(double));
    if(D==NULL){
        printf("Error in main: Couldn't allocate memory for D");
        exit(-1);        
    }
    int* sampleIndex = sampleSet(X,D,sampleSize,n,d);
    printf("Finally:\n D=\n");
    for(int k=0;k<sampleSize*d;++k){
            printf("%lf ",D[k]);
        }
    printf("\n sampleIndex=\n");
    for(int k=0;k<sampleSize;++k){
            printf("%d ",sampleIndex[k]);
        }
    int vp =selectVP(X,n,d);
    printf("vp=%d\n",vp);
    */

    int B = 2;

    // VPT Process
    Node *root = makeVPT(X, n, d, indexes, B);
    if (root == NULL)
    {
        printf("Error in main: Couldn't allocate memory for root");
        exit(-1);
    }

    printf("VPT has been created.\n");
    //printf("root: p = %d , mu = %lf ", root->p, root->mu);

    // Search VPT Process
    queryPoint *p;
    p = malloc(sizeof(queryPoint));
    if (p == NULL)
    {
        printf("failed to allocate memory for query point\n");
        return 0;
    }

    double *tempCoord = (double *)malloc(d * sizeof(double));
    for (int i = 0; i < d; i++)
    {
        tempCoord[i] = X[i + 2];
    }

    p->coord = tempCoord;
    p->d = 2;
    p->k = 5;
    p->flag = 0;

    printf("Searching the folowwing query point:\n ");
    printf("x = %lf ", p->coord[0]);
    printf(", y = %lf \n\n", p->coord[1]);
    searchVPT(root, p, X);

    free(indexes);
    free(root);
    free(p);
    free(X);
}
