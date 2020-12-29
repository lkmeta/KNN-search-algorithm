#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

typedef struct queryPoint queryPoint;
typedef struct Node Node;

// Definition of the query point struct
struct queryPoint
{
    double *coord; //!< d coords for query point               [1-by-d]
    int *nidx;     //!< Indices (0-based) of nearest neighbors [m-by-k]
    double *ndist; //!< Distance of nearest neighbors          [m-by-k]
    int k;         //!< Number of nearest neighbors            [scalar]
};

struct Node
{
    Node *inside;
    Node *outside;
    int p;
    int *nidx;
    double mu;
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

int *sampleSet(double *X, double *sample, int sampleSize, int n, int d)
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

        if (sampledAlready(sampleIndex, sampleSize, index) == 0)
        {
            //printf("Didn't find it\n");
            sampleIndex[count] = index;
            for (int j = 0; j < d; ++j)
            {
                sample[count * d + j] = X[index * d + j];
            }
            count++;
        }
        /*        
        printf("sampleIndex=\n");
        for(int k=0;k<sampleSize;++k){
            printf("%d ",sampleIndex[k]);
        }
        printf("\n");
        i++;
*/
    }

    return sampleIndex;
}

double findMedian(double *sampleDistances, int *indexes, int sampleSize)
{
    double mu;

    insertionSort(sampleDistances, indexes, sampleSize);

    // printf("\nInside findMedian: Sorted sampleDistances=\n");
    // for(int j=0;j<sampleSize;++j){
    //     printf("%lf ",sampleDistances[j]);
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

int selectVP(double *X, int n, int d)
{

    int sampleSize = 10; //na dw an yparxei isws allh kalyterh epilogh

    srand(time(NULL));

    while (sampleSize >= n)
    {
        sampleSize /= 2;
    }

    // printf("samplesize=%d\n", sampleSize);

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

    int *sampleIndex = sampleSet(X, P, sampleSize, n, d);

    // printf("\nP=\n");
    // for(int i=0;i<sampleSize*d;++i){
    //     printf("%lf ",P[i]);
    // }

    int middle;
    double mu;
    double bestSpread = 0.0;
    double spread;
    int bestP = 0;

    for (int i = 0; i < sampleSize; ++i)
    {
        // printf("\ni=%d, p=%lf %lf\n",i,P[i*d], P[i*d+1]);
        int *sampleIndex2 = sampleSet(X, D, sampleSize, n, d);
        // printf("\nD=\n");
        // for(int j=0;j<sampleSize*d;++j){
        //     printf("%lf ", D[j]);
        // }

        for (int j = 0; j < sampleSize; ++j)
        {
            sampleDistances[j] = 0;
            for (int k = 0; k < d; ++k)
            {
                sampleDistances[j] += pow(P[i * d + k] - D[j * d + k], 2);
            }
            sampleDistances[j] = sqrt(sampleDistances[j]);
        }
        // printf("\nsampleDistances=\n");
        // for(int j=0;j<sampleSize;++j){
        //     printf("%lf ",sampleDistances[j]);
        // }

        mu = findMedian(sampleDistances, sampleIndex2, sampleSize);
        free(sampleIndex2);
        //printf("mu=%lf\n", mu);
        spread = 0;
        for (int j = 0; j < sampleSize; ++j)
        {
            spread += pow(sampleDistances[j] - mu, 2);
        }
        spread /= sampleSize;
        //printf("spread=%lf\n",spread);
        if (spread > bestSpread)
        {
            bestSpread = spread;
            bestP = sampleIndex[i];
        }
        //printf("bestspread=%lf, bestP=%d\n",bestSpread,bestP);
    }

    free(P);
    free(D);
    free(sampleDistances);
    free(sampleIndex);

    return bestP;
}

Node *makeVPT(double *S, int n, int d)
{
    Node *nd;

    // printf("\n\nn = %d\nS = ", n);
    // for(int i = 0; i < n*d; i++)
    // {
    //     printf(" %lf ", S[i]);
    // }

    if (n <= 1)
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

    nd->p = selectVP(S, n, d);

    // printf("\nselectVP = %d ", nd->p);

    double *VPDistances = (double *)malloc(n * sizeof(double));
    if (VPDistances == NULL)
    {
        printf("Error in makeVPT: Couldn't allocate memory for VPDistances");
        exit(-1);
    }

    int *VPIndexes = (int *)malloc(n * sizeof(int));
    if (VPIndexes == NULL)
    {
        printf("Error in makeVPT: Couldn't allocate memory for VPIndexes");
        exit(-1);
    }

    for (int i = 0; i < n; ++i)
    {
        VPIndexes[i] = i;
        VPDistances[i] = 0;
        for (int j = 0; j < d; ++j)
        {
            VPDistances[i] += pow(S[i * d + j] - S[(nd->p) * d + j], 2);
        }
        VPDistances[i] = sqrt(VPDistances[i]);
    }

    // printf("\nbefore median VPDistances=\n");
    // for(int i=0;i<n;++i){
    //     printf("%lf ", VPDistances[i]);
    // }

    // printf("\nbefore median VPIndexes=\n");
    // for(int i=0;i<n;++i){
    //     printf("%d ", VPIndexes[i]);
    // }

    nd->mu = findMedian(VPDistances, VPIndexes, n);

    // printf("\nafter median VPDistances=\n");
    // for(int i=0;i<n;++i){
    //     printf("%lf ", VPDistances[i]);
    // }

    // printf("\nafter median VPIndexes=\n");
    // for(int i=0;i<n;++i){
    //     printf("%d ", VPIndexes[i]);
    // }

    /* Make left and right child */
    if (n % 2 == 1)
    {
        double *insideS = malloc((n * d) / 2 * sizeof(double));
        if (!insideS)
        {
            /* handle error */
        }

        double *outsideS = malloc(((n * d) / 2 + 1) * sizeof(double));
        if (!outsideS)
        {
            /* handle error */
        }

        for (int i = 0; i < n / 2; i++)
        {
            for (int j = 0; j < d; j++)
            {
                insideS[i * d + j] = S[2 * VPIndexes[i] + j];
            }
        }

        for (int i = 0; i < n / 2 + 1; i++)
        {
            for (int j = 0; j < d; j++)
            {
                outsideS[i * d + j] = S[2 * VPIndexes[i + n / 2 - 1] + j];
            }
        }

        /* Print out left and right child */
        // printf("\ninsideS: ");
        // for(int i = 0; i < n/2; i++)
        // {
        //     for(int j = 0; j < d; j++)
        //     {
        //         printf(" %lf ", insideS[i*d + j]);
        //     }

        // }

        // printf("\noutsideS: ");
        // for(int i = 0; i < n/2 +1; i++)
        // {
        //     for(int j = 0; j < d; j++)
        //     {
        //         printf(" %lf ", outsideS[i*d + j]);
        //     }

        // }

        nd->inside = makeVPT(insideS, n / 2, d);
        nd->outside = makeVPT(outsideS, n / 2 + 1, d);
    }
    else
    {

        double *insideS = malloc((n * d) / 2 * sizeof(double));
        if (!insideS)
        {
            /* handle error */
        }

        double *outsideS = malloc((n * d) / 2 * sizeof(double));
        if (!outsideS)
        {
            /* handle error */
        }

        for (int i = 0; i < n / 2; i++)
        {
            for (int j = 0; j < d; j++)
            {
                insideS[i * d + j] = S[2 * VPIndexes[i] + j];
            }
        }

        for (int i = 0; i < n / 2; i++)
        {
            for (int j = 0; j < d; j++)
            {
                outsideS[i * d + j] = S[2 * VPIndexes[i + n / 2] + j];
            }
        }
        /* Print out left and right child */
        // printf("\ninsideS: ");
        // for(int i = 0; i < n/2; i++)
        // {
        //     for(int j = 0; j < d; j++)
        //     {
        //         printf(" %lf ", insideS[i*d + j]);
        //     }

        // }

        // printf("\noutsideS: ");
        // for(int i = 0; i < n/2; i++)
        // {
        //     for(int j = 0; j < d; j++)
        //     {
        //         printf(" %lf ", outsideS[i*d + j]);
        //     }
        // }

        nd->inside = makeVPT(insideS, n / 2, d);
        nd->outside = makeVPT(outsideS, n / 2, d);
    }

    return nd;
}

queryPoint *searchVPT(Node *root, queryPoint *p, double *S)
{
    // todo work(root->vp, p);
    Node *currentNode = root;
    queryPoint *currentP;
    double tau = INFINITY;

    if (currentNode->inside == NULL && currentNode->outside == NULL)
    {
        return currentP;
    }
    else
    {

        // printf("\ntempDistances=\n");
        // for (int i = 0; i < p.k; ++i)
        //     printf("%lf ", tempDistances[i]);
        printf("\ntempDist ");
        printf("\ncurrentNode->p = %d \n", currentNode->p);
        //printf("\np.coord[i] = %lf \n", (*p).coord[0]);

        double tempDist = 0;

        // for (int i = 0; i < p->k; ++i)
        // {

        //     tempDist += pow(S[currentNode->p] -  p->coord[i], 2);
        // }
        // tempDist = sqrt(tempDist);

        printf("\ntempDist = %lf ", tempDist);
        printf("\ncurrentNode->mu = %lf \n", currentNode->mu);

        // if (tempDist < currentNode->mu)
        // {
        //     currentP = searchVTP(currentNode->inside, p, S);
        // }



        /*
        *   TODO: 
        *   debug query point and calculate distance from vp
        *   searchVPT right and left
        *   
        */




    }   
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
    Node *root = makeVPT(X, n, d);
    printf("VPT has been created.\n");

    queryPoint *p = NULL;

    //Assign memory to the pointer
    p = malloc(sizeof(queryPoint));
    if (p == NULL)
    {
        printf("fail to allocate memory for query point\n");
        return 0;
    }

    // p->coord[0] = 4.0;
    // p->coord[1] = 2.0;
    // p->k = 2;

    //printf("\np.coord[1] = %lf \n", p->coord[1]);
    searchVPT(root, p, X);

    free(X);
}