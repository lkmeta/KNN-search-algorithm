#include <stdio.h>
#include <string.h>

#include "mpi.h"
#include "test/tester.c"
#include "reader.c"

// Compute distributed all-kNN of points in X
/*

  \param  X      Data points                     [n-by-d]
  \param  n      Number of data points           [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]

  \return  The kNN result
*/

//global for easier use both in main and distrAllkNN
int numtasks, rank;

/**
 * Function that handles how each process receives and sends the data and computes the knn in between
 * communications. At first we send a part of X to each process (in such a way that each process either n/numtasks
 * or n/numtasks+1 elements). We compute the knn for the local points and then the processes participate in a
 * communication ring where the send their local points to their next process and receive the local points from their
 * previous. This is done "numtasks-1" times so that each local points set has passed from every process. Eventually,
 * each process has checked all points from the corpus set and has found the knn of its local points set.
 * Input: 
 *      double* X: nxd matrix containing the corpus data points (here it is both the corpus and the query points set)
 *      int n: number of corpus points
 *      int d: number of dimensions
 *      int k: number of nearest neighbors we are looking for
 * Output:
 *      knnresult: structure containing the info about the knn of each point of the process' local set
**/
knnresult distrAllkNN(double *X, int n, int d, int k)
{

    knnresult result;
    int count;      //number of elements received
    int originRank; //the process' rank for which the particular set of points is considered the local points set
    int originSize; //the number of elements of this set of points (not all processes have the same number of local points)

    MPI_Status stats[2]; //for Isend first, then Irecv
    MPI_Request reqs[2]; //for Isend first, then Irecv

    int *finalIndexes;      //array containing the indexes of the k-NNs for every point of X
    double *finalDistances; //array containing the distances between every point of X and its k-NNs

    //the struct containing the kNNs of all points will only be at the root process
    if(rank==0){
        finalIndexes = (int *)malloc(n * k * sizeof(int));
        if (finalIndexes == NULL)
        {
            printf("Error in distrAllkNN: Couldn't allocate memory for finalIndexes in process %d", rank);
            exit(-1);
        }

        finalDistances = (double *)malloc(n * k * sizeof(double));
        if (finalDistances == NULL)
        {
            printf("Error in distrAllkNN: Couldn't allocate memory for finalDistances in process %d", rank);
            exit(-1);
        }
    }
    else
    {
        finalIndexes = NULL;
        finalDistances = NULL;
    }
    

    //define the 'neighbors' of each process in the communication ring
    int prev = rank - 1;
    int next = rank + 1;
    if (rank == 0)
        prev = numtasks - 1;
    if (rank == (numtasks - 1))
        next = 0;

    /*
    * Note: In order to achieve greater balnce when spliting the initial X matrix we applied the following:
    * a)numtasks-(n%numtasks) processes have n/numtasks points
    * b)n%numtasks processes have n/numtasks + 1 points
    * This is always better or equally good than the classic split of n/numtasks elements in each process and the last
    * one gets all the remaining points too.
    * i.e if we have 4 processes and 11 points our split gives 3 processes with 3 points and 1 process with 2 points
    * while the classic split would give 3 processes with 2 points and one process with 5 points which is claerly more
    * unbalanced.
    */

    //array containing from which point of X does the local set of each process begin
    int *offsets = (int *)malloc(numtasks * sizeof(int));
    if (offsets == NULL)
    {
        printf("Error in distrAllkNN: Couldn't allocate memory for offsets in process %d", rank);
        exit(-1);
    }

    offsets[0] = 0;
    for (int i = 1; i < numtasks; ++i)
    {
        if (i <= numtasks - (n % numtasks))
        {
            offsets[i] = offsets[i - 1] + n / numtasks * d;
        }
        else
        {
            offsets[i] = offsets[i - 1] + (n / numtasks + 1) * d;
        }
    }

    //array containing the new corpus set points
    //size set n/numtasks+1 points because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    double *Y = (double *)malloc((n / numtasks + 1) * d * sizeof(double));
    if (Y == NULL)
    {
        printf("Error in distrAllkNN: Couldn't allocate memory for Y in process %d", rank);
        exit(-1);
    }

    //array where we temporarily store the newly arrived data to each process
    //size set n/numtasks+1 points because we are not able to know whether n/numtasks or n/numtasks + 1 points
    //are about to arrive
    double *Z = (double *)malloc((n / numtasks + 1) * d * sizeof(double));
    if (Y == NULL)
    {
        printf("Error in distrAllkNN: Couldn't allocate memory for Z in process %d", rank);
        exit(-1);
    }

    //array containing how many points of X scatterv will send to each process
    //also works as an array containing the number of local points of each process
    int *sendCounts = (int *)malloc(numtasks * sizeof(int));
    for (int i = 0; i < numtasks; ++i)
    {
        if (i < numtasks - (n % numtasks))
        {
            sendCounts[i] = (n / numtasks) * d;
        }
        else
        {
            sendCounts[i] = (n / numtasks + 1) * d;
        }
    }

    //send the local points to each process, 0 works as the root process here
    MPI_Scatterv(X, sendCounts, offsets, MPI_DOUBLE, Z, sendCounts[rank], MPI_DOUBLE, 0, MPI_COMM_WORLD);

    //if we are in the root process
    if (rank == 0)
    {

        memcpy(Y, Z, n / numtasks * d * sizeof(double));

        knnresult newResult;

        originRank = rank;

        for (int i = 0; i < numtasks - 1; ++i)
        {
            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n / numtasks + 1) * d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration compute the knn of the local points from the local points
            if (i == 0)
            {
                result = kNN(X, X, n / numtasks, n / numtasks, d, k);
            }

            //wait for the new points to arrive
            MPI_Wait(&reqs[1], &stats[1]);

            //check how many points arrived
            MPI_Get_count(&stats[1], MPI_DOUBLE, &count);

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0], &stats[0]);

            memcpy(Y, Z, count * sizeof(double));

            //apply kNN using the newly arrived points as the new corpus set and the local points as the query set
            if (count == (n / numtasks) * d)
            {
                newResult = kNN(Y, X, n / numtasks, n / numtasks, d, k);
            }
            else
            {
                newResult = kNN(Y, X, n / numtasks + 1, n / numtasks, d, k);
            }

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[1].MPI_SOURCE - i) % numtasks;

            mergeLists(result, newResult, n / numtasks, k, offsets[originRank] / d);

        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<n/numtasks;++i){
            insertionSort(result.ndist + i*k, result.nidx + i*k, k);
        }

        free(newResult.ndist);
        free(newResult.nidx);
    }

    //if the process contains n/numtasks local points
    else if (rank < numtasks - (n % numtasks))
    {

        X = (double *)malloc((n / numtasks) * d * sizeof(double));
        if (X == NULL)
        {
            printf("In process %d X is NULL\n", rank);
            exit(-1);
        }

        memcpy(X, Z, (n / numtasks) * d * sizeof(double));
        memcpy(Y, Z, (n / numtasks) * d * sizeof(double));

        knnresult newResult;

        originRank = rank;

        for (int i = 0; i < numtasks - 1; ++i)
        {
            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n / numtasks + 1) * d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration compute the knn of the local points from the local points
            //additionaly add the offset these elements have in the original X matrix
            if (i == 0)
            {
                result = kNN(X, X, n / numtasks, n / numtasks, d, k);
                for (int i = 0; i < n / numtasks * k; ++i)
                {
                    result.nidx[i] += offsets[rank] / d;
                }
            }

            //wait for the new points to arrive
            MPI_Wait(&reqs[1], &stats[1]);

            //check how many points arrived
            MPI_Get_count(&stats[1], MPI_DOUBLE, &count);

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0], &stats[0]);

            memcpy(Y, Z, count * sizeof(double));

            //apply kNN using the newly arrived points as the new corpus set and the local points as the query set
            if (count == (n / numtasks) * d)
            {
                newResult = kNN(Y, X, n / numtasks, n / numtasks, d, k);
            }
            else
            {
                newResult = kNN(Y, X, n / numtasks + 1, n / numtasks, d, k);
            }

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[1].MPI_SOURCE - i) % numtasks;

            mergeLists(result, newResult, n / numtasks, k, offsets[originRank] / d);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<n/numtasks;++i){
            insertionSort(result.ndist + i*k, result.nidx + i*k, k);
        }
        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    //if the process contains n/numtasks + 1 local points
    else
    {

        X = (double *)malloc((n / numtasks + 1) * d * sizeof(double));
        if (X == NULL)
        {
            printf("In process %d X is NULL\n", rank);
        }

        memcpy(X, Z, (n / numtasks + 1) * d * sizeof(double));
        memcpy(Y, Z, (n / numtasks + 1) * d * sizeof(double));

        knnresult newResult;

        originRank = rank;

        for (int i = 0; i < numtasks - 1; ++i)
        {
            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n / numtasks + 1) * d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 11, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration compute the knn of the local points from the local points
            //additionaly add the offset these elements have in the original X matrix
            if (i == 0)
            {
                result = kNN(X, X, n / numtasks + 1, n / numtasks + 1, d, k);
                for (int i = 0; i < (n / numtasks + 1) * k; ++i)
                {
                    result.nidx[i] += offsets[rank] / d;
                }
            }

            //wait for the new points to arrive
            MPI_Wait(&reqs[1], &stats[1]);

            //check how many points arrived
            MPI_Get_count(&stats[1], MPI_DOUBLE, &count);

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0], &stats[0]);

            memcpy(Y, Z, count * sizeof(double));

            //apply kNN using the newly arrived points as the new corpus set and the local points as the query set
            if (count == (n / numtasks) * d)
            {
                newResult = kNN(Y, X, n / numtasks, n / numtasks + 1, d, k);
            }
            else
            {
                newResult = kNN(Y, X, n / numtasks + 1, n / numtasks + 1, d, k);
            }

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[1].MPI_SOURCE - i) % numtasks;

            mergeLists(result, newResult, n / numtasks + 1, k, offsets[originRank] / d);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<n/numtasks+1;++i){
            insertionSort(result.ndist + i*k, result.nidx + i*k, k);
        }
        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    //offsets is now used to define in which position will the ndist and nidx elements of each process go when we'll
    //gather them in the root process
    for (int i = 1; i < numtasks; ++i)
    {
        if (i <= numtasks - (n % numtasks))
        {
            offsets[i] = offsets[i - 1] + n / numtasks * k;
        }
        else
        {
            offsets[i] = offsets[i - 1] + (n / numtasks + 1) * k;
        }
    }

    //similarly, sendCounts now defines how many elements of ndist and nidx will be sent from each process
    //to the root process where we'll gather them together
    for (int i = 0; i < numtasks; ++i)
    {
        if (i < numtasks - (n % numtasks))
        {
            sendCounts[i] = (n / numtasks) * k;
        }
        else
        {
            sendCounts[i] = (n / numtasks + 1) * k;
        }
    }

    MPI_Gatherv(result.nidx, sendCounts[rank], MPI_INT, finalIndexes, sendCounts, offsets, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Gatherv(result.ndist, sendCounts[rank], MPI_DOUBLE, finalDistances, sendCounts, offsets, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0)
    {
        //deallocate the arrays containing the kNNs of the local points
        free(result.nidx);
        free(result.ndist);

        //now the knn struct in the root process contains the kNNs of all elements of X
        result.nidx = finalIndexes;
        result.ndist = finalDistances;

    }

    free(Y);
    free(Z);
    free(offsets);
    free(sendCounts);

    return result;
}

//Main has either 2 or 3 command line arguments
//Case of 3: The points are read from a file
//Argument 1: the filename, Argument 2: value of k
//Case of 4: The points are created randomly
//Argument 1: value of n, Argument 2: value of d, Argument 3: value of k
int main(int argc, char *argv[])
{
    int n,d,k;
    
    double *X = NULL;

    knnresult processResult;

    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized)
    {
        MPI_Init(&argc, &argv);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //only the root process reads the whole X matrix and then scatters it to all processes
    if (rank == 0)
    {
        if(argc<3){
            printf("Not enough command line arguments\n");
            exit(-1);
        }

        else if(argc==3){
            char* s=argv[1];
            uint nameLength = strlen(s);    //length of the name of the file 
            if((s[nameLength-1]=='c') && (s[nameLength-2]=='s') && (s[nameLength-3]=='a') && (s[nameLength-4]=='.')){
                printf("Your argument is a .asc file\n");
                X = readASCGZ(s,&n,&d);
            }
            else{
                printf("Not a .asc file!\n");
                exit(-1);
            }
        }

        else if(argc==4){
            srand(time(NULL));

            n = atoi(argv[1]);
            d = atoi(argv[2]);
            k = atoi(argv[3]);

            X = (double*)malloc(n*d*sizeof(double));
            if (X == NULL)
            {
                printf("Error in main: Couldn't allocate memory for X");
                exit(-1);
            }

            createRandomMatrix(X,n*d);
        }

        else{
            printf("Too many command line arguments\n");
            exit(-1);
        }
    }

    if(argc==3){

        //send n,d to each process
        MPI_Bcast(&n,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&d,1,MPI_INT,0,MPI_COMM_WORLD);

        //define k
        k=atoi(argv[2]);
    }
    else if(argc==4){

        //define n,d,k
        n = atoi(argv[1]);
        d = atoi(argv[2]);
        k = atoi(argv[3]);
    }

    printf("argc=%d, n=%d, d=%d, k=%d\n",argc,n,d,k);
    
    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init);

    processResult = distrAllkNN(X, n, d, k);

    //End timer
    struct timespec last;   
    clock_gettime(CLOCK_MONOTONIC, &last);

    long ns;
    uint seconds;
    if(last.tv_nsec <init.tv_nsec){
        ns=init.tv_nsec - last.tv_nsec;
        seconds= last.tv_sec - init.tv_sec -1;
    }

    if(last.tv_nsec >init.tv_nsec){
        ns= last.tv_nsec -init.tv_nsec ;
        seconds= last.tv_sec - init.tv_sec ;
    }

    if(rank==0){

        printf("argc=%d, n = %d, d = %d, k = %d, numoftasks = %d\n", argc, n, d, k, numtasks);
        printf("For V1 the seconds elapsed are %u and the nanoseconds are %ld\n", seconds, ns);
        
        free(X);
    }

    free(processResult.nidx);
    free(processResult.ndist);

    MPI_Finalized(&finalized);
    if (!finalized)
    {
        MPI_Finalize();
    }

    return 0;
}