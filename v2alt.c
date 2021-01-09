#include <stdio.h>
#include <string.h>
#include "mpi.h"
#include "test/tester2.c"
#include "reader.c"

//global for easier use both in main and distrAllkNN
int numtasks, rank;


// Compute distributed all-kNN of points in X
/*
  \param  X      Data points                     [n-by-d]
  \param  n      Number of data points           [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]
  \return  The kNN result
*/


/**
 * Function that handles how each process receives and sends the data and computes the knn in between
 * communications. distrAllkNN works in a similar fashion to that in the v1.c file but there are some main differences.
 * First of all, each process creates a VPT from its local points and then computes for the nearest neighbors for
 * these local points. Then, each process, apart from sending the points it has to the next process also sends two 
 * arrays accompanying these points: an array containing the neighbors for each of these points and an array
 * containing the indexes of those neighbors. The kNN algorithm is executed again for these points but now for
 * a different VPT since they are in another process and their neighbors are updated. In the end all the neighbors
 * are gathered in the root process.
 * Input: 
 *      double* X: nxd matrix containing the corpus data points (here it is both the corpus and the query points set)
 *      int n: number of corpus points
 *      int d: number of dimensions
 *      int k: number of nearest neighbors we are looking for
 * Output:
 *      knnresult: structure containing the info about the knn of each point of the process' local set (except the
 *      structure in the root process which contains info for all points and not only the local)
**/
knnresult distrAllkNN(double * X, int n, int d, int k){

    knnresult result;   //struct to be returned
    int count;          //number of elements received
    int originRank;     //the process' rank for which the particular set of points is considered the local points set
    int originSize;     //the number of elements of this set of points (not all processes have the same number of local points)
    Node *root = NULL;  //the root of the VPT for each process
    int B=4;            //na to dw auto
    
    int* indexes = NULL;    //array containing the indexes the local points have in the array that contains all the points
 
    MPI_Status stats[6];   //First 3 for the Isends , last 3 for the Irecvs
    MPI_Request reqs[6];   //First 3 for the Isends , last 3 for the Irecvs

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
    int prev = rank-1;
    int next = rank+1;
    if (rank == 0)  prev = numtasks - 1;
    if (rank == (numtasks - 1))  next = 0;

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
    int* offsets = (int*)malloc(numtasks*sizeof(int));
    if(offsets==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for offsets in process %d", rank);
        exit(-1);    
    }

    offsets[0] = 0;
    for(int i=1;i<numtasks;++i){
        if(i<=numtasks-(n%numtasks)){
            offsets[i] = offsets[i-1] + n/numtasks*d;
        }
        else{
            offsets[i] = offsets[i-1] + (n/numtasks+1)*d;
        }
    }

    //array containing the new corpus set points
    //size set n/numtasks+1 points because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    double* Y = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Y==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Y in process %d", rank);
        exit(-1);    
    }

    //array containing the distances between each of the newly arrived points in Y and their neighbors
    //size set (n/numtasks+1)*k because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    double* Ydists = (double*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Ydists==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Ydists in process %d", rank);
        exit(-1);    
    }

    //array containing the indexes of the neighbors of the newly arrived points in Y
    //size set (n/numtasks+1)*k because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    int* Yidx = (int*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Yidx==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Yidx in process %d", rank);
        exit(-1);    
    }

    //array where we temporarily store the newly arrived data to each process
    //size set n/numtasks+1 points because we are not able to know whether n/numtasks or n/numtasks + 1 points
    //are about to arrive
    double* Z = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Z==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Z in process %d", rank);
        exit(-1);    
    }

    //array where we temporarily store the distances between each of the newly arrived points and their neighbors
    //size set (n/numtasks+1)*k because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    double* Zdists = (double*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Zdists==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Zdists in process %d", rank);
        exit(-1);    
    }

    //array where we temporarily store the indexes of the neighbors of the newly arrived points
    //size set (n/numtasks+1)*k because we are not able to know whether n/numtasks or n/numtasks +1 points
    //are about to arrive
    int* Zidx = (int*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Zidx==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Zidx in process %d", rank);
        exit(-1);    
    }

    //array containing how many points of X scatterv will send to each process
    //also works as an array containing the number of local points of each process
    int* sendCounts = (int*)malloc(numtasks*sizeof(int));
    for(int i=0;i<numtasks;++i){
        if(i<numtasks-(n%numtasks)){
            sendCounts[i] = (n/numtasks)*d;
        }
        else{
            sendCounts[i] = (n/numtasks + 1)*d;
        } 
    }

    //send the local points to each process, 0 works as the root process here
    MPI_Scatterv(X,sendCounts,offsets,MPI_DOUBLE,Z,sendCounts[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rank==0){

        //set the correct indexes for the points taking into account the offset these local points have from the
        //beginning of the matrix that contains all points
        indexes = (int *)malloc(n/numtasks * sizeof(int));
        if (indexes == NULL)
        {
            printf("Error in process %d: Couldn't allocate memory for indexes",rank);
            exit(-1);
        }

        for (int i = 0; i < n/numtasks; ++i)
        {
            indexes[i] = i + offsets[rank]/d;
        }

        memcpy(Y,Z,n/numtasks*d*sizeof(double));

        knnresult newResult;

        originRank = rank;

        for(int i=0;i<numtasks-1;++i){

            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration create the local VPT and compute the knn of the local points from the local points
            if(i==0){
                root = makeVPT(X, n/numtasks, d, indexes, B);
                result = kNN(X,X,root,n/numtasks,n/numtasks,d,k, offsets[rank]);
                memcpy(Ydists,result.ndist,n/numtasks*k*sizeof(double));
                memcpy(Yidx,result.nidx,n/numtasks*k*sizeof(int));
            }

            //send the arrays containing the neighbors and the neighbors' indexes for the sent points in Y
            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            //wait for the new points to arrive and check how many points arrived
            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            //apply kNN using the newly arrived points as the new corpus set and the local points as those in the VPT
            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks,d,k,offsets[rank]);
            }
            else{
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks+1,d,k,offsets[rank]);
            }

            //wait for the neighbors to arrive and get their count
            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);

            //wait for the neighbors of the previous points to be sent
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
            
            //wait for the neighbors' indexes to arrive and get their count
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);

            //wait for the neighbors' indexes of the previous points to be sent
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));

            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }

        free(newResult.ndist);
        free(newResult.nidx);
    }

    //if the process contains n/numtasks local points
    else if(rank<numtasks-(n%numtasks)){

        X = (double*)malloc((n/numtasks)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
            exit(-1);
        }

        //set the correct indexes for the points taking into account the offset these local points have from the
        //beginning of the matrix that contains all points
        indexes = (int *)malloc(n/numtasks * sizeof(int));
        if (indexes == NULL)
        {
            printf("Error in process %d: Couldn't allocate memory for indexes",rank);
            exit(-1);
        }

        for (int i = 0; i < n/numtasks; ++i)
        {
            indexes[i] = i + offsets[rank]/d;
        }
        printf("\n");

        memcpy(X,Z,(n/numtasks)*d*sizeof(double));
        memcpy(Y,Z,(n/numtasks)*d*sizeof(double));

        knnresult newResult;

        originRank = rank;

        for(int i=0;i<numtasks-1;++i){

            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration create the local VPT and compute the knn of the local points from the local points
            if(i==0){
                root = makeVPT(X, n/numtasks, d, indexes, B);
                result = kNN(X,X,root,n/numtasks,n/numtasks,d,k,offsets[rank]);
                memcpy(Ydists,result.ndist,n/numtasks*k*sizeof(double));
                memcpy(Yidx,result.nidx,n/numtasks*k*sizeof(int));
            }

            //send the arrays containing the neighbors and the neighbors' indexes for the sent points in Y
            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            //wait for the new points to arrive and check how many points arrived
            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            //apply kNN using the newly arrived points as the new corpus set and the local points as those in the VPT
            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks,d,k,offsets[rank]);

            }
            else{
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks+1,d,k,offsets[rank]);

            }

            //wait for the neighbors to arrive and get their count
            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);

            //wait for the neighbors of the previous points to be sent
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
            
            //wait for the neighbors' indexes to arrive and get their count
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);

            //wait for the neighbors' indexes of the previous points to be sent
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));

            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }

        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    //if the process contains n/numtasks + 1 local points
    else{
        
        X = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
        }

        //set the correct indexes for the points taking into account the offset these local points have from the
        //beginning of the matrix that contains all points
        indexes = (int *)malloc((n/numtasks +1) * sizeof(int));
        if (indexes == NULL)
        {
            printf("Error in process %d: Couldn't allocate memory for indexes",rank);
            exit(-1);
        }

        for (int i = 0; i < n/numtasks+1; ++i)
        {
            indexes[i] = i + offsets[rank]/d;
        }

        memcpy(X,Z,(n/numtasks+1)*d*sizeof(double));
        memcpy(Y,Z,(n/numtasks+1)*d*sizeof(double));

        knnresult newResult;

        originRank = rank;

        for(int i=0;i<numtasks-1;++i){

            //size of the data that previously arrived and now have to be sent to the next process
            originSize = sendCounts[originRank];

            //non blocking sends and receives so that we can compute knn for the local set while waiting for the new
            //points to arrive
            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            //in the first iteration create the local VPT and compute the knn of the local points from the local points
            if(i==0){
                root = makeVPT(X, n/numtasks+1, d, indexes, B);
                result = kNN(X,X,root,n/numtasks+1,n/numtasks+1,d,k,offsets[rank]);
                memcpy(Ydists,result.ndist,(n/numtasks+1)*k*sizeof(double));
                memcpy(Yidx,result.nidx,(n/numtasks+1)*k*sizeof(int));
            }

            //send the arrays containing the neighbors and the neighbors' indexes for the sent points in Y
            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            //wait for the new points to arrive and check how many points arrived
            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            //wait for the previous points to be sent
            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));

            //check from which process these points originally came from
            //(in which process they are considered local points)
            //we add numtasks and then take the remainder of division with numtasks so that no negative ranks occur
            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            //apply kNN using the newly arrived points as the new corpus set and the local points as those in the VPT
            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks+1,n/numtasks,d,k,offsets[rank]);
            }
            else{
                newResult = kNN(X,Y,root,n/numtasks+1,n/numtasks+1,d,k,offsets[rank]);
            }

            //wait for the neighbors to arrive and get their count
            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);

            //wait for the neighbors of the previous points to be sent
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
            
            //wait for the neighbors' indexes to arrive and get their count
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);

            //wait for the neighbors' indexes of the previous points to be sent
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));

            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }

        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    result.nidx = Yidx;
    result.ndist = Ydists;
    result.k = k;
    result.m = sendCounts[originRank]/d;

    //each process now contains the arrays of neighbors for the local points of their next process
    //this is why we are using this temp variable to swap the places of the offsets in the offsets array
    int temp;

    //offsets is now used to define in which position will the ndist and nidx elements of each process go when we'll
    //gather them in the root process
    for (int i = 0; i < numtasks; ++i)
    {
        offsets[i] = offsets[i]/d * k;
    }

    for(int i=0; i < numtasks; ++i){
        if(i==0) temp = offsets[i];
        if(i==(numtasks-1)){
            offsets[i] = temp; 
            break;   
        }
        offsets[i] = offsets[i+1];

    }

    //similarly, sendCounts now defines how many elements of ndist and nidx will be sent from each process
    //to the root process where we'll gather them together
    for(int i=0; i<numtasks; ++i){
        sendCounts[i] = sendCounts[i]/d * k;
    }

    for(int i=0; i<numtasks; ++i){
        if(i==0) temp = sendCounts[i];
        if(i==(numtasks-1)){
            sendCounts[i] = temp; 
            break;   
        }
        sendCounts[i] = sendCounts[i+1];
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

    free(sendCounts);
    free(indexes);
    free(Y);
    free(Zdists);
    free(Zidx);
    free(Z);
    free(offsets);
    freeNode(root);
    free(root);

    return result;
}

//Main has either 2 or 3 command line arguments
//Case of 3: The points are read from a file
//Argument 1: the filename, Argument 2: value of k
//Case of 4: The points are created randomly
//Argument 1: value of n, Argument 2: value of d, Argument 3: value of k
int main(int argc, char* argv[]){

    srand(time(NULL));

    int n,d,k;

    double *X = NULL;

    knnresult processResult;

    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(&argc, &argv);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //only the root process reads the whole X matrix and then scatters it to all processes
    if(rank==0){
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

    processResult = distrAllkNN(X,n,d,k);

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

        printf("For V2 the seconds elapsed are %u and the nanoseconds are %ld\n",seconds, ns);

        //comfirm the validity of our results using the tester provided
        checkResult(processResult,X,X,n,n,d,k);
        
        free(X);
    }

    free(processResult.nidx);
    free(processResult.ndist);

    MPI_Finalized(&finalized);
    if (!finalized){
        MPI_Finalize();
    }

}