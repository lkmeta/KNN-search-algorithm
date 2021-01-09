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

knnresult distrAllkNN(double * X, int n, int d, int k){

    knnresult result;
    int count;
    int originRank;
    int originSize;
    Node *root = NULL;  //na dw pws tha kanw free to dentro
    int B=4;            //na to dw auto
    int* indexes = NULL;

    //auta tha theloun allagh logika meta
    MPI_Status stats[2*numtasks];   //Isend first, then Irecv
    MPI_Request reqs[2*numtasks];   //Isend first, then Irecv

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

    int prev = rank-1;
    int next = rank+1;
    if (rank == 0)  prev = numtasks - 1;
    if (rank == (numtasks - 1))  next = 0;

    //pinakas pou apothikeuei apo poio stoixeio tou X 3ekinaei to local set kathe process
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

    double* Y = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Y==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Y in process %d", rank);
        exit(-1);    
    }

    double* Ydists = (double*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Ydists==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Ydists in process %d", rank);
        exit(-1);    
    }

    int* Yidx = (int*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Yidx==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Yidx in process %d", rank);
        exit(-1);    
    }

    double* Z = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Z==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Z in process %d", rank);
        exit(-1);    
    }

    double* Zdists = (double*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Zdists==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Zdists in process %d", rank);
        exit(-1);    
    }

    int* Zidx = (int*)malloc((n/numtasks + 1)*k*sizeof(double));
    if(Zidx==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Zidx in process %d", rank);
        exit(-1);    
    }

    int* sendCounts = (int*)malloc(numtasks*sizeof(int));
    for(int i=0;i<numtasks;++i){
        if(i<numtasks-(n%numtasks)){
            sendCounts[i] = (n/numtasks)*d;
        }
        else{
            sendCounts[i] = (n/numtasks + 1)*d;
        } 
    }

    MPI_Scatterv(X,sendCounts,offsets,MPI_DOUBLE,Z,sendCounts[rank],MPI_DOUBLE,0,MPI_COMM_WORLD);

    if(rank==0){

        // printf("In process %d, size = %d\n", rank, (n/numtasks));
        // printf("In process %d X=\n", rank);
        // printMatrix(X,n*d);

        // printf("In process %d Z=\n", rank);
        // printMatrix(Z,n/numtasks*d);

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
            originSize = sendCounts[originRank];

            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            if(i==0){
                root = makeVPT(X, n/numtasks, d, indexes, B);
                result = kNN(X,X,root,n/numtasks,n/numtasks,d,k, offsets[rank]);
                memcpy(Ydists,result.ndist,n/numtasks*k*sizeof(double));
                memcpy(Yidx,result.nidx,n/numtasks*k*sizeof(int));
/*
                printf("\nIn process %d 1st iter ndist = \n", rank);
                printMatrix(result.ndist, sendCounts[rank]*k/d,k);

                printf("In process %d 1st iter nidx = \n",rank);
                for(int j=0;j<sendCounts[rank]*k/d;++j){
                    if(j%k==0 && j!=0){
                        printf("\n");
                    }   
                    printf("%d ", result.nidx[j]);
                }
                printf("\n");
*/
            }

            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));

            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            //printf("In process %d : now came %d elements from process %d\n",rank, count, stats[numtasks+rank].MPI_SOURCE);

            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks,d,k,offsets[rank]);
            }
            else{
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks+1,d,k,offsets[rank]);
            }

            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
            
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));

            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        //for each point in the local set of points, sort its neighbors based on their distance from it
        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }
//bgale ta sxolia apo katw gia testing se auto to branch

        // printf("\nndist = \n");
        // printMatrix(result.ndist, n/numtasks*k);
        // printf("nidx = \n");
        // for(int i=0;i<n/numtasks*k;++i){
        //     printf("%d ", result.nidx[i]);
        // }
        // printf("\n");
        free(newResult.ndist);
        free(newResult.nidx);
    }

    else if(rank<numtasks-(n%numtasks)){

        X = (double*)malloc((n/numtasks)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
            exit(-1);
        }

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

        // printf("In process %d, size = %d\n", rank, (n/numtasks));
        // printf("In process %d X=\n", rank);
        // printMatrix(X,(n/numtasks)*d);

        // printf("In process %d Z=\n", rank);
        // printMatrix(Z,n/numtasks*d);

        knnresult newResult;

        originRank = rank;

        for(int i=0;i<numtasks-1;++i){
            originSize = sendCounts[originRank];

            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            if(i==0){
                root = makeVPT(X, n/numtasks, d, indexes, B);
                result = kNN(X,X,root,n/numtasks,n/numtasks,d,k,offsets[rank]);
                memcpy(Ydists,result.ndist,n/numtasks*k*sizeof(double));
                memcpy(Yidx,result.nidx,n/numtasks*k*sizeof(int));
/*
                printf("\nIn process %d 1st iter ndist = \n", rank);
                printMatrix(result.ndist, sendCounts[rank]*k/d,k);

                printf("In process %d 1st iter nidx = \n",rank);
                for(int j=0;j<sendCounts[rank]*k/d;++j){
                    if(j%k==0 && j!=0){
                        printf("\n");
                    }   
                    printf("%d ", result.nidx[j]);
                }
                printf("\n");
*/
            }

            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));
/*
            printf("In process %d : now came %d elements from process %d\n",rank, count, stats[3].MPI_SOURCE);
            printMatrix(Y,count,d);
*/
            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks,d,k,offsets[rank]);
/*
                printf("\nIn process %d 2nd iter ndist = \n", rank);
                printMatrix(newResult.ndist, count*k/d,k);

                printf("In process %d 2nd iter nidx = \n",rank);
                for(int j=0;j<count*k/d;++j){
                    if(j%k==0 && j!=0){
                        printf("\n");
                    }   
                    printf("%d ", newResult.nidx[j]);
                }
                printf("\n");
*/
            }
            else{
                newResult = kNN(X,Y,root,n/numtasks,n/numtasks+1,d,k,offsets[rank]);
/*
                printf("\nIn process %d 2nd iter ndist = \n", rank);
                printMatrix(newResult.ndist, count*k/d,k);

                printf("In process %d 2nd iter nidx = \n",rank);
                for(int j=0;j<count*k/d;++j){
                    if(j%k==0 && j!=0){
                        printf("\n");
                    }   
                    printf("%d ", newResult.nidx[j]);
                }
                printf("\n");
*/
            }

            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
/*
            printf("\nYdists=\n");
            printMatrix(Ydists,count,k);
*/            
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));
/*
            printf("\nYidx=\n");
            for(int j=0;j<count;++j){
                if(j%k==0 && j!=0){
                    printf("\n");
                }
            printf("%d ", Yidx[j]);
            }
            printf("\n");
*/
            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }

//bgale ta sxolia apo katw gia testing se auto to branch

        // printf("\nIn process %d ndist = \n", rank);
        // printMatrix(result.ndist, n/numtasks*k);
        // printf("In process %d nidx = \n",rank);
        // for(int i=0;i<n/numtasks*k;++i){
        //     printf("%d ", result.nidx[i]);
        // }
        // printf("\n");

        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    else{
        
        X = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
        }

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

        //printf("In process %d, size = %d\n", rank, (n/numtasks + 1));
        //printf("In process %d X=\n", rank);
        //printMatrix(X,(n/numtasks + 1)*d);

        //printf("In process %d Z=\n", rank);
        //printMatrix(Z,(n/numtasks + 1)*d);

        knnresult newResult;

        originRank = rank;

        for(int i=0;i<numtasks-1;++i){
            originSize = sendCounts[originRank];

            MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, 0, MPI_COMM_WORLD, &reqs[3]);
            MPI_Irecv(Zdists,(n/numtasks + 1)*k,MPI_DOUBLE,prev,10,MPI_COMM_WORLD, &reqs[4]);
            MPI_Irecv(Zidx,(n/numtasks + 1)*k,MPI_INT,prev,20,MPI_COMM_WORLD,&reqs[5]);
            MPI_Isend(Y, originSize, MPI_DOUBLE, next, 0, MPI_COMM_WORLD, &reqs[0]);

            if(i==0){
                root = makeVPT(X, n/numtasks+1, d, indexes, B);
                result = kNN(X,X,root,n/numtasks+1,n/numtasks+1,d,k,offsets[rank]);
                memcpy(Ydists,result.ndist,(n/numtasks+1)*k*sizeof(double));
                memcpy(Yidx,result.nidx,(n/numtasks+1)*k*sizeof(int));
            }

            MPI_Isend(Ydists,originSize/d*k,MPI_DOUBLE,next,10,MPI_COMM_WORLD, &reqs[1]);
            MPI_Isend(Yidx,originSize/d*k,MPI_INT,next,20,MPI_COMM_WORLD,&reqs[2]);

            MPI_Wait(&reqs[3], &stats[3]);
            MPI_Get_count(&stats[3], MPI_DOUBLE, &count);  

            MPI_Wait(&reqs[0],&stats[0]);         
            memcpy(Y,Z,count*sizeof(double));

            // printf("In process %d : now came %d elements from process %d\n",rank, count, stats[numtasks+rank].MPI_SOURCE);

            originRank = (numtasks + stats[3].MPI_SOURCE - i) % numtasks;

            if(count==(n/numtasks)*d){
                newResult = kNN(X,Y,root,n/numtasks+1,n/numtasks,d,k,offsets[rank]);
            }
            else{
                newResult = kNN(X,Y,root,n/numtasks+1,n/numtasks+1,d,k,offsets[rank]);
            }

            MPI_Wait(&reqs[4],&stats[4]);
            MPI_Get_count(&stats[4],MPI_DOUBLE,&count);
            MPI_Wait(&reqs[1],&stats[1]);
            memcpy(Ydists,Zdists,count*sizeof(double));
            
            MPI_Wait(&reqs[5],&stats[5]);
            MPI_Get_count(&stats[5],MPI_INT,&count);
            MPI_Wait(&reqs[2],&stats[2]);
            memcpy(Yidx,Zidx,count*sizeof(int));

            mergeLists(Ydists,Yidx,newResult.ndist,newResult.nidx,count/k,k,0);
        }

        for(int i=0;i<count/k;++i){
            insertionSort(Ydists + i*k, Yidx + i*k, k);
        }

//bgale ta sxolia apo katw gia testing se auto to branch
/*
        printf("\nIn process %d ndist = \n", rank);
        printMatrix(result.ndist, (n/numtasks+1)*k);
        printf("In process %d nidx = \n",rank);
        for(int i=0;i<(n/numtasks+1)*k;++i){
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
*/
        free(newResult.ndist);
        free(newResult.nidx);
        free(X);
    }

    result.nidx = Yidx;
    result.ndist = Ydists;
    result.k = k;
    result.m = sendCounts[originRank]/d;
/*
    if(rank==0){
        printf("\nIn process %d ndist = \n", rank);
        printMatrix(result.ndist, sendCounts[rank]*k/d,k);

        printf("In process %d nidx = \n",rank);
        for(int i=0;i<sendCounts[rank]*k/d;++i){
            if(i%k==0 && i!=0){
            printf("\n");
            }
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
    }
*/

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
/*
        //for testing, remove later
        printf("\nIn process %d ndist = \n", rank);
        printMatrix(result.ndist, n*k,k);

        printf("In process %d nidx = \n", rank);
        for (int i = 0; i < n*k; ++i)
        {
            if(i%k==0 && i!=0){
            printf("\n");
            }   
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
*/    
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

int main(int argc, char* argv[]){

    srand(time(NULL));

    int n , d, k = 4;

    double *X = NULL;

    knnresult processResult;

    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0){

        X = readASCGZ("CoocTexture.asc",&n,&d);

/*
        X = (double *)malloc(n*d*sizeof(double));
        if(X==NULL){
            printf("Error in main: Couldn't allocate memory for X in process %d", rank);
            exit(-1);        
        }

        createRandomMatrix(X,n*d);
*/
/*
        X[0]=1.0;
        X[1]=3.0;
        X[2]=4.0;
        X[3]=2.0;
        X[4]=-2.0;
        X[5]=0.0;
        X[6]=5.0;
        X[7]=-1.0;
        X[8]=7.0;
        X[9]=-3.0;
        X[10]=-4.0;
        X[11]=6.0;
        X[12]=-8.0;
        X[13]=1.0;
        X[14]=3.0;
        X[15]=4.0;
        X[16]=-7.0;
        X[17]=5.0;
        X[18]=-1.0;
        X[19]=-2.0;
*/
        printf("Numtasks = %d\n", numtasks);
    }

    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init);

    processResult = distrAllkNN(X,n,d,k);

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