#include <stdio.h>
#include "mpi.h"
#include "v0alt.c"

//Find the k smallest elements of the two lists
void mergeLists(knnresult old, knnresult new, int m, int k, int offset){

    double* ndistComb = (double*)malloc(2*k*sizeof(double));
    if(ndistComb==NULL){
        printf("Error in mergeLists: Couldn't allocate memory for ndistComb");
        exit(-1);    
    }

    int* nidxComb = (int*)malloc(2*k*sizeof(int));
    if(nidxComb==NULL){
        printf("Error in mergeLists: Couldn't allocate memory for nidxComb");
        exit(-1);    
    }

    for(int i=0;i<m;++i){
        for(int j=0;j<k;++j){
            ndistComb[j] = old.ndist[i*k+j];
            ndistComb[k+j] = new.ndist[i*k+j];
            nidxComb[j] = old.nidx[i*k+j];
            nidxComb[k+j] = new.nidx[i*k+j]+offset; //giati ta metraei apo to 0 kai den lamvanei ypopsin oti den einai sthn arxh tou d
        }
        printf("ndistComb=\n");
        printMatrix(ndistComb,2*k);

        printf("nidxComb = \n");
        for(int j=0;j<2*k;++j){
            printf("%d ", nidxComb[j]);
        }
        printf("\n");

        double kElem = quickSelect(ndistComb,nidxComb,0,2*k-1,k-1); //to kElem tha fugei meta

        printf("kElem=%lf\n", kElem);

        for(int j=0;j<k;++j){
            old.ndist[i*k+j] = ndistComb[j];
            old.nidx[i*k+j] = nidxComb[j];
        }
    }

    printf("Finally:\n");
    printf("ndist = \n");
    printMatrix(old.ndist, m*k);

    printf("nidx = \n");
    for(int i=0;i<m*k;++i){
        printf("%d ", old.nidx[i]);
    }
    printf("\n");


    free(ndistComb);
    free(nidxComb);
}


// Compute distributed all-kNN of points in X
/*

  \param  X      Data points                     [n-by-d]
  \param  n      Number of data points           [scalar]
  \param  d      Number of dimensions            [scalar]
  \param  k      Number of neighbors             [scalar]

  \return  The kNN result
*/
void distrAllkNN(double * X, int n, int d, int k){

    int numtasks, rank, count;
    int initialized, finalized;

    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    //auta tha theloun allagh logika meta
    MPI_Status stats[2*numtasks];   //Isend first, then Irecv
    MPI_Request reqs[2*numtasks];   //Isend first, then Irecv

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
            offsets[i] = offsets[i-1] + n/numtasks;
        }
        else{
            offsets[i] = offsets[i-1] + n/numtasks+1;
        }
    }

    double* Y = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Y==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Y in process %d", rank);
        exit(-1);    
    }

    double* Z = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
    if(Y==NULL){
        printf("Error in distrAllkNN: Couldn't allocate memory for Z in process %d", rank);
        exit(-1);    
    }

    if(rank==0){
        //srand(time(NULL));
        X = (double *)malloc(n*d*sizeof(double));
        if(X==NULL){
            printf("Error in distrAllkNN: Couldn't allocate memory for X");
            exit(-1);        
        }
        //createRandomMatrix(X,n*d);
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

        printf("Numtasks = %d\n", numtasks);

        printf("In process %d X=\n", rank);
        printMatrix(X,n*d);

        printf("In process %d, size = %d\n", rank, (n/numtasks));

        int dest = 1;

        //send local sets
        for(int i=1;i<numtasks;++i){
            if(dest<numtasks-(n%numtasks)){
                MPI_Isend(&X[offsets[i]*d], (n/numtasks)*d, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD, &reqs[i]);
            }
            else{
                MPI_Isend(&X[offsets[i]*d], (n/numtasks + 1)*d, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD, &reqs[i]);
            }
            dest++;
        }

//        knnresult result = kNN(X,X,n/numtasks,n/numtasks,d,k);
        knnresult newResult;

        MPI_Waitall(numtasks-1,&reqs[1],&stats[1]);

        MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[numtasks+rank]);
        MPI_Isend(X, (n/numtasks)*d, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &reqs[rank]);

        MPI_Wait(&reqs[numtasks+rank], &stats[numtasks+rank]);

        MPI_Get_count(&stats[numtasks+rank], MPI_DOUBLE, &count);
        Y = Z;

        printf("In process %d : now came %d elements from process %d\n",rank, count, stats[numtasks+rank].MPI_SOURCE);
/*
        if(count==(n/numtasks)*d){
            newResult = kNN(Y,X,n/numtasks,n/numtasks,d,k);
        }
        else{
            newResult = kNN(Y,X,n/numtasks+1,n/numtasks,d,k);
        }

        mergeLists(result,newResult,n/numtasks,k,offsets[stats[numtasks+rank].MPI_SOURCE]);
*/
/*
        printf("\nndist = \n");
        printMatrix(result.ndist, n/numtasks*k);

        printf("nidx = \n");
        for(int i=0;i<n/numtasks*k;++i){
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
*/
        free(X);
    }

    else if(rank<numtasks-(n%numtasks)){
        X = (double*)malloc((n/numtasks)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
            exit(-1);
        }
    
        printf("In process %d, size = %d\n", rank, (n/numtasks));

        int source = 0;

        MPI_Recv(X, (n/numtasks)*d, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &stats[numtasks+rank]);

        MPI_Get_count(&stats[numtasks+rank], MPI_DOUBLE, &count);
        printf("Process %d: Received %d elems from process %d with tag %d \n", rank, count, stats[numtasks+rank].MPI_SOURCE, stats[numtasks+rank].MPI_TAG);

        printf("In process %d X=\n", rank);
        printMatrix(X,(n/numtasks)*d);

        MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[numtasks+rank]);
        MPI_Isend(X, (n/numtasks)*d, MPI_DOUBLE, next, 10, MPI_COMM_WORLD, &reqs[rank]);
/*
        knnresult result = kNN(X,X,n/numtasks,n/numtasks,d,k);

        for(int i=0;i<n/numtasks*k;++i){
            result.nidx[i] += offsets[rank];
        }
*/
        knnresult newResult;

        MPI_Wait(&reqs[numtasks+rank], &stats[numtasks+rank]);

        MPI_Get_count(&stats[numtasks+rank], MPI_DOUBLE, &count);
        Y = Z;

        printf("In process %d : now came %d elements from process %d\n",rank, count, stats[numtasks+rank].MPI_SOURCE);
/*
        if(count==(n/numtasks)*d){
            newResult = kNN(Y,X,n/numtasks,n/numtasks,d,k);
        }
        else{
            newResult = kNN(Y,X,n/numtasks+1,n/numtasks,d,k);
        }

        mergeLists(result,newResult,n/numtasks,k,offsets[stats[numtasks+rank].MPI_SOURCE]);
*/
/*
        printf("\nIn process %d ndist = \n", rank);
        printMatrix(result.ndist, n/numtasks*k);

        printf("In process %d nidx = \n",rank);
        for(int i=0;i<n/numtasks*k;++i){
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
*/

        free(X);
    }

    else{
        X = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
        }

        printf("In process %d, size = %d\n", rank, (n/numtasks + 1));

        int source = 0;

        MPI_Recv(X, (n/numtasks + 1)*d, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &stats[numtasks+rank]);

        MPI_Get_count(&stats[numtasks+rank], MPI_DOUBLE, &count);
        printf("Process %d: Received %d elems from process %d with tag %d \n", rank, count, stats[numtasks+rank].MPI_SOURCE, stats[numtasks+rank].MPI_TAG);

        printf("In process %d X=\n", rank);
        printMatrix(X,(n/numtasks + 1)*d);

        MPI_Irecv(Z, (n/numtasks + 1)*d, MPI_DOUBLE, prev, MPI_ANY_TAG, MPI_COMM_WORLD, &reqs[numtasks+rank]);
        MPI_Isend(X, (n/numtasks+1)*d, MPI_DOUBLE, next, 11, MPI_COMM_WORLD, &reqs[rank]);

        knnresult result = kNN(X,X,n/numtasks+1,n/numtasks+1,d,k);
        for(int i=0;i<(n/numtasks+1)*k;++i){
            result.nidx[i] += offsets[rank];
        }

        knnresult newResult;

        MPI_Wait(&reqs[numtasks+rank], &stats[numtasks+rank]);

        MPI_Get_count(&stats[numtasks+rank], MPI_DOUBLE, &count);
        Y = Z;

        printf("In process %d : now came %d elements from process %d\n",rank, count, stats[numtasks+rank].MPI_SOURCE);

        if(count==(n/numtasks)*d){
            newResult = kNN(Y,X,n/numtasks,n/numtasks+1,d,k);
        }
        else{
            newResult = kNN(Y,X,n/numtasks+1,n/numtasks+1,d,k);
        }

        mergeLists(result,newResult,n/numtasks+1,k,offsets[stats[numtasks+rank].MPI_SOURCE]);

/*
        printf("\nIn process %d ndist = \n", rank);
        printMatrix(result.ndist, (n/numtasks+1)*k);

        printf("In process %d nidx = \n",rank);
        for(int i=0;i<(n/numtasks+1)*k;++i){
            printf("%d ", result.nidx[i]);
        }
        printf("\n");
*/

        free(X);
    }

    MPI_Finalized(&finalized);
    if (!finalized){
        MPI_Finalize();
    }
}

int main(int argc, char* argv[]){

    int n = 10, d = 2, k = 2;

    double *X = NULL;

    distrAllkNN(X,n,d,k);

}