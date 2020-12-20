#include <stdio.h>
#include "mpi.h"
#include "v0alt.c"

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

    MPI_Status stat;

    MPI_Initialized(&initialized);
    if (!initialized){
        MPI_Init(NULL, NULL);
    }

    MPI_Comm_size(MPI_COMM_WORLD, &numtasks);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if(rank==0){
        //srand(time(NULL));
        X = (double *)malloc(n * d * sizeof(double));
        if(X==NULL){
            printf("Error in main: Couldn't allocate memory for X");
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

        int offset=n/numtasks * d;
        printf("offset=%d\n", offset);
        int dest = 1;

        for(int i=1;i<numtasks;++i){
            if(dest<numtasks-(n%numtasks)){
                MPI_Send(&X[offset], (n/numtasks)*d, MPI_DOUBLE, dest, 1, MPI_COMM_WORLD);
                offset += (n/numtasks)*d;
                printf("offset=%d\n", offset);
            }
            else{
                MPI_Send(&X[offset], (n/numtasks + 1)*d, MPI_DOUBLE, dest, 2, MPI_COMM_WORLD);
                offset += (n/numtasks + 1)*d;
                printf("offset=%d\n", offset);
            }
            dest++;
        }

        knnresult result = kNN(X,X,n,n,d,k);

        printf("\nndist = \n");
        printMatrix(result.ndist, n*k);

        printf("nidx = \n");
        for(int i=0;i<n*k;++i){
            printf("%d ", result.nidx[i]);
        }
        printf("\n");

    }
    else if(rank<numtasks-(n%numtasks)){
        X = (double*)malloc((n/numtasks)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
            exit(-1);
        }
    
        printf("In process %d, size = %d\n", rank, (n/numtasks));

        int source = 0;

        MPI_Recv(X, (n/numtasks)*d, MPI_DOUBLE, source, 1, MPI_COMM_WORLD, &stat);

        MPI_Get_count(&stat, MPI_DOUBLE, &count);
        printf("Process %d: Received %d elems from process %d with tag %d \n", rank, count, stat.MPI_SOURCE, stat.MPI_TAG);

        printf("In process %d X=\n", rank);
        printMatrix(X,(n/numtasks)*d);
    }
    else{
        X = (double*)malloc((n/numtasks + 1)*d*sizeof(double));
        if(X==NULL){
            printf("In process %d X is NULL\n", rank);
        }

        printf("In process %d, size = %d\n", rank, (n/numtasks + 1));

        int source = 0;

        MPI_Recv(X, (n/numtasks + 1)*d, MPI_DOUBLE, source, 2, MPI_COMM_WORLD, &stat);

        MPI_Get_count(&stat, MPI_DOUBLE, &count);
        printf("Process %d: Received %d elems from process %d with tag %d \n", rank, count, stat.MPI_SOURCE, stat.MPI_TAG);

        printf("In process %d X=\n", rank);
        printMatrix(X,(n/numtasks + 1)*d);
    }

    MPI_Finalized(&finalized);
    if (!finalized){
        MPI_Finalize();
    }
}

int main(int argc, char* argv[]){

    int n = 10, d = 2, k = 3;

    double *X = NULL;

    distrAllkNN(X,n,d,k);

}