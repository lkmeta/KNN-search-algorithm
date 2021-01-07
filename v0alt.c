#include <stdio.h>
#include <stdlib.h>
#include <time.h>
    
#include "test/tester.c" 

int main(int argc, char* argv[])
{   
    //currently change them by hand, will fix later
    int n = 10000, m = 4, d = 10, k = 10;

    srand(time(NULL));

    // Allocate memory for X and Y lists
    double *X = (double *)malloc(n * d * sizeof(double));
    double *Y = (double *)malloc(m * d * sizeof(double));

    if(X==NULL){
        printf("Error in main: Couldn't allocate memory for X");
        exit(-1);
    }

    if(Y==NULL){
        printf("Error in main: Couldn't allocate memory for Y");
        exit(-1);
    }

    //create set X and y
    createRandomMatrix(X,n*d);
//    createRandomMatrix(Y,m*d);
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
    X[20]=-6.0;
    X[21]=-4.0;

    Y[0] = 1.0;
    Y[1] = 5.0;
    Y[2] = -2.0;
    Y[3] = 4.0;
    Y[4] = -3.0;
    Y[5] = -4.0;
    Y[6] = 2.0;
    Y[7] = -6.0;
*/
//    printf("X= \n");
//    printMatrix(X, n*d);

//    printf("Y= \n");
//    printMatrix(Y, m*d);

    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init);

    knnresult result = kNN(X,X,n,n,d,k);

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
    printf("For V0 the seconds elapsed are %u and the nanoseconds are %ld\n",seconds, ns);
/*
    printf("\nndist = \n");
    printMatrix(result.ndist, n*k,k);

    printf("nidx = \n");
    for(int i=0;i<n*k;++i){
        if(i%k==0 && i!=0){
            printf("\n");
        }
        printf("%d ", result.nidx[i]);
    }
    printf("\n");
*/
    //comfirm the validity of our results using the tester provided
    checkResult(result,X,X,n,n,d,k);

    //deallocate used memory 
    free(X);
    free(Y);

    free(result.nidx);
    free(result.ndist);

    return 0;
}