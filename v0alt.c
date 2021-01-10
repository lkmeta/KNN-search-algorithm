#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>
    
#include "test/tester.c" 
#include "reader.c"

int main(int argc, char* argv[])
{   
    int n,d,k;

    double *X = NULL;

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

        k=atoi(argv[2]);
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

    printf("argc=%d, n=%d, d=%d, k=%d\n",argc,n,d,k);

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

    //comfirm the validity of our results using the tester provided
    //checkResult(result,X,X,n,n,d,k);

    //deallocate used memory 
    free(X);

    free(result.nidx);
    free(result.ndist);

    return 0;
}