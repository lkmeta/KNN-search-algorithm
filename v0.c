#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

#include "knn.h"
#include "reader.c"


//Main has either 2 or 3 command line arguments
//Case of 2: The points are read from a file
//Argument 1: the filename, Argument 2: value of k
//Case of 3: The points are created randomly
//Argument 1: value of n, Argument 2: value of d, Argument 3: value of k
int main(int argc, char *argv[])
{
    int n, d, k;

    double *X = NULL;

    if (argc < 3)
    {
        printf("Not enough command line arguments\n");
        exit(-1);
    }

    else if (argc == 3)
    {

        char *s = argv[1];
        uint nameLength = strlen(s); //length of the name of the file

        if (strstr(s, "Color") != NULL || strstr(s, "Cooc") != NULL || strstr(s, "Layout") != NULL)
        {
            printf("Your argument matrix is %s file\n", s);
            X = readCOL(s, &n, &d);
        }

        else if (strstr(s, "features") != NULL)
        {
            printf("Your argument matrix is %s file\n", s);
            X = readFEAT(s, &n, &d);
        }

        else if (strstr(s, "Mini") != NULL)
        {
            printf("Your argument matrix is %s file\n", s);
            X = readMINI(s, &n, &d);
        }

        else if (strstr(s, "BBC") != NULL || strstr(s, "CNN.") != NULL || strstr(s, "CNNI") != NULL || strstr(s, "NDT") != NULL || strstr(s, "TIME") != NULL)
        {
            printf("Your argument matrix is %s file\n", s);
            X = readTV(s, &n, &d);
        }

        else
        {
            printf("Cannot read this file!\n");
            exit(-1);
        }

        k = atoi(argv[2]);
    }

    else if (argc == 4)
    {

        n = atoi(argv[1]);
        d = atoi(argv[2]);
        k = atoi(argv[3]);

        X = (double *)malloc(n * d * sizeof(double));
        if (X == NULL)
        {
            printf("Error in main: Couldn't allocate memory for X");
            exit(-1);
        }

        createRandomMatrix(X, n * d);
    }

    else
    {
        printf("Too many command line arguments\n");
        exit(-1);
    }

    printf("argc=%d, n = %d, d = %d, k = %d\n", argc, n, d, k);

    //Start timer
    struct timespec init;
    clock_gettime(CLOCK_MONOTONIC, &init);

    knnresult result = kNN(X, X, n, n, d, k);

    //End timer
    struct timespec last;
    clock_gettime(CLOCK_MONOTONIC, &last);

    long ns;
    uint seconds;
    if (last.tv_nsec < init.tv_nsec)
    {
        ns = init.tv_nsec - last.tv_nsec;
        seconds = last.tv_sec - init.tv_sec - 1;
    }

    if (last.tv_nsec > init.tv_nsec)
    {
        ns = last.tv_nsec - init.tv_nsec;
        seconds = last.tv_sec - init.tv_sec;
    }

    printf("For V0 the seconds elapsed are %u and the nanoseconds are %ld\n", seconds, ns);

    //deallocate used memory
    free(X);

    free(result.nidx);
    free(result.ndist);

    return 0;
}
