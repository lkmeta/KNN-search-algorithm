#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/**
* Function to read the .asc matrices given. At first, we define the number the number of points and dimensions
* of the matrix, allocate the needed amount of memory and read the matrix element by element.
* Input:
*       char* filename: string containing the name of the matrix
*       int* n: pointer to the number of points
*       int* d: pointer to the number of dimensions
* Output:
*       double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
**/
double *readCOL(char *filename, int *n, int *d)
{

    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    double doubleValue;

    int nCount = 1;
    int dCount = 0;

    int tempCounter = 0;
    int currPos = -1;

    //define number of dimensions
    while (fscanf(stream, "%lf ", &doubleValue))
    {
        //skip the first element of each row which is just the number of the row
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) != 0)
        {
            if (tempCounter == 0)
            {
                tempCounter++;
                continue;
            }
            else
                break;
        }
        else
        {
            dCount++;
        }
    }

    //define number of points
    rewind(stream);
    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {
        if (fabs(doubleValue - nCount) < 1e-6)
        {
            nCount++;
        }
    }

    nCount--;

    double *X = (double *)malloc(nCount * dCount * sizeof(double));

    rewind(stream);

    //store the values in X matrix
    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {
        //skip the first element of each row
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) != 0)
        {
            continue;
        }
        else
        {
            currPos++;
            X[currPos] = 10 * doubleValue;  //x10 so as not to have so small values and have problems of precision due to working with small doubles
        }
    }

    fclose(stream);

    *n = nCount;
    *d = dCount;

    return X;
}


/**
* Function to read the features matrix given. We skip the first 4 lines because they do not contain
* any points and then proceed to read the actual points of the matrix.
* Input:
*       char* filename: string containing the name of the matrix
*       int* n: pointer to the number of points
*       int* d: pointer to the number of dimensions
* Output:
*       double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
**/
double *readFEAT(char *filename, int *n, int *d)
{
    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    double doubleValue;

    int nCount = 106574;
    int dCount = 518;
    double *X = (double *)malloc(nCount * dCount * sizeof(double));

    int currPos = -1;

    //skip the first 4 lines
    char *line = (char *)malloc(1024 * 1024);
    for (int skip = 0; skip < 4; skip++)
    {
        int got = fscanf(stream, "%s\n", line);
    }
    free(line);

    //store the values in X matrix
    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) != 0)
        {
            continue;
        }
        else
        {
            currPos++;
            X[currPos] = doubleValue;
        }
    }

    fclose(stream);

    *n = nCount;
    *d = dCount;

    return X;
}


/**
* Function to read the MiniBooNE_PID matrix given. We skip the first 2 elements which do not represent
* actual points and then proceed to read the actual points of the matrix.
* Input:
*       char* filename: string containing the name of the matrix
*       int* n: pointer to the number of points
*       int* d: pointer to the number of dimensions
* Output:
*       double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
**/
double *readMINI(char *filename, int *n, int *d)
{
    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    double doubleValue;

    int nCount = 130064;
    int dCount = 50;
    double *X = (double *)malloc(nCount * dCount * sizeof(double));

    int currPos = -1;

    //skip the first 2 elements read
    int temp;
    for (int i = 0; i < 2; i++)
    {
        int got = fscanf(stream, "%d", &temp);
    }

    //store the values in X matrix
    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {
        currPos++;
        X[currPos] = doubleValue;
    }

    fclose(stream);

    printf("currPos=%d\n", currPos);

    *n = nCount;
    *d = dCount;

    return X;
}


/**
* Function to read the matrices contained in the tv file given. The number of points and dimensions is set 
* specifically for each of the matrices.
* Input:
*       char* filename: string containing the name of the matrix
*       int* n: pointer to the number of points
*       int* d: pointer to the number of dimensions
* Output:
*       double* X: nxd matrix containing the corpus data points (n points with d coordinates each)
**/
double *readTV(char *filename, int *n, int *d)
{
    int counter = -1;

    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    int nCount;
    
    //set the number of dimensions
    int dCount = 17;

    //set the number of points for each matrix
    if (strstr(filename, "BBC") != NULL)
    {
        nCount = 17720;
    }
    else if (strstr(filename, "CNN.") != NULL)
    {
        nCount = 22545;
    }
    else if (strstr(filename, "CNNI") != NULL)
    {
        nCount = 33117;
    }
    else if (strstr(filename, "NDTV") != NULL)
    {
        nCount = 17051;
    }
    else if (strstr(filename, "TIMES") != NULL)
    {
        nCount = 39252;
    }

    double *X = (double *)malloc(nCount * dCount * sizeof(double));

    double doubleValue;
    int intValue;

    //skip the first number
    fseek(stream, 1, SEEK_SET);

    //store the actual values in X
    for (int i = 0; i < nCount; i++)
    {
        //read only the 17 first in each line
        for (int j = 0; j < dCount; j++)
        {
            if (fscanf(stream, "%d:%lf ", &intValue, &doubleValue) != 2)
            {
                continue;
            }
            counter++;
            X[counter] = doubleValue;
        }
        //skip the rest in this line
        while (1)
        {
            if (fscanf(stream, "%d:%lf", &intValue, &doubleValue) != 2)
            {
                break;
            }
        }
    }

    fclose(stream);

    *n = nCount;
    *d = dCount;

    return X;
}
