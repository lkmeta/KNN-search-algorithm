#include <stdio.h>
#include <stdlib.h>
#include <math.h>


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

    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) != 0)
        {
            continue;
        }
        else
        {
            currPos++;
            X[currPos] = 10 * doubleValue;
        }
    }

    fclose(stream);

    *n = nCount;
    *d = dCount;

    return X;
}

double *readFEAT(char *filename, int *n, int *d)
{
    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    double doubleValue;

    int nCount = 1038;
    int dCount = 518;
    double *X = (double *)malloc(nCount * dCount * sizeof(double));

    int currPos = -1;

    char *line = (char *)malloc(1024 * 1024);
    for (int skip = 0; skip < 4; skip++)
    {
        int got = fscanf(stream, "%s\n", line);
    }
    free(line);

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

    int temp;
    for (int i = 0; i < 2; i++)
    {
        int got = fscanf(stream, "%d", &temp);
    }

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
    int dCount = 17;

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

    fseek(stream,1,SEEK_SET);

    for (int i = 0; i < nCount; i++)
    {
        for (int j = 0; j < dCount; j++)
        {
            if (fscanf(stream, "%d:%lf ", &intValue, &doubleValue) != 2)
            {
                continue;
            }
            counter++;
            X[counter] = doubleValue;
        }
        while (1)
        {
            if (fscanf(stream, "%d:%lf", &intValue, &doubleValue) != 2)
            {
                break;
            }
        }
    }

    fclose(stream);

    printf("counter=%d\n",counter);

    *n = nCount;
    *d = dCount;

    return X;
}
