#include <stdio.h>
#include <stdlib.h>
#include <math.h>
    
double* readASCGZ(char* filename, int* n, int* d){

    FILE *stream = fopen(filename, "r");
    if (stream == NULL)
    {
        printf("Error opening file");
        exit(-1);
    }

    double doubleValue;

    int nCount = 1;
    int dCount = 0;

    int tempCounter=0;
    int currPos=-1;

    //define number of dimensions
    while(fscanf(stream, "%lf ", &doubleValue)){
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) !=0) 
        {
            if(tempCounter==0){
                tempCounter++;
                continue;
            }
            else break;
        }
        else{
            dCount++;
        }
    }

    //define number of points
    rewind(stream);
    while (fscanf(stream, "%lf ", &doubleValue) == 1){
        if (fabs(doubleValue-nCount)<1e-6) 
        {
            nCount++;
        }
    }

    nCount--;

    printf("Dimension is %d\n", dCount);
    printf("Points are %d\n", nCount);

    double* X = (double*)malloc(nCount*dCount*sizeof(double));

    rewind(stream);

    while (fscanf(stream, "%lf ", &doubleValue) == 1)
    {   
        if (floorf(doubleValue) == ceilf(doubleValue) && floorf(doubleValue) !=0) 
        {
            continue;
        }
        else
        {   
            currPos++;
            X[currPos] = 10*doubleValue;
        }

    }
    printf("aaaaaaa");
    fclose(stream);
    printf("bbbbbbbb");
    printf("\n points = %d, total = %d ", nCount, dCount);

    *n = nCount;
    *d = dCount;

    return X;
}