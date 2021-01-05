#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "tester_helper.h"

void checkResult(knnresult result, double* corpus, double* query, int n, int m, int d, int k){

    int isValidR = validateResult(result, corpus, query, n, m, d, k, ROWMAJOR);

    printf("Tester validation: %s NEIGHBORS\n", STR_CORRECT_WRONG[isValidR]);

    if(!isValidR){
        exit(-1);
    }
}