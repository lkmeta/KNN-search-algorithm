#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <assert.h>

#include "tester_helper2.h"

/**
 * Function that checks if the results are correct and prints an according message. If the result is wrong
 * we get an error exit code. It is based on the tester functions given by the professors.
 * Input:
 *      knnresult result: the knn structure created by the method applied to our data
 *      double* corpus: matrix containing the corpus set points
 *      double* query: matrix containing the query set points
 *      int n: number of corpus points
 *      int m: number of query points
 *      int d: number of dimensions
 *      int k: number of nearest neighbors we are looking for
 * Output:
 *      None
**/
void checkResult(knnresult result, double* corpus, double* query, int n, int m, int d, int k){

    //we store the data in row major format, so there is no point checking the col major edition
    int isValidR = validateResult(result, corpus, query, n, m, d, k, ROWMAJOR);

    printf("Tester validation: %s NEIGHBORS\n", STR_CORRECT_WRONG[isValidR]);

    if(!isValidR){
        exit(-1);
    }
}