#include <stdio.h>
#include "matrix.h"

#define PD(D) printf("%s = %lf\n", #D, D)

int main(){

    Matrix m1, m2;

    m1 = allocIdentity(3);
    m2 = allocMatrix(4);

    //setMatrix(m2, 0, 0, 4);setMatrix(m2, 0, 1, 12);setMatrix(m2, 0, 2, -16);
    //setMatrix(m2, 1, 0, 12);setMatrix(m2, 1, 1, 37);setMatrix(m2, 1, 2, -43);
    //setMatrix(m2, 2, 0, -16);setMatrix(m2, 2, 1, -43);setMatrix(m2, 2, 2, 98);


    setMatrix(m2, 0, 0, 1);setMatrix(m2, 0, 1, 2);setMatrix(m2, 0, 2, 1);setMatrix(m2, 0, 3, 0);
    setMatrix(m2, 1, 0, 0);setMatrix(m2, 1, 1, 3);setMatrix(m2, 1, 2, 1);setMatrix(m2, 1, 3, 1);
    setMatrix(m2, 2, 0,-1);setMatrix(m2, 2, 1, 0);setMatrix(m2, 2, 2, 3);setMatrix(m2, 2, 3, 1);
    setMatrix(m2, 3, 0, 3);setMatrix(m2, 3, 1, 1);setMatrix(m2, 3, 2, 2);setMatrix(m2, 3, 3, 0);
    printMatrix(m2);


    setMatrix(m1, 0, 0, 1);setMatrix(m1, 0, 1, 1);setMatrix(m1, 0, 2, 0);
    setMatrix(m1, 1, 0, 0);setMatrix(m1, 1, 1, 1);setMatrix(m1, 1, 2, 1);
    setMatrix(m1, 2, 0, 0);setMatrix(m1, 2, 1, 0); setMatrix(m1, 2, 2, 1);

    double t = determinant(&m1);
    double t2 = laplaceDeterminant(m1.data, 3);
    printf("%lf  %lf \n", t, t2);


    //setMatrix(m1, 0, 0, 1);setMatrix(m1, 0, 1, -1);setMatrix(m1, 0, 2, -1);
    //setMatrix(m1, 1, 0, 3);setMatrix(m1, 1, 1, -4);setMatrix(m1, 1, 2, -2);
    //setMatrix(m1, 2, 0, 2);setMatrix(m1, 2, 1, -3);setMatrix(m1, 2, 2, -2);

    //Matrix result = GramSchmidt(m1);
    //printMatrix(result);

    freeMatrix(&m1);
    freeMatrix(&m2);

    return 0;
}
