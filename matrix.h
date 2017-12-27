#ifndef _MAT_
#define _MAT_


//#define printMatrix(A) printf("%s = \n", #A); print_matrix(A)
#define setMatrix(M,I,J,V) M.data[I][J] = V
#define setMatrixPtr(M,I,J,V) M->data[I][J] = V

typedef struct{
    int dimension; //
    double ** data;
} Matrix;

Matrix allocMatrix(int n);
Matrix allocIdentity(int n);
void freeMatrix(Matrix* m);

// impress�o
void printMatrix(Matrix m);

// Opera��es elementares
void swapLines(Matrix* m, int l1, int l2);
void lineByScalar(Matrix* m, int l1, double s);
void sumLineByScalar(Matrix* m, int l1, int l2, double s);

// triangula��o
void gaussianElimination(Matrix* m);

// decomposi��o LU
Matrix LU(Matrix* m);
Matrix* LDU(Matrix* m);

// ortogonaliza��o de Gram-Schimidt
Matrix GramSchmidt(Matrix m);


// retorna QR sem alterar a original
Matrix* QR(Matrix* m);

// Cholesky�Banachiewicz algorithm
Matrix* cholesky(Matrix* m);
double laplaceDeterminant(double** m, int n);
double determinant(Matrix* m);

#endif // _MAT_
