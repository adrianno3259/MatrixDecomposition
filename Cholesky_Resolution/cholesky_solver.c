#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define setMatrix(M,I,J,V) M.data[I][J] = V
#define setMatrixPtr(M,I,J,V) M->data[I][J] = V
#define PI(I) printf("%s = %d\n", #I, (I))
#define PD(D) printf("%s = %lf\n", #D, (D))
#define INFO(S) printf("%s\n", #S)

typedef struct{
    int dimension; //
    double ** data;
} Matrix;

Matrix allocMatrix(int n){
    Matrix ret; int i;
    ret.dimension = n;
    ret.data = (double**) malloc(n*sizeof(double*));
    for(i = 0; i < n; i++){
        ret.data[i] = (double*) malloc(n*sizeof(double));
    }
    return ret;
}

Matrix allocIdentity(int n){
    Matrix ret; int i, j;
    ret.dimension = n;
    ret.data = (double**) malloc(n*sizeof(double*));
    for(i = 0; i < n; i++){
        ret.data[i] = (double*) malloc(n*sizeof(double));
        for(j = 0; j < n; j++){
            if(i != j)  ret.data[i][j] = 0;
            else        ret.data[i][j] = 1;
        }
    }
    return ret;
}

void freeMatrix(Matrix* m){
    int i;
    for(i = 0; i < m->dimension; i++){
        free(m->data[i]);
    }
    free(m->data);
}

void printMatrix(Matrix m){
    int i, j;
    for(i = 0; i < m.dimension; i++){
        for(j = 0; j < m.dimension; j++){
            printf("%lf ", m.data[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void printSystem(Matrix m, double* b){
    int i, j;
    for(i = 0; i < m.dimension; i++){
        for(j = 0; j < m.dimension; j++){
            if(j > 0) printf("%c ", (m.data[i][j] < 0 ? '-' : '+'));
            printf("%.2lf x", m.data[i][j] < 0 ? -m.data[i][j] : m.data[i][j]);
            printf("%d ",j);
        }
        printf("= %lf\n", b[i]);
    }
    printf("\n");
}

// Operações elementares
void swapLines(Matrix* m, int l1, int l2){
    double *tmp = m->data[l1];
    m->data[l1] = m->data[l2];
    m->data[l2] = tmp;
}

void lineByScalar(Matrix* m, int l1, double s){
    int i;
    for(i = 0; i < m->dimension; i++)
        m->data[l1][i] *= s;
}

void sumLineByScalar(Matrix* m, int l1, int l2, double s){
    int i;
    for(i = 0; i < m->dimension; i++)
        m->data[l1][i] += s*m->data[l2][i];
}

// Cholesky-Banachiewicz algorithm
Matrix* cholesky(Matrix* m){

    // verificação Cholesky
    int i, j, k, n = m->dimension;

    for(i = 0; i < n; i++)
    for(j = i; j < n; j++)
        if((m->data[i][j] != m->data[j][i])){
            fprintf(stderr, "Erro: Matrix nao e simétrica...\n");
            return NULL;
        } else if(i == j && (m->data[i][j] < 0)){
            fprintf(stderr, "Erro: Diagonal nao positiva...\n");
            return NULL;
        }

    Matrix* LL = (Matrix*) malloc(sizeof(Matrix)*2);
    LL[0] = allocIdentity(n);
    LL[1] = allocIdentity(n);

    //printf("teste");

    for(j = 0; j < n; j++){
        //PI(j);
        for(i = j; i < n; i++){
            //PI(i);
            double sum = 0;
            double lij;
            if( i > j){
                for(k = 0; k < j; k++)
                    sum += LL[0].data[i][k]*LL[0].data[j][k];
                lij = (m->data[i][j] - sum)/LL[0].data[j][j];

            } else if (i == j) {
                for(k = 0; k < j; k++)
                    sum += LL[0].data[j][k]*LL[0].data[j][k];
                lij = sqrt(m->data[k][i] - sum);
            }
            //PD(LL[0].data[i][j]);
            LL[0].data[i][j] = LL[1].data[j][i] = lij;
        }
    }
    //printf("teste\n\n");
    return LL;
}

double laplaceDeterminant(double** m, int n){
    if(n == 1){
        double det = m[0][0];
        return det;
    }
    else{

        double det = 0;
        int col, i;
        for(col = 0; col < n; col++){
            if(m[0][col] != 0){
                int j;
                Matrix tmp_mat = allocMatrix(n-1);
                for(i = 1; i < n; i++)
                    for(j = 0; j < n; j++)
                        if(j != col){
                            if(j > col) tmp_mat.data[i-1][j-1] = m[i][j];
                            else tmp_mat.data[i-1][j] = m[i][j];
                        }
                double tmp = laplaceDeterminant(tmp_mat.data, n-1);
                det+= pow(-1, col+1)*m[0][col]*tmp;
                freeMatrix(&tmp_mat);
            }
        }
        return det;
    }
}

double determinant(Matrix* m){
    printMatrix(*m);
    double r = laplaceDeterminant(m->data, m->dimension);
    PD(r);
    return r;
}

void solveCholesky(Matrix* m, double* b){

    int i, j, k, n = m->dimension;
    //printMatrix(*m);
    //printf("b = [ ");
    //for(i = 0; i < n; i++) printf("%lf ", b[i]);
    //printf("]\n\n");

    Matrix* res = cholesky(m);
    PI(res==NULL);
    if(res != NULL){
        printMatrix(*m);
        printMatrix(res[0]);
        printMatrix(res[1]);

        // L * y = b
        // retro-substituição de "cima para baixo"
        double *ys = (double*) malloc(sizeof(double)*m->dimension);
        for(i = 0; i < n; i++){
            double sum = 0;
            for(j = 0; j < i; j++)
                sum += res[0].data[i][j]*ys[j];
            ys[i] = (1/res[0].data[i][i]) * (b[i] - sum);
        }

        for(i = 0; i < n; i++){
            printf("y%d = %lf\n", i, ys[i]);
        }
        printf("\n\n\n");

        // Lt * x = y
        double* xs = (double*) malloc(sizeof(double)*n);
        // retro-substituição
        for(i = n-1; i >= 0; i--){
            double sum = 0;
            for(j = n-1; j > i; j--)
                sum += res[1].data[i][j]*xs[j];
            xs[i] = (1/res[1].data[i][i]) * (ys[i] - sum);
        }

        for(i = 0; i < n; i++){
            printf("x%d = %lf\n", i, xs[i]);
        }
        printf("\n\n");
        free(xs);
        free(ys);
        freeMatrix(&res[0]);
        freeMatrix(&res[1]);
        free(res);
    }
}

int main(){

    Matrix m1;
    double b[3];
    m1 = allocMatrix(3);
    setMatrix(m1, 0, 0, 4);setMatrix(m1, 0, 1, 12);setMatrix(m1, 0, 2,-16); b[0] = 5;
    setMatrix(m1, 1, 0, 12);setMatrix(m1, 1, 1, 37);setMatrix(m1, 1, 2,-43); b[1] = 1;
    setMatrix(m1, 2, 0, -16);setMatrix(m1, 2, 1,-43);setMatrix(m1, 2, 2, 98); b[2] = 10;

    int size, ans;
    printf("Deseja rodar o exemplo ou um sistema especifico(0 - exemplo/ 1 - digitar matriz)? Exemplo:\n\n");
    printSystem(m1,b);
    printf("\nresposta: ");
    scanf("%d", &ans);

    if(!ans)
        solveCholesky(&m1, b);
    else{
        printf("\nQual o tamanho da matriz: ");
        scanf("%d", &size);
        freeMatrix(&m1);
        m1 = allocMatrix(size);
        double tmp[size];
        int i, j;
        for(i = 0; i < size; i++){
            for(j = 0; j < size; j++){
                printf("\nA[%d][%d] = ", i, j);
                scanf("%lf", &m1.data[i][j]);
            }
            printf("\nb[%d] = ", i);
            scanf("%lf", tmp+i);
        }
        printSystem(m1, tmp);
        if(determinant(&m1))
            solveCholesky(&m1, tmp);
        else
            printf("matriz não pode ser escalonada");
    }
    freeMatrix(&m1);
    return 0;
}

