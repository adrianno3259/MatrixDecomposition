#include "matrix.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

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

void transposeMatrix(Matrix* m){
    int i, j;
    for(i = 0; i < m->dimension; i++)
    {
        for(j = 0; j < i; j++){
            double tmp = m->data[i][j];
            m->data[i][j] = m->data[j][i];
            m->data[j][i] = tmp;
        }
    }
}

void printMatrixLine(Matrix m, int i){
    int j;
    for(j = 0; j < m.dimension; j++){
        printf("%lf ", m.data[i][j]);
    }
    printf("\n\n");
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

#define PI(I) printf("%s = %d\n", #I, (I))
#define PD(D) printf("%s = %lf\n", #D, (D))
#define INFO(S) printf("%s\n", #S)

void gaussianElimination(Matrix* m){
    int i, j, col = 0, n = m->dimension;

    for(i = 0; i < n; i++){ // linha por linha
        PI(i);
        // verifica se pivot é 0
        if(m->data[i][col] == 0){
            for(j = i; j < n; j++){
                PI(j);
                if(m->data[j][col] != 0){
                    swapLines(m, i, j);
                }
            }
            if(j == n){ col++; continue;}
        }
        double pivo = m->data[i][col];
        PD(pivo);
        for(j = i+1; j < n; j++){
            double pivo2 = m->data[j][col];
            PD(pivo2);
            sumLineByScalar(m, j, i, -(pivo2/pivo));
            printMatrix(*m);
        }
        col++;
    }
}

// altera a matriz original para U e retorna L
Matrix LU(Matrix* m){
    Matrix L = allocIdentity(m->dimension);
    int i, j, col = 0, n = m->dimension;

    for(i = 0; i < n; i++){ // linha por linha
        //PI(i);
        // verifica se pivot é 0
        if(m->data[i][col] == 0){
            for(j = i; j < n; j++){
                //PI(j);
                if(m->data[j][col] != 0){
                    swapLines(m, i, j);
                }
            }
            if(j == n){ col++; continue;}
        }
        double pivo = m->data[i][col];
        //PD(pivo);
        for(j = i+1; j < n; j++){
            double pivo2 = m->data[j][col];
            //PD(pivo2);
            double mij = (pivo2/pivo);
            setMatrix(L,j,i,mij);
            sumLineByScalar(m, j, i, -mij);
            //printMatrix(*m);
        }
        col++;
    }
    return L;
}

Matrix* LU_NoChanges(Matrix* m){
    Matrix* LU = (Matrix*) malloc(sizeof(Matrix)*2);
    LU[0] =  allocIdentity(m->dimension);
    LU[1] = allocMatrix(m->dimension);
    int i, j, col = 0, n = m->dimension;

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            setMatrix(LU[1], i, j, m->data[i][j]);
        }
    }
    for(i = 0; i < n; i++){ // linha por linha
        //PI(i);
        // verifica se pivot é 0
        if(LU[1].data[i][col] == 0){
            for(j = i; j < n; j++){
                //PI(j);
                if(LU[1].data[j][col] != 0){
                    swapLines(m, i, j);
                }
            }
            if(j == n){ col++; continue;}
        }
        double pivo = LU[1].data[i][col];
        //PD(pivo);
        for(j = i+1; j < n; j++){
            double pivo2 = LU[1].data[j][col];
            //PD(pivo2);
            double mij = (pivo2/pivo);
            setMatrix(LU[0],j,i,mij);
            sumLineByScalar(LU+1, j, i, -mij);
            //printMatrix(*m);
        }
        col++;
    }
    return LU;
}


// altera a matriz original para U e retorna L
Matrix* LDU(Matrix* m){
    Matrix* LD = malloc(2*sizeof(Matrix));
    LD[0] = allocIdentity(m->dimension);
    LD[1] = allocIdentity(m->dimension);
    int i, j, col = 0, n = m->dimension;

    for(i = 0; i < n; i++){ // linha por linha
        //PI(i);
        // verifica se pivot é 0
        if(m->data[i][col] == 0){
            for(j = i; j < n; j++){
                //PI(j);
                if(m->data[j][col] != 0){
                    swapLines(m, i, j);
                }
            }
            if(j == n){ col++; continue;}
        }
        double pivo = m->data[i][col];
        //PD(pivo);
        for(j = i+1; j < n; j++){
            double pivo2 = m->data[j][col];
            //PD(pivo2);
            double mij = (pivo2/pivo);
            setMatrix(LD[0],j,i,mij);
            sumLineByScalar(m, j, i, -mij);
            //printMatrix(*m);
        }
        col++;
    }

    // decompondo U em D e U
    for(i = 0; i < n; i++){
        double piv = m->data[i][i];
        setMatrix(LD[1],i,i,piv);
        lineByScalar(m, i, 1/piv);
    }

    return LD;
}

// NÃO altera a matriz original. Retorna L D e U no array
Matrix* LDU_NoChanges(Matrix* m){
    Matrix* LD = malloc(3*sizeof(Matrix));
    LD[0] = allocIdentity(m->dimension);
    LD[1] = allocIdentity(m->dimension);
    LD[2] = allocIdentity(m->dimension);
    int i, j, col = 0, n = m->dimension;

    for(i = 0; i < n; i++){
        for(j = 0; j < n; j++){
            setMatrix(LD[2], i, j, m->data[i][j]);
        }
    }

    for(i = 0; i < n; i++){ // linha por linha
        // verifica se pivot é 0
        if(LD[2].data[i][col] == 0){
            for(j = i; j < n; j++){
                if(LD[2].data[j][col] != 0){
                    swapLines(LD+2, i, j);
                }
            }
            if(j == n){ col++; continue;}
        }
        double pivo = LD[2].data[i][col];
        //PD(pivo);
        for(j = i+1; j < n; j++){
            double pivo2 = LD[2].data[j][col];
            //PD(pivo2);
            double mij = (pivo2/pivo);
            setMatrix(LD[0],j,i,mij);
            sumLineByScalar(LD+2, j, i, -mij);
        }
        col++;
    }

    // decompondo U em D e U
    for(i = 0; i < n; i++){
        double piv = LD[2].data[i][i];
        setMatrix(LD[1],i,i,piv);
        lineByScalar(LD+2, i, 1/piv);
    }

    return LD;
}

Matrix GramSchmidt(Matrix m){
    printMatrix(m);
    int i, j, k, n = m.dimension;
    Matrix w = allocMatrix(n);
    printMatrix(w);
    for(i = 0; i < n; i++){ // obtendo cada vetor Wi ortogonal
        // definindo o primeiro vetor ortogonal igual
        for(j = 0; j < n; j++){setMatrix(w, i, j, m.data[i][j]);} //wi = vi
        //PI(i);
        //vecPrint("vi", &m.data[i], n);
        //vecPrint("wi", &w.data[i], n);
        //printf("antes: \n");
        //printMatrix(w);

        for(j = 0; j < i; j++){ //cada vetor anterior
            //PI(j);
            // produtos internos <vi,wj> e <wj, wj>
            int col; double res_vw = 0, res_ww = 0;
            for(col = 0; col < n; col++){
                //printf("---------------------------------\n");
                //PI(i); PI(j); PI(col);
                //printMatrix(w);
                //PD(m.data[i][col]);
                //PD(w.data[j][col]);
                //PD(m.data[i][col] * w.data[j][col]);
                //PD(w.data[j][col] * w.data[j][col]);
                res_vw += m.data[i][col] * w.data[j][col],  //  <vi,wi> . vi
                res_ww += w.data[j][col] * w.data[j][col];  //  <wi,wi>
            }

            //double vi_dot_wj = vecInnerProduct(m.data+i, w.data+j, n);
            //double vi2 = vecInnerProduct(w.data+j, w.data+j, n);
            double coef = res_vw/res_ww;
            //PD(res_vw); PD(res_ww); PD(coef);
            for(col = 0; col < n; col++){
                //printf("\n\n");
                //PD(w.data[i][col]); PD(coef); PD(w.data[j][col]);
                setMatrix(w, i, col, w.data[i][col] -(coef*w.data[j][col]));
                //PD(w.data[i][col]);
            }
        }
        //printf("depois: \n");
        //printMatrix(w);
    }

    //printf("finalizada: \n");
    //printMatrix(w);

    // normalização
    for(i = 0; i < n; i++){ // para cada vetor wi
        double size = 0;
        for(j = 0; j < n; j++) size += w.data[i][j]*w.data[i][j];
        //PD(size);
        //PD(sqrt(size));
        size = sqrt(size);
        for(j = 0; j < n; j++) w.data[i][j] /= size ;
    }
    //printf("normalizada: \n");
    //printMatrix(w);
    return w;
}

// retorna QR sem alterar a original
Matrix* QR(Matrix* m){
    transposeMatrix(m);
    int i, j, k, n = m->dimension;
    Matrix* QR = (Matrix*) malloc(sizeof(Matrix)*2);
    QR[0] = allocIdentity(n);
    QR[1] = allocIdentity(n);
    for(i = 0; i < n; i++){ // obtendo cada vetor Wi ortogonal
        // definindo o primeiro vetor ortogonal igual
        for(j = 0; j < n; j++){setMatrix(QR[0], i, j, m->data[i][j]);} //wi = vi
        for(j = 0; j < i; j++){ //cada vetor anterior
            // produtos internos <vi,wj> e <wj, wj>
            int col; double res_vw = 0, res_ww = 0;
            for(col = 0; col < n; col++)
                res_vw += m->data[i][col] * QR[0].data[j][col],  //  <vi,wi> . vi
                res_ww += QR[0].data[j][col] * QR[0].data[j][col];  //  <wi,wi>
            double coef = res_vw/res_ww;
            for(col = 0; col < n; col++)
                setMatrix(QR[0], i, col, QR[0].data[i][col] -(coef*QR[0].data[j][col]));
        }
    }
    // normalização
    for(i = 0; i < n; i++){ // para cada vetor wi
        double size = 0;
        for(j = 0; j < n; j++) size += QR[0].data[i][j]*QR[0].data[i][j];
        size = sqrt(size);
        for(j = 0; j < n; j++) QR[0].data[i][j] /= size ;
    }

    int col;
    // preenchendo R
    for(i = 0; i < n; i++){
        for(j = i; j < n; j++){ // R[i,j] <- <ei, aj>
            double res_ei_aj = 0;
            for(col = 0; col < n; col++)
                res_ei_aj += QR[0].data[i][col]*m->data[j][col];
            setMatrix(QR[1], i, j, res_ei_aj);
        }
    }
    transposeMatrix(m); // retornando a matriz ao estado original
    return QR;
}


// Cholesky–Banachiewicz algorithm
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

    printf("teste");

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
    printf("teste\n\n");
    return LL;
}


double laplaceDeterminant(double** m, int n){
    //PI(n);
    if(n == 1){
        double det = m[0][0];
        //PD(det);
        //printf("-------------end %d----------------\n",n);
        return det;
    }
    else{

        double det = 0;
        int col, i;
        for(col = 0; col < n; col++){
            if(m[0][col] != 0){
                int j;
                Matrix tmp_mat = allocMatrix(n-1);
                for(i = 1; i < n; i++){
                    for(j = 0; j < n; j++){
                        if(j != col){
                            if(j > col){
                                tmp_mat.data[i-1][j-1] = m[i][j];
                            }
                            else{
                                tmp_mat.data[i-1][j] = m[i][j];
                            }
                        }
                    }
                }
                /*
                printf("teste\n");
                printf("matrix: \n");
                for(i = 0; i < n-1; i++){
                    for(j = 0; j < n-1; j++){
                        printf("%lf ", tmp_mat.data[i][j]);
                    }
                    printf("\n");
                }
                printf("------------ %d -----------------\n", n);
                */
                double tmp = laplaceDeterminant(tmp_mat.data, n-1);
                det+= pow(-1, col+1)*m[0][col]*tmp;
                freeMatrix(&tmp_mat);
            }

        }
        //PD(det);
        //printf("-------------end %d----------------\n",n);
        return det;
    }
}

double determinant(Matrix* m){
    printMatrix(*m);
    double r = laplaceDeterminant(m->data, m->dimension);
    PD(r);
    return r;
}
