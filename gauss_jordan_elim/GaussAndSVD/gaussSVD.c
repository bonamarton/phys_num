#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <float.h>

/*
    A program c99 szabvany szerint keszult.

    A mukodeshez szuksegesek a Lapacke header, dll es lib fajlok!!!

    Random generalt szingularis matrixok gyartasa:
        A "randGen = 1" ertek beallitasaval lehet elerni.
        A fuggveny ket parametere:
            "N": meret
            "Const": numerikus szingularitas beallitasa
*/

#include "lapacke.h"

void readArg(int argc);

struct tensor {

    double* body;
    int nrow;
    int ncolumn;
    int* swappedCols;
};

void init_tensor(struct tensor* ts);

void destruct_tensor(struct tensor* ts);

void makeTensor(FILE*, struct tensor* ts);

struct tensor copyTensor(struct tensor *M);

void printTensor(struct tensor *ts);

void matrixMulti(struct tensor* ts1, struct tensor* ts2, struct tensor* result);

int quadraticCheck(struct tensor *ts);

void makeUnitMatrix(struct tensor *ts, int N);

double maxAbsValRow(struct tensor *ts, int row);

int maxAbsValSubMatrix(struct tensor *ts, int row, int column);

void multiScalRow(struct tensor *ts, int row, double scal);

void swapRow(struct tensor *ts, int row1, int row2);

void subRow(struct tensor *ts, int row, int subRow, double coef);

void swapColumn(struct tensor *ts, int col1, int col2);

void preCond(struct tensor *M, struct tensor *B);

void pivot(struct tensor *M, struct tensor *B, int idx);

int normWithPivotValue(struct tensor *M, struct tensor *B, int idx);

void eliminationStage1(struct tensor *M, struct tensor *B, int idx);

void eliminationStage2(struct tensor *M, struct tensor *B, int idx);

void backSwap(struct tensor *M, struct tensor *B);

struct tensor linearSolve(struct tensor *M, struct tensor *B);

double norm1(struct tensor *ts);

double normInf(struct tensor *ts);

double cond1(struct tensor *M, struct tensor *invM);

double condInf(struct tensor *M, struct tensor *invM);

double errorUnitMatrix(struct tensor *UM);

void quadraticTranspose(struct tensor *M);

struct tensor SVD(struct tensor *M);

void TEST(struct tensor *M, struct tensor *invM);

struct tensor read2(FILE* data);

struct tensor read1(FILE* data);

struct tensor randomMake(int N, int Const);

void GaussVsSVD(FILE* dataOut, int maxDim);

///////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

    readArg(argc);

    FILE* data;

    struct tensor B;
    init_tensor(&B);
    struct tensor M;

    int randGen = 0;

    if(randGen){

        M = randomMake(100,10000);
    }

    else{

        init_tensor(&M);

        data = fopen(argv[1],"r");

        makeTensor(data,&M);

        fclose(data);
        if(quadraticCheck(&M) == 1){ exit(0); }
    }

    clock_t ido;
    ido = clock();

    makeUnitMatrix(&B,M.nrow);
    struct tensor result = linearSolve(&M,&B);

    double ido1 = (double)(clock()-ido)/(double)CLOCKS_PER_SEC;
    ido = clock();

    struct tensor resultSVD = SVD(&M);

    double ido2 = (double)(clock()-ido)/(double)CLOCKS_PER_SEC;

    printf("\n #################### Gauss-Jordan #################\n\n");
    TEST(&M, &result);
    printf("\n ######################### SVD #####################\n\n");
    TEST(&M, &resultSVD);

    printf("\n Gauss-Jordan futasi ido: %f sec\n", ido1);
    printf("\n SVD futasi ido: %f sec\n\n\n\n", ido2);

    //// Process Time:  ///////////////////
    /*
     data = fopen("gaussSVD.txt","a");
     GaussVsSVD(data,3000);
     fclose(data);
    */
    ///////////////////////////////////////

    destruct_tensor(&M);
    destruct_tensor(&B);
    destruct_tensor(&result);
    destruct_tensor(&resultSVD);

return 0;
}
///////////////////////////////////////////////////////////////////////////

void readArg(int argc){

    if(argc != 2){

        printf("\n\tNem megfelelo parancssori argumentumok!\n");
        printf("\tProbald ujra igy: <program> <File>\n\n");
        exit(0);
    }
}

void init_tensor(struct tensor* ts){

    ts->nrow = 0;
    ts->ncolumn = 0;
    ts->swappedCols = (int*)calloc(1, sizeof(int));
}

void destruct_tensor(struct tensor* ts){

    ts->nrow = 0;
    ts->ncolumn = 0;
    free(ts->body);
    free(ts->swappedCols);
}

void makeTensor(FILE* data, struct tensor* ts){

    int preNColumn = 0;
    int len = 0;
    int idx = 0;
    int size = 16;

    char buffer[512];
    char ch;

    ts->body = (double*)malloc(size * sizeof(double));

    if(data == NULL) {

        printf("\n Az adatfilet nem sikerult megnyitni!\n\n");
        exit(0);
    }

    while(!feof(data)){

        idx = 0;

        while( fscanf(data,"%c",&ch)==1 && ch!=' ' && ch!='\t' && ch!='\r' && ch!='\n'){

            buffer[idx++] = ch;
        }
        buffer[idx++] = '\0';

        if(idx > 1){

            if(len >= size){
                ts->body = (double*)realloc(ts->body, size * 2 * sizeof(double));
                size *= 2;
            }
            ts->body[len++] = atof(buffer);
            ts->ncolumn++;
        }

        if(ch == '\n' && ts->ncolumn > 0){

            if( preNColumn > 0 && ts->ncolumn != preNColumn ){

                printf("\n Adathiba!\n A sorok elemszama eltero!\n\n");
                exit(0);
            }
            preNColumn = ts->ncolumn;
            ts->ncolumn = 0;
        }
    }

    ts->ncolumn = preNColumn;

    if(ts->ncolumn == 0){printf("\n Hiba:  ncolumn = 0!\n\n"); exit(0);}

    ts->nrow = len / ts->ncolumn;

    if(ts->nrow * ts->ncolumn != len){printf("\n Hibas meret!\n\n"); exit(0);}
}

struct tensor copyTensor(struct tensor *M){

    struct tensor Mc;
    init_tensor(&Mc);
    Mc.nrow = M->nrow;
    Mc.ncolumn = M->ncolumn;
    Mc.body = (double*)malloc(Mc.nrow * Mc.ncolumn * sizeof(double));
    memcpy(Mc.body, M->body, (Mc.nrow * Mc.ncolumn * sizeof(double)));

return Mc;
};

void printTensor(struct tensor *ts){

    int i, j;

    if(ts->nrow==0 || ts->ncolumn==0){
        printf("\n  A tensor ures, merete: %dx%d\n\n",ts->nrow,ts->ncolumn);
        return;
    }
    printf("\n  A tensor meret: %dx%d\n",ts->nrow,ts->ncolumn);
    printf("\n  A tensor elemei:\n\n");

    for(i=0; i < ts->nrow; i++){

        for(j=0; j < ts->ncolumn; j++){

            printf("  %.3f", ts->body[(i*ts->ncolumn)+j]);
        }
        printf("\n\n");
    }
    printf("\n");
}

void matrixMulti(struct tensor* ts1, struct tensor* ts2, struct tensor* result){

    int i, j, k;
    int idxMax = ts1->ncolumn;
    double sum = 0;

    if(ts1->ncolumn != ts2->nrow) {printf("\n A muvelet nem hajthato vegre!\n\n");exit(0);}

    result->nrow = ts1->nrow;
    result->ncolumn = ts2->ncolumn;

    result->body = (double*)malloc(result->ncolumn * result->nrow * sizeof(double));

    for(i=0; i < result->nrow; i++){

        for(j=0; j < result->ncolumn; j++){

            for(k=0; k < idxMax; k++){

                sum += ts1->body[(i*ts1->ncolumn)+k] * ts2->body[(k*ts2->ncolumn)+j];
            }
            result->body[(i*result->ncolumn)+j] = sum;
            sum = 0;
        }
    }
}

int quadraticCheck(struct tensor *ts){

    if(ts->nrow != ts->ncolumn){
        printf("\n  A matrix nem negyzetes, merete: %dx%d\n\n",ts->nrow,ts->ncolumn);
        return 1;
    }
return 0;
}

void makeUnitMatrix(struct tensor *ts, int N){

    int i;

    ts->body = (double*)calloc(N*N, sizeof(double));

    for(i=0; i < N; i++){

        ts->body[(i*N)+i] = 1;
    }

    ts->nrow = N;
    ts->ncolumn = N;
}

double maxAbsValRow(struct tensor *ts, int row){

    if(row >= ts->nrow){
        printf("\n Hiba! A legnagyobb elem megtalalasanal!\n");
        printf("\n A keresett sor kilog a tartomanybol!\n");
        printf("\n Kereses itt: %d Sorok szama: %d\n\n",row+1,ts->nrow);
        return -1;
    }

    int i;
    double val = 0;
    double maxVal = 0;

    for(i=0; i < ts->ncolumn; i++){

        val = fabs( ts->body[(row*ts->ncolumn)+i] );

        if( val > maxVal){ maxVal = val;}
    }

return maxVal;
}

int maxAbsValSubMatrix(struct tensor *ts, int row, int column){

    if(row >= ts->nrow){
        printf("\n Hiba! A legnagyobb elem megtalalasanal!\n");
        printf("\n A keresett tartomany nem megfelelo!\n");
        printf("\n Meret: %dx%d Kereses itt: %dx%d\n\n",ts->nrow,ts->ncolumn,row+1,column+1);
        return -1;
    }

    int i, j;
    double val = 0;
    double maxVal = -1;
    int maxValIdx = 0;

    for(i=row; i < ts->nrow; i++){

        for(j=column; j < ts->ncolumn; j++){

            val = fabs( ts->body[(j*ts->ncolumn)+i] );

            if( val > maxVal){

                maxVal = val;
                maxValIdx = (j*ts->ncolumn) + i;
            }
        }
    }

return maxValIdx;
}

void multiScalRow(struct tensor *ts, int row, double scal){

    if(row >= ts->nrow){
        printf("\n Hiba! A keresett sor kilog a tartomanybol!\n");
        printf("\n Hozzaferes ide: %d Sorok szama: %d\n\n",row+1,ts->nrow);
        exit(0);
    }

    int i;

    for(i=0; i < ts->ncolumn; i++){

        ts->body[(row*ts->ncolumn)+i] *= scal;
    }
}

void swapRow(struct tensor *ts, int row1, int row2){

    if(row1 >= ts->nrow || row2 >= ts->nrow){
        printf("\n Hiba! A keresett sor kilog a tartomanybol!\n");
        printf("\n Hozzaferes ide: %d es %d Sorok szama: %d\n\n",row1+1,row2+1,ts->nrow);
        exit(0);
    }

    int i;
    int idx1 = 0;
    int idx2 = 0;
    double buff = 0;

    for(i=0; i < ts->ncolumn; i++){

        idx1 = (row1*ts->ncolumn)+i;
        idx2 = (row2*ts->ncolumn)+i;

        buff = ts->body[idx1];
        ts->body[idx1] = ts->body[idx2];
        ts->body[idx2] = buff;
    }
}

void subRow(struct tensor *ts, int row, int subRow, double coef){

    if(row >= ts->nrow || subRow >= ts->nrow){
        printf("\n Hiba! A keresett sor kilog a tartomanybol!\n");
        printf("\n Hozzaferes ide: %d es %d Sorok szama: %d\n\n",row+1,subRow+1,ts->nrow);
        exit(0);
    }

    int i;

    for(i=0; i < ts->ncolumn; i++){

        ts->body[(row*ts->ncolumn)+i] -= ts->body[(subRow*ts->ncolumn)+i] * coef;
    }
}

void swapColumn(struct tensor *ts, int col1, int col2){

    if(col1 >= ts->ncolumn || col2 >= ts->ncolumn){
        printf("\n Hiba! A keresett oszlop kilog a tartomanybol!\n");
        printf("\n Hozzaferes ide: %d es %d Oszlopok szama: %d\n\n",col1+1,col2+1,ts->ncolumn);
        exit(0);
    }

    int i;
    int idx1 = 0;
    int idx2 = 0;
    double buff = 0;

    for(i=0; i < ts->nrow; i++){

        idx1 = (i*ts->ncolumn)+col1;
        idx2 = (i*ts->ncolumn)+col2;

        buff = ts->body[idx1];
        ts->body[idx1] = ts->body[idx2];
        ts->body[idx2] = buff;
    }

    buff = ts->swappedCols[col1];
    ts->swappedCols[col1] = ts->swappedCols[col2];
    ts->swappedCols[col2] = buff;
}

void preCond(struct tensor *M, struct tensor *B){

    int i;
    double scal = 0;

    for(i=0; i < M->nrow; i++){

        scal = maxAbsValRow(M,i);

        if(fabs(scal) < DBL_MIN){
            printf("\n\n A szamolas nem vegezheto el!\n");
            printf("\n A matrix szingularis!\n\n");
            exit(0);
        }

        scal = 1 / scal;
        multiScalRow(M, i, scal);
        multiScalRow(B, i, scal);
    }
}

void pivot(struct tensor *M, struct tensor *B, int idx){

    int pivotIdx = maxAbsValSubMatrix(M, idx, idx);
    int pivotRow = (int)floor( (double)pivotIdx / (double)M->ncolumn );
    int pivotCol = pivotIdx - (pivotRow * M->ncolumn);

    if(pivotRow != idx){
        swapRow(M, idx, pivotRow);
        swapRow(B, idx, pivotRow);
    }

    if(pivotCol != idx){
        swapColumn(M, idx, pivotCol);
    }
}

int normWithPivotValue(struct tensor *M, struct tensor *B, int idx){

    double pivotVal = M->body[(idx*M->ncolumn)+idx];

    if(fabs(pivotVal) < DBL_MIN){
        printf("\n\n A szamolas nem vegezheto el!\n");
        printf("\n A matrix szingularis!\n\n");
        return 1;
    }

    pivotVal = 1 / pivotVal;
    multiScalRow(M, idx, pivotVal);
    multiScalRow(B, idx, pivotVal);

return 0;
}

void eliminationStage1(struct tensor *M, struct tensor *B, int idx){

    int i;
    double val = 0;
    double diagVal = M->body[(idx*M->ncolumn)+idx];

    if( (diagVal - 1.0) < DBL_EPSILON ){

        for(i=idx+1; i < M->nrow; i++){

            val = M->body[(i*M->ncolumn)+idx];

            if(fabs(val) < DBL_MIN){continue;}

            subRow(M, i, idx, val);
            subRow(B, i, idx, val);
        }
    }

    else{

        for(i=idx+1; i < M->nrow; i++){

            val = M->body[(i*M->ncolumn)+idx];

            if(fabs(val) < DBL_MIN){continue;}

            val /= diagVal;

            subRow(M, i, idx, val);
            subRow(B, i, idx, val);
        }
    }
}

void eliminationStage2(struct tensor *M, struct tensor *B, int idx){

    int i;
    double val = 0;
    double diagVal = M->body[(idx*M->ncolumn)+idx];

    if( (diagVal - 1.0) < DBL_EPSILON ){

        for(i=idx-1; i > -1; i--){

            val = M->body[(i*M->ncolumn)+idx];

            if(fabs(val) < DBL_MIN){continue;}

            subRow(M, i, idx, val);
            subRow(B, i, idx, val);
        }
    }

    else{

        for(i=idx-1; i > -1; i--){

            val = M->body[(i*M->ncolumn)+idx];

            if(fabs(val) < DBL_MIN){continue;}

            val /= diagVal;

            subRow(M, i, idx, val);
            subRow(B, i, idx, val);
        }
    }
}

void backSwap(struct tensor *M, struct tensor *B){

    int i, j;
    int row1 = 0;
    int row2 = 0;

    free(M->body);
    M->body = (double*)malloc(B->nrow * B->ncolumn * sizeof(double));
    M->nrow = B->nrow;
    M->ncolumn = B->ncolumn;

    for(i=0; i < M->nrow; i++){

        row1 = M->swappedCols[i] * M->ncolumn;
        row2 = i * M->ncolumn;

        for(j=0; j < M->ncolumn; j++){

            M->body[row1+j] = B->body[row2+j];
        }
    }
}

struct tensor linearSolve(struct tensor *M, struct tensor *B){

    struct tensor Mc = copyTensor(M);
    struct tensor Bc = copyTensor(B);

    int i;
    int err = 0;

    Mc.swappedCols = (int*)realloc(Mc.swappedCols, Mc.ncolumn * sizeof(int));
    for(i=0; i < Mc.ncolumn; i++){ Mc.swappedCols[i] = i; }

    preCond(&Mc,&Bc);

    for(i=0; i < Mc.nrow; i++){

        pivot(&Mc, &Bc, i);

        err = normWithPivotValue(&Mc, &Bc, i);
        if(err){ break; }

        eliminationStage1(&Mc, &Bc, i);
    }

    if(err){

        destruct_tensor(&Bc);
        destruct_tensor(&Mc);
        init_tensor(&Mc);
        Mc.body = (double*)malloc(1* sizeof(double));
    }

    else{

        for(i=Mc.nrow-1; i > -1; i--){

            eliminationStage2(&Mc, &Bc, i);
        }

        backSwap(&Mc,&Bc);

        destruct_tensor(&Bc);

    }

return Mc;
}

double norm1(struct tensor *ts){

    int i, j;
    double sumCol = 0;
    double maxSumCol = 0;

    for(i=0; i < ts->ncolumn; i++){

        for(j=0; j < ts->nrow; j++){

            sumCol += fabs( ts->body[(j*ts->ncolumn)+i] );
        }

        if( sumCol > maxSumCol){

            maxSumCol = sumCol;
        }

        sumCol = 0;
    }

return maxSumCol;
}

double normInf(struct tensor *ts){

    int i, j;
    double sumRow = 0;
    double maxSumRow = 0;

    for(i=0; i < ts->nrow; i++){

        for(j=0; j < ts->ncolumn; j++){

            sumRow += fabs( ts->body[(i*ts->ncolumn)+j] );
        }

        if( sumRow > maxSumRow){

            maxSumRow = sumRow;
        }

        sumRow = 0;
    }

return maxSumRow;
}

double cond1(struct tensor *M, struct tensor *invM){

    double cond = norm1(M) * norm1(invM);

return log10(cond);
}

double condInf(struct tensor *M, struct tensor *invM){

    double cond = normInf(M) * normInf(invM);

return log10(cond);
}

double errorUnitMatrix(struct tensor *UM){

    double err = 0;

    int i, j;
    int row = 0;

    for(i=0; i < UM->nrow; i++){

        row = i*UM->ncolumn;

        for(j=0; j < UM->ncolumn; j++){

            if(i == j){ err += fabs( 1 - UM->body[row+j] ); }
            else { err += fabs( UM->body[row+j] ); }
        }
    }

return log10(err);
}

void quadraticTranspose(struct tensor *M){

    if(quadraticCheck(M) == 1){ exit(0); }

    int i, j;
    double buff = 0;

    for(i=0; i < M->nrow; i++){

        for(j=i+1; j < M->ncolumn; j++){

            buff = M->body[(i*M->ncolumn)+j];
            M->body[(i*M->ncolumn)+j] = M->body[(j*M->ncolumn)+i];
            M->body[(j*M->ncolumn)+i] = buff;
        }
    }
}

struct tensor SVD(struct tensor *M){

    if(quadraticCheck(M) == 1){ exit(0); }
    int N = M->nrow;

    struct tensor Mc = copyTensor(M);
    struct tensor temp;
    init_tensor(&temp);
    struct tensor result;
    init_tensor(&result);

    double* s = (double*)malloc(N * sizeof(double));

    struct tensor U;
    init_tensor(&U);
    U.body = (double*)malloc(N*N * sizeof(double));
    U.nrow = N;
    U.ncolumn = N;

    struct tensor vt;
    init_tensor(&vt);
    vt.body = (double*)malloc(N*N * sizeof(double));
    vt.nrow = N;
    vt.ncolumn = N;

    double* superb = (double*)malloc((N-1) * sizeof(double));

    LAPACKE_dgesvd(LAPACK_ROW_MAJOR, 'A', 'A', N, N, Mc.body, N, s, U.body, N, vt.body, N, superb);

    destruct_tensor(&Mc);
    free(superb);

    struct tensor SMinv;
    init_tensor(&SMinv);
    SMinv.body = (double*)calloc(N*N, sizeof(double));
    SMinv.nrow = N;
    SMinv.ncolumn = N;
    int i;
    for(i=0; i < N; i++){
        if(fabs(s[i]) >= DBL_MIN){
            SMinv.body[(i*N)+i] = 1 / s[i];
        }
    }
    free(s);

    quadraticTranspose(&U);
    quadraticTranspose(&vt);

    matrixMulti(&vt,&SMinv,&temp);
    destruct_tensor(&vt);
    destruct_tensor(&SMinv);

    matrixMulti(&temp,&U,&result);
    destruct_tensor(&temp);
    destruct_tensor(&U);

return result;
}

void TEST(struct tensor *M, struct tensor *invM){

    struct tensor result;
    init_tensor(&result);

    if(M->nrow == 0){ printf("\n A teszt nem vegezheto el!\n\n\n\n"); }

    if(M->nrow < 7){
        printf("\n A matrix:\n\n");
        printTensor(M);
        printf("\n A matrix inverse:\n\n");
        printTensor(invM);
    }

    printf("\n  (log10)Cond1: %.3f (log10)CondInf: %.3f\n\n\n", cond1(M,invM), condInf(M,invM));

    matrixMulti(M, invM, &result);

    if(M->nrow < 7){
        printf("\n A matrix es inverzenek szorzata:\n\n");
        printTensor(&result);
    }

    printf("\n  A egyseg matrixtol valo (log10)elteres: %.3f\n\n", errorUnitMatrix(&result));
    printf(" __________________________________________________________________\n\n\n\n\n");

    destruct_tensor(&result);
}

struct tensor read2(FILE* data){

    if(data == NULL) {

        printf("\n Az adatfilet nem sikerult megnyitni!\n\n");
        exit(0);
    }

    struct tensor ts;
    init_tensor(&ts);

    int row = 0;
    int col = 0;
    int buff = 0;

    double element = 0;

    fscanf(data, "%d %d %d ", &ts.nrow, &ts.ncolumn, &buff);

    if(ts.nrow == 0){printf("\n Hiba:  nrow = 0!\n\n"); exit(0);}
    if(ts.ncolumn == 0){printf("\n Hiba:  ncolumn = 0!\n\n"); exit(0);}

    ts.body = (double*)calloc(ts.nrow * ts.ncolumn, sizeof(double));

    while( fscanf(data,"%d %d %lf ",&row,&col,&element)==3 ){

        ts.body[(row * ts.ncolumn)+ col] = element;
    }

    if(feof(data) == 0) {

        printf("\n Az adatfilet nem sikerult beolvasni!\n\n");
        exit(0);
    }

return ts;
}

struct tensor read1(FILE* data){

    if(data == NULL) {

        printf("\n Az adatfilet nem sikerult megnyitni!\n\n");
        exit(0);
    }

    struct tensor ts;
    init_tensor(&ts);

    int row = 0;
    int col = 0;
    int buff = 0;

    fscanf(data, "%d %d %d ", &ts.nrow, &ts.ncolumn, &buff);

    if(ts.nrow == 0){printf("\n Hiba:  nrow = 0!\n\n"); exit(0);}
    if(ts.ncolumn == 0){printf("\n Hiba:  ncolumn = 0!\n\n"); exit(0);}

    ts.body = (double*)calloc(ts.nrow * ts.ncolumn, sizeof(double));

    while( fscanf(data,"%d %d ",&col,&row)==2 ){

        ts.body[(row * ts.ncolumn)+ col] = (double)col;
    }

    if(feof(data) == 0) {

        printf("\n Az adatfilet nem sikerult beolvasni!\n\n");
        exit(0);
    }

return ts;
}

struct tensor randomMake(int N, int Const){

    if(N < 1){printf("\n Hiba:  N = 0!\n\n"); exit(0);}

    struct tensor M;
    init_tensor(&M);
    M.body = (double*)malloc(N*N * sizeof(double));
    M.nrow = N;
    M.ncolumn = N;

    int i, j;
    double element = 0;

    for(i=0; i < N; i++){

        for(j=0; j < N; j++){

            element = (double)rand() / ((double)RAND_MAX * Const * 100000);
            M.body[(i*N)+j] = element + Const;
            if((rand()%100) < 1){ M.body[(i*N)+j] = 0; }
        }
    }

return M;
}

void GaussVsSVD(FILE* dataOut, int maxDim){

    int i;
    clock_t ido;
    struct tensor M;

    printf("\n\n\n Osszehasonlitas progressbar:\n\n");

    for(i=30; i < maxDim; i++){

        if(i < 1000) {i = round((double)i * 1.05); }
        else if(i < 3000) {i = round((double)i * 1.3); }

        else{   printf("\n\n Az osszehasonlitas leallt N ~= 3000-nel!\n");
                printf(" N > 3000-re sok idot vennne igenybe!\n\n");
                break;
            }

        if (i > maxDim) { break;}

        printf(" %d / %d \n", i, maxDim);

        M = randomMake(i,10000);

        ido = clock();

        struct tensor B;
        init_tensor(&B);
        makeUnitMatrix(&B,M.nrow);
        struct tensor result = linearSolve(&M,&B);

        double ido1 = (double)(clock()-ido)/(double)CLOCKS_PER_SEC;
        ido = clock();

        struct tensor resultSVD = SVD(&M);

        double ido2 = (double)(clock()-ido)/(double)CLOCKS_PER_SEC;

        fprintf(dataOut, "%d %f %f\n", i, ido1, ido2);

        destruct_tensor(&B);
        destruct_tensor(&M);
        destruct_tensor(&result);
        destruct_tensor(&resultSVD);
    }
}
