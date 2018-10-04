#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

/*
    A program c99 szabvany szerint keszult.

    Program futtatasa: <program neve> <valtozok szama> <illesztes fokszama> <file>
*/

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

void TEST(struct tensor *M, struct tensor *invM);

/////////////////// New functions /////////////////////

int inputCheck(struct tensor *ts, int nvariable);

int factor(int n);

struct tensor makeIdxSequence(int nvariable, int degree);

void makeLinearSystem(struct tensor *Data, struct tensor *IS, struct tensor *M, struct tensor *B);

void printCoef(struct tensor *IdxSequence, struct tensor *Coef);

struct tensor predict(struct tensor *Data, struct tensor *IS, struct tensor *Coef);

void writeOutput(FILE* Out, struct tensor *True, struct tensor *Predict);

double loss(struct tensor *True, struct tensor *Predict);

/////////////////////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

    readArg(argc);

    const int nvariable = atoi(argv[1]);
    const int degree = atoi(argv[2]);

    FILE* data;

    struct tensor Data;
    struct tensor M;
    struct tensor B;
    init_tensor(&Data);
    init_tensor(&M);
    init_tensor(&B);

    data = fopen(argv[3],"r");
    makeTensor(data,&Data);
    fclose(data);
    if(inputCheck(&Data,nvariable) == 1){ exit(0); }

    struct tensor IS = makeIdxSequence(nvariable,degree);
    makeLinearSystem(&Data,&IS,&M,&B);

    struct tensor Coef = linearSolve(&M,&B);

    struct tensor Pred = predict(&Data,&IS,&Coef);

    data = fopen("predict.txt","a");
    writeOutput(data,&Data,&Pred);
    fclose(data);

    printCoef(&IS,&Coef);
    printf("\n LOSS = %f\n\n\n",loss(&Data,&Pred));

    destruct_tensor(&Data);
    destruct_tensor(&M);
    destruct_tensor(&B);
    destruct_tensor(&IS);
    destruct_tensor(&Coef);

return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////

void readArg(int argc){

    if(argc != 4){

        printf("\n\n Nem megfelelo parancssori argumentumok!\n");
        printf(" Probald ujra igy: <program neve> <valtozok szama> <illesztes fokszama> <file>\n\n");
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
}

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

    if(fabs(pivotVal) < DBL_MIN){ return 1; }

    pivotVal = 1 / pivotVal;
    multiScalRow(M, idx, pivotVal);
    multiScalRow(B, idx, pivotVal);

return 0;
}

void eliminationStage1(struct tensor *M, struct tensor *B, int idx){

    int i;
    double val = 0;
    double diagVal = M->body[(idx*M->ncolumn)+idx];

    if( (diagVal - 1) <= DBL_EPSILON  ){

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

    if( (diagVal - 1) <= DBL_EPSILON  ){

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

        if(err){
            printf("\n\n A szamolas nem vegezheto el!\n");
            printf("\n A matrix szingularis!\n\n");
            break;
        }

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

////////////////// Beginning of the new code ///////////////////////

int inputCheck(struct tensor *ts, int nvariable){

    if(ts->ncolumn-2 != nvariable){
        printf("\n\n  Az adatfile nem felel meg a bekert valtozok szamanak!\n\n");
        return 1;
    }
return 0;
}

int factor(int n){

    if (n < 0) {
        printf("\n\n  ERROR: Factorial number is negative!\n\n");
        return -1;
    }

    int i, fact = 1;

    for(i=2; i <= n; i++){

        fact *= i;
    }

return fact;
}

struct tensor makeIdxSequence(int nvariable, int degree){

    struct tensor IS;
    init_tensor(&IS);
    IS.ncolumn = nvariable;
    IS.nrow = factor(degree+nvariable) / ( factor(nvariable) * factor(degree) );
    IS.body = (double*)malloc(IS.nrow * IS.ncolumn * sizeof(double));

    int* idxSeq = (int*)calloc(IS.ncolumn, sizeof(int));

    int i, j;
    int idxMax = round(pow(degree+1,nvariable));
    int col, IS_row=0, sumDegree, end=0;

    for(j=0; j < nvariable; j++){

        IS.body[j] = 0;
    }
    IS_row++;

    for(i=0; i < idxMax-1; i++){

        col = 0;
        sumDegree = 0;
        end = 0;

        while( (end == 0) && (col < nvariable) ){

            if (idxSeq[col] < degree){

                idxSeq[col]++;
                end = 1;
            }
            else{

                idxSeq[col] = 0;
                col++;
            }
        }

        for(j=0; j < nvariable; j++){

            sumDegree += idxSeq[j];
        }

        if (sumDegree <= degree){

            for(j=0; j < nvariable; j++){

                IS.body[(IS_row * nvariable) + j] = (double)idxSeq[j];
            }
            IS_row++;
        }
    }

return IS;
}

void makeLinearSystem(struct tensor *Data, struct tensor *IS, struct tensor *M, struct tensor *B){

    int nvariable = IS->ncolumn;
    int N = IS->nrow;

    M->ncolumn = N;
    M->nrow = N;
    M->body = (double*)calloc(N*N, sizeof(double));

    B->ncolumn = 1;
    B->nrow = N;
    B->body = (double*)calloc(N, sizeof(double));

    double* value = (double*)malloc(N * sizeof(double));

    int i, j, k;
    int row=0, rowData=0;
    double absError = 0;
    double b = 0;

    for(k=0; k < Data->nrow; k++){

        rowData = k * Data->ncolumn;

        absError = Data->body[rowData+nvariable+1];

        for(j=0; j < N; j++){

            value[j] = 1;

            for(i=0; i < nvariable; i++){

                value[j] *= pow(Data->body[rowData+i], IS->body[(j*nvariable)+i]);
            }
            value[j] /= absError;
        }

        for(j=0; j < N; j++){

            row = j * N;

            for(i=0; i < N; i++){

                M->body[row+i] += value[j] * value[i];
            }
        }

        b = Data->body[rowData+nvariable] / absError;

        for(i=0; i < N; i++){

            B->body[i] += value[i] * b;
        }
    }

    free(value);
}

void printCoef(struct tensor *IdxSequence, struct tensor *Coef){

    if(IdxSequence->nrow != Coef->nrow) { printf("Different datasets in printCoef()"); exit(0); }

    int i, j;

    printf("\n\n Az illesztes parameterei:\n\n\n\n");

    for(i=0; i < IdxSequence->ncolumn; i++){
        printf(" x%d",i+1);
    }
    printf("\n\n\n");

    for(i=0; i < Coef->nrow; i++){

        for(j=0; j < IdxSequence->ncolumn; j++){

            printf("  %d", (int)IdxSequence->body[(i*IdxSequence->ncolumn)+j]);
        }
        printf(" :   %f\n\n", Coef->body[i]);
    }
    printf("\n\n");
}

struct tensor predict(struct tensor *Data, struct tensor *IS, struct tensor *Coef){

    if(IS->nrow != Coef->nrow) { printf("Different datasets in predict()"); exit(0); }

    int nvariable = IS->ncolumn;
    int N = IS->nrow;

    struct tensor Pred;
    init_tensor(&Pred);
    Pred.nrow = Data->nrow;
    Pred.ncolumn = 1;
    Pred.body = (double*)calloc(Pred.nrow, sizeof(double));

    int i, j, k;
    int rowData = 0;
    double value = 1;

    for(k=0; k < Data->nrow; k++){

        rowData = k * Data->ncolumn;

        for(j=0; j < N; j++){

            for(i=0; i < nvariable; i++){

                value *= pow(Data->body[rowData+i], IS->body[(j*nvariable)+i]);
            }

            Pred.body[k] += Coef->body[j] * value;

            value = 1;
        }
    }

return Pred;
}

void writeOutput(FILE* Out, struct tensor *True, struct tensor *Predict){

    if(True->nrow != Predict->nrow) { printf("Different datasets in writeOutput()"); exit(0); }

    int i;
    int col = True->ncolumn-2;

    for(i=0; i < True->nrow; i++){

        fprintf(Out, "%f %f\n", True->body[(i*True->ncolumn)+col], Predict->body[i]);
    }
}

double loss(struct tensor *True, struct tensor *Predict){

    int i;
    int col = True->ncolumn-2;
    double diff = 0;

    for(i=0; i < True->nrow; i++){

        diff += fabs( True->body[(i*True->ncolumn)+col] - Predict->body[i] ) / True->body[(i*True->ncolumn)+col+1];
    }

return diff / 1000.0;
}
