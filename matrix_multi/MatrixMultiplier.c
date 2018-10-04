#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/*
    A program ket file-al dolgozik.

    A bemeneti fajlokban tarolt adat terbeli elhelyeszkedesere epit a program,
    ez alapjan az oszlop vektor oszlopba rendezett adatokat jelent.

    A program c99 szabvany szerint keszult.
*/

void readArg(int argc);

struct tensor {

    double* body;
    int nrow;
    int ncolumn;
};

void init_tensor(struct tensor* ts);

void destruct_tensor(struct tensor* ts);

void makeTensor(FILE*, struct tensor* ts);

void printTensor(struct tensor *ts);

void scalScalMult(struct tensor* s1, struct tensor* s2, struct tensor* result);

void scalarMult(struct tensor* scal, struct tensor* nonScal, struct tensor* result);

void dotProduct(struct tensor* vs, struct tensor* vo, struct tensor* result);

void dyadProduct(struct tensor* vo, struct tensor* vs, struct tensor* result);

void elementMulti(struct tensor* vec1, struct tensor* vec2, struct tensor* result);

void matrixMulti(struct tensor* ts1, struct tensor* ts2, struct tensor* result);

void operation(struct tensor* ts1, struct tensor* ts2, struct tensor* result);

///////////////////////////////////////////////////////////////////////////

int main(int argc, char *argv[]){

    readArg(argc);

    FILE* data;

    struct tensor T1;
    struct tensor T2;
    struct tensor Te;

    init_tensor(&T1);
    init_tensor(&T2);
    init_tensor(&Te);

    data = fopen(argv[1],"r");
    makeTensor(data,&T1);

    freopen(argv[2],"r",data);
    makeTensor(data,&T2);
    fclose(data);

    printTensor(&T1);
    printTensor(&T2);

    operation(&T1,&T2,&Te);
    printTensor(&Te);

    destruct_tensor(&T1);
    destruct_tensor(&T2);
    destruct_tensor(&Te);

return 0;
}
///////////////////////////////////////////////////////////////////////////

void readArg(int argc){

    if(argc != 3){

        printf("\n\tNem megfelelo parancssori argumentumok!\n");
        printf("\tProbald ujra igy: <program> <File1> <File2>\n\n");
        exit(0);
    }
}

void init_tensor(struct tensor* ts){

    ts->nrow = 0;
    ts->ncolumn = 0;
}

void destruct_tensor(struct tensor* ts){

    free(ts->body);
}

void makeTensor(FILE* data, struct tensor* ts){

    int preNColumn = 0;
    int len = 0;
    int idx = 0;
    int size = 16;

    char buffer[512];
    char ch;

    ts->body = (double*)malloc(size * sizeof(double));

    if(data == NULL) { // EZ az ellenõrzés ne itt legyen, hanem a megnyitással együtt!

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

//////////////
/*
    while(!feof(data)){

        ch = fgetc(data);

        if(ch!=' ' && ch!='\t' && ch!='\r' && ch!='\n'){

            buffer[idx++] = ch;
        }

        else {

            buffer[idx++] = '\0';

            if(idx > 1){

                if(len >= size){
                    ts->body = (double*)realloc(ts->body, size * 2 * sizeof(double));
                    size *= 2;
                }
                ts->body[len++] = atof(buffer);
                (ts->ncolumn)++;
            }

            if(ch == '\n' && ts->ncolumn > 0){

                if( preNColumn > 0 && ts->ncolumn != preNColumn ){

                    printf("\n Adathiba!\n A sorok elemszama eltero!\n\n");
                    exit(0);
                }
                preNColumn = ts->ncolumn;
                ts->ncolumn = 0;
            }

            idx = 0;
        }
    }
*/
//////////////

    ts->ncolumn = preNColumn;

    if(ts->ncolumn == 0){printf("\n Hiba:  ncolumn = 0!\n\n"); exit(0);}

    ts->nrow = len / ts->ncolumn;

    if(ts->nrow * ts->ncolumn != len){printf("\n Hibas meret!\n\n"); exit(0);}
}

void printTensor(struct tensor *ts){

    int i, j;

    printf("\n  A tensor meret: %dx%d\n",ts->nrow,ts->ncolumn);
    printf("\n  A tensor elemei:\n\n");

    for(i=0; i < ts->nrow; i++){

        for(j=0; j < ts->ncolumn; j++){

            printf("  %.2f", ts->body[(i*ts->ncolumn)+j]);
        }
        printf("\n\n");
    }
    printf("\n");
}

void scalScalMult(struct tensor* s1, struct tensor* s2, struct tensor* result){

    result->nrow = 1;
    result->ncolumn = 1;
    result->body = (double*)malloc(sizeof(double));
    result->body[0] = s1->body[0] * s2->body[0];
}

void scalarMult(struct tensor* scal, struct tensor* nonScal, struct tensor* result){

    int i, j;
    int idx = 0;
    double scalar = scal->body[0];

    result->nrow = nonScal->nrow;
    result->ncolumn = nonScal->ncolumn;
    result->body = (double*)malloc(result->ncolumn * result->nrow * sizeof(double));

    for(i=0; i < result->nrow; i++){

        for(j=0; j < result->ncolumn; j++){

            idx = (i*result->ncolumn)+j;
            result->body[idx] = nonScal->body[idx] * scalar;
        }
    }
}

void dotProduct(struct tensor* vs, struct tensor* vo, struct tensor* result){

    int i;
    double sum = 0;

    if(vs->ncolumn != vo->nrow){printf("\n A ket vektor hossza nem azonos!\n\n");exit(0);}

    result->nrow = 1;
    result->ncolumn = 1;
    result->body = (double*)malloc(sizeof(double));

    for(i=0; i < vs->ncolumn; i++){

        sum += vs->body[i] * vo->body[i];
    }
    result->body[0] = sum;
}

void dyadProduct(struct tensor* vo, struct tensor* vs, struct tensor* result){

    int i, j;

    result->nrow = vo->nrow;
    result->ncolumn = vs->ncolumn;
    result->body = (double*)malloc(result->ncolumn * result->nrow * sizeof(double));

    for(i=0; i < result->nrow; i++){

        for(j=0; j < result->ncolumn; j++){

            result->body[(i*result->ncolumn)+j] = vo->body[i] * vs->body[j];
        }
    }
}

void elementMulti(struct tensor* vec1, struct tensor* vec2, struct tensor* result){

    int i;
    int len = vec1->ncolumn * vec1->nrow;

    if(vec1->nrow != vec2->nrow){printf("\n A ket vektor hossza nem azonos!\n\n");exit(0);}
    if(vec1->ncolumn != vec2->ncolumn){printf("\n A ket vektor hossza nem azonos!\n\n");exit(0);}

    result->nrow = vec1->nrow;
    result->ncolumn = vec1->ncolumn;
    result->body = (double*)malloc(len * sizeof(double));

    for(i=0; i < len; i++){

        result->body[i] = vec1->body[i] * vec2->body[i];
    }
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

                sum += ts1->body[(i*idxMax)+k] * ts2->body[(j*idxMax)+k];
            }
            result->body[(i*result->ncolumn)+j] = sum;
            sum = 0;
        }
    }
}

void operation(struct tensor* ts1, struct tensor* ts2, struct tensor* result){

    if(ts1->nrow==1 && ts1->ncolumn==1 && ts2->nrow==1 && ts2->ncolumn==1){ scalScalMult(ts1,ts2,result); }

    else if(ts1->nrow==1 && ts1->ncolumn==1){ scalarMult(ts1,ts2,result); }

    else if(ts2->nrow==1 && ts2->ncolumn==1){ scalarMult(ts2,ts1,result); }

    else if((ts1->nrow==1 || ts1->ncolumn==1) && (ts2->nrow==1 || ts2->ncolumn==1)){

        if(ts1->nrow==1 && ts2->ncolumn==1){ dotProduct(ts1,ts2,result); }
        else if(ts1->ncolumn==1 && ts2->nrow==1){ dyadProduct(ts1,ts2,result); }
        else {elementMulti(ts1,ts2,result);}
    }
    else { matrixMulti(ts1,ts2,result); }
}
