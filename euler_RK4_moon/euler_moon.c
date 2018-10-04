#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

double df(double* param, double* var, int Switch, int tStep);

void eulerStep(double* param, double* var, int nVar, int tStep);

double energy(double* param, double* var);

int printIdx(int tStep, int tEnd, int N);

/////////////////////////////////////////////////////////////////////////////////////////

int main(){

    clock_t exTime;
    exTime = clock();

    const int tEnd = 10000000;
    const int tStep = 1000;
    const int nVar = 4;

    double param[11];
    param[0] = 5.9736e24;     // massEarth
    param[1] = 7.349e22;      // massMoon
    param[2] = 4055e5;        // distApo
    param[3] = 3633e5;        // distPeri
    param[4] = 964;           // velApo
    param[5] = 1076;          // velPeri
    param[6] = -6.67384e-11;  // -G
    param[7] = param[6] * param[0]; // const1
    param[8] = param[7] * param[1]; // const2
    param[9] = 0;             // distance
    param[10] = 0;            // modified distance

    double var[2*nVar];   // elso oszlop az aktualis ertek, masoik az eulerban kiszamlot uj lesz
    var[0] = 0;      // x
    var[1] = 3633e5; // y
    var[2] = 1076;   // v_x
    var[3] = 0;      // v_y

    FILE* data;
    data = fopen("moonEuler.txt","w");

    int i, p=0;
    int print = printIdx(tStep,tEnd,100000);

    for(i=0; i < tEnd; i+=tStep){

        param[9] = pow((pow(var[0],2) + pow(var[1],2)), 0.5);
        param[10] = pow(param[9], 3);

        if(++p == print){

            fprintf(data, "%d %f %f %f\n", i, var[0], var[1], energy(param,var));
            p = 0;
        }

        eulerStep(param, var, nVar, tStep);
    }

    data = freopen("moonEulerConfig.txt","w",data);

    double Time = (double)(clock()-exTime)/(double)CLOCKS_PER_SEC;

    fprintf(data,"\n tStep = %d\n", tStep);
    fprintf(data,"\n tEnd = %d\n", tEnd);
    fprintf(data,"\n Execution time = %f sec\n", Time);

    fclose(data);

return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////

double df(double* param, double* var, int f, int tStep){

    double value = 0;

    switch(f){

        default: { printf("Hiba! Nem letezo switch case!"); break; }

        case 0 : { value = var[2] + tStep * param[7] * var[0] / (2 * param[10]); break; } // modification in progress !!!!!!!!

        case 1 : { value = var[3] + tStep * param[7] * var[1] / (2 * param[10]); break; }

        case 2 : { value = param[7] * var[0] / param[10]; break; }

        case 3 : { value = param[7] * var[1] / param[10]; break; }
    }

return value;
}

void eulerStep(double* param, double* var, int nVar, int tStep){

    int i;

    for(i=0; i < nVar; i++){

        var[i+nVar] = var[i] + (tStep * df(param, var, i, tStep));
    }

    for(i=0; i < nVar; i++){

        var[i] = var[i+nVar];
    }
}

double energy(double* param, double* var){

return (param[8] / param[9]) + (0.5 * param[1] * (pow(var[2],2) + pow(var[3],2)));
}

int printIdx(int tStep, int tEnd, int N){

    if(N >= (tEnd / tStep)){ return 1; }

    else { return tEnd / tStep / N; }
}
