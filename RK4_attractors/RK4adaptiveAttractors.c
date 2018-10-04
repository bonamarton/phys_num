#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

void load(double* param, double* var);

typedef double (*df)(double*, double*, double, int); // param, var, t, idxF

double lorenz(double* param, double* var, double t, int idxF);

double aizawa(double* param, double* var, double t, int idxF);

double chanLee(double* param, double* var, double t, int idxF);

double luChan(double* param, double* var, double t, int idxF);

double thomas(double* param, double* var, double t, int idxF);

void RK4basic(double* param, double* var, double* result, int n, double t, double tStep, df f);

void RK4StepAdaptive(double* param, double* var, int n, double* t, double* tStep, double delta0, df f);

void RK4_Adaptive_Attractors(FILE* data, FILE* config, clock_t exTime,
    double* param, double* var, double tStepStart, double tEnd, double delta0, df f);

/////////////////////////////////////////////////////////////////////////////////////////

int main(){

    FILE *data, *config;
    clock_t exTime;

    const double tEnd = 1000;
    const double tStepStart = 0.001;
    const double delta0 = 1e-5;

    double param[21], var[15]; // don't forget resize!
    load(param,var);

    data = fopen("lorenzRK4.txt","w");
    config = fopen("lorenzRK4config.txt","w");
    exTime = clock();
    RK4_Adaptive_Attractors(data, config, exTime, &param[0], &var[0], tStepStart, tEnd, delta0, lorenz);

    data = freopen("aizawaRK4.txt","w",data);
    config = freopen("aizawaRK4config.txt","w",config);
    exTime = clock();
    RK4_Adaptive_Attractors(data, config, exTime, &param[4], &var[3], tStepStart, tEnd, delta0, aizawa);

    data = freopen("chanLeeRK4.txt","w",data);
    config = freopen("chanLeeRK4config.txt","w",config);
    exTime = clock();
    RK4_Adaptive_Attractors(data, config, exTime, &param[11], &var[6], tStepStart, tEnd, delta0, chanLee);

    data = freopen("luChanRK4.txt","w",data);
    config = freopen("luChanRK4config.txt","w",config);
    exTime = clock();
    RK4_Adaptive_Attractors(data, config, exTime, &param[15], &var[9], tStepStart, tEnd, delta0, luChan);

    data = freopen("thomasRK4.txt","w",data);
    config = freopen("thomasRK4config.txt","w",config);
    exTime = clock();
    RK4_Adaptive_Attractors(data, config, exTime, &param[19], &var[12], tStepStart, tEnd, delta0, thomas);

    fclose(data);
    fclose(config);

return 0;
}
/////////////////////////////////////////////////////////////////////////////////////////

void load(double* param, double* var){

    param[0] = 3;       // l.number of variables
    param[1] = 10;      // l.sigma
    param[2] = 28;      // l.rho
    param[3] = 8.0/3.0; // l.beta

    param[4] = 3;       // a.number of variables
    param[5] = 0.25;    // a.epsilon
    param[6] = 0.95;    // a.alpha
    param[7] = 0.6;     // a.gamma
    param[8] = 3.5;     // a.delta
    param[9] = 0.7;     // a.beta
    param[10] = 0.1;    // a.zeta

    param[11] = 3;      // c.number of variables
    param[12] = 5;      // c.alpha
    param[13] = -10;    // c.beta
    param[14] = -0.38;  // c.delta

    param[15] = 3;      // lu.number of variables
    param[16] = -10;    // lu.alpha
    param[17] = -4;     // lu.beta
    param[18] = 18.1;   // lu.delta

    param[19] = 3;      // t.number of variables
    param[20] = 0.19;   // t.beta
    // resize param[]

    var[0] = 0.1;   // l.x
    var[1] = 0.1;   // l.y
    var[2] = 0.1;   // l.z

    var[3] = 0.1;   // a.x
    var[4] = 0.0;   // a.y
    var[5] = 0.0;   // a.z

    var[6] = -2;    // c.x
    var[7] = 0.5;   // c.y
    var[8] = 4.5;   // c.z

    var[9] = 0.1;   // lu.x
    var[10] = 0.1;  // lu.y
    var[11] = 5;    // lu.z

    var[12] = 0.1;  // t.x
    var[13] = 0;    // t.y
    var[14] = 4;    // t.z
    // resize var[]
}

double lorenz(double* param, double* var, double t, int idxF){

    double value = 0;

    switch(idxF){

        default: { printf("Error in ODE df!"); break; }

        case 0 : { value = param[1] * (var[1] - var[0]); break; }

        case 1 : { value = var[0] * (param[2] - var[2]) - var[1]; break; }

        case 2 : { value = (var[0] * var[1]) - (param[3] * var[2]); break; }
    }

return value;
}

double aizawa(double* param, double* var, double t, int idxF){

    double value = 0;

    switch(idxF){

        default: { printf("Error in ODE df!"); break; }

        case 0 : { value = ((var[2]-param[5]) * var[0]) - (param[4]*var[1]); break; }

        case 1 : { value = ((var[2]-param[5]) * var[1]) + (param[4]*var[0]); break; }

        case 2 : { value = param[3] + (param[2] * var[2]) - (pow(var[2],3) / 3.0);
                   value -= (pow(var[0],2) + pow(var[1],2)) * (1 + (param[1] * var[2]));
                   value += param[6] * var[2] * pow(var[0],3);
                   break; }
    }

return value;
}

double chanLee(double* param, double* var, double t, int idxF){

    double value = 0;

    switch(idxF){

        default: { printf("Error in ODE df!"); break; }

        case 0 : { value = (param[1] * var[0]) - (var[1] * var[2]); break; }

        case 1 : { value = (param[2] * var[1]) + (var[0] * var[2]); break; }

        case 2 : { value = (param[3] * var[2]) + (var[0] * var[1] / 3.0); break; }
    }

return value;
}

double luChan(double* param, double* var, double t, int idxF){

    double value = 0;

    switch(idxF){

        default: { printf("Error in ODE df!"); break; }

        case 0 : { value = -(param[1] * param[2] * var[0] / (param[1] + param[2])) - (var[1] * var[2]) + param[3]; break; }

        case 1 : { value = (param[1] * var[1]) + (var[0] * var[2]); break; }

        case 2 : { value = (param[2] * var[2]) + (var[0] * var[1]); break; }
    }

return value;
}

double thomas(double* param, double* var, double t, int idxF){

    double value = 0;

    switch(idxF){

        default: { printf("Error in ODE df!"); break; }

        case 0 : { value = sin(var[1]) - (param[1] * var[0]); break; }

        case 1 : { value = sin(var[2]) - (param[1] * var[1]); break; }

        case 2 : { value = sin(var[0]) - (param[1] * var[2]); break; }
    }

return value;
}

void RK4basic(double* param, double* var, double* result, int n, double t, double tStep, df f){

    int i;
    double k1[n], k2[n], k3[n], k4[n], temp[n];

    for(i=0; i < n; i++){
        k1[i] = tStep * f(param, var, t, i);
    }
    for(i=0; i < n; i++){
        temp[i] = var[i] + 0.5 * k1[i];
    }

    for(i=0; i < n; i++){
        k2[i] = tStep * f(param, temp, t, i);
    }
    for(i=0; i < n; i++){
        temp[i] = var[i] + 0.5 * k2[i];
    }

    for(i=0; i < n; i++){
        k3[i] = tStep * f(param, temp, t, i);
    }
    for(i=0; i < n; i++){
        temp[i] = var[i] + k3[i];
    }

    for(i=0; i < n; i++){
        k4[i] = tStep * f(param, temp, t, i);
    }

    for(i=0; i < n; i++){
        temp[i] = var[i] + ((k1[i] + 2*k2[i] + 2*k3[i] + k4[i]) / 6.0);
    }
    for(i=0; i < n; i++){
        result[i] = temp[i];
    }
}

void RK4StepAdaptive(double* param, double* var, int n, double* t, double* tStep, double delta0, df f){

    int i;
    double result1[n], result2[n], result3[n];
    double delta[n], deltaMax=0;

    RK4basic(param, var, result1, n, *t, *tStep, f);
    RK4basic(param, var, result2, n, *t, *tStep/2.0, f);
    RK4basic(param, result2, result3, n, *t, *tStep/2.0, f);

    for(i=0; i < n; i++){

        var[i] = result3[i];

        delta[i] = fabs(result3[i] - result1[i]);
        if(delta[i] > deltaMax) { deltaMax = delta[i]; }
    }

    *t += *tStep;

    *tStep *= pow( delta0/deltaMax , 0.2 );
}

void RK4_Adaptive_Attractors(FILE* data, FILE* config, clock_t exTime,
double* param, double* var, double tStepStart, double tEnd, double delta0, df f){

    int lines = 0;
    double t = 0;
    double tStep = tStepStart;

    while(t < tEnd){

        fprintf(data, "%f %f %f %f\n", t, var[0], var[1], var[2]);

        RK4StepAdaptive(param, var, (int)param[0], &t, &tStep, delta0, f);

        lines++;
    }

    double Time = (double)(clock()-exTime)/(double)CLOCKS_PER_SEC;

    fprintf(config,"\n delta0 = %f\n", delta0);
    fprintf(config,"\n tStepStart = %f\n", tStepStart);
    fprintf(config,"\n tEnd = %f\n", tEnd);
    fprintf(config,"\n Written lines = %d\n", lines);
    fprintf(config,"\n Execution time = %f sec\n", Time);
}
