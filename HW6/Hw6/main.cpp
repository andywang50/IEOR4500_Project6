#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <algorithm>

int solve(double, double, int, int);

int main(int argc, char **argv)
{
    int N, T;
    int retcode = 0;
    double alpha, pi;

    
    if(argc != 5){
        printf("use: trade alpha pi N T\n"); retcode = 1; return retcode;
    }
    
    alpha = atof(argv[1]);
    pi = atof(argv[2]);
    N = atoi(argv[3]);
    T = atoi(argv[4]);

    printf("alpha = %g pi = %g N = %d T =%d\n", alpha, pi, N, T);
    
    retcode = solve(alpha, pi, N, T);
    

    
    printf("\nran with code %d\n", retcode);
    return retcode;
}

int solve(double alpha, double pi, int N, int T){
    int retcode = 0;
    int j,t;
    double *optimal0, *optimal1, candidate, bestone0, bestone1, *shift0, *shift1;
    int *execution0, *execution1, *best_executions;
    int bestk0, bestk1;
    double *prices;
    optimal0 = (double *)calloc((N + 1)*T, sizeof(double));
    optimal1 = (double *)calloc((N + 1)*T, sizeof(double));
    execution0 = (int *)calloc((N + 1)*T, sizeof(int));
    execution1 = (int *)calloc((N + 1)*T, sizeof(int));
    
    
    if (!optimal0 || !optimal1 || !execution0 || !execution1){
        printf("cannot allocate large matrix\n"); retcode = 2; return retcode;
    }
    
    
    shift0 = (double *)calloc(N + 1, sizeof(double));
    shift1 = (double *)calloc(N + 1, sizeof(double));
    
    if (!shift0 || !shift1){
        printf("cannot allocate large matrix\n"); retcode = 2; return retcode;
    }
    
    prices = (double *)calloc(T, sizeof(double));
    best_executions = (int *)calloc(T, sizeof(int));
    
    if (!prices || !best_executions){
        printf("cannot allocate large matrix\n"); retcode = 2; return retcode;
    }
    
    for (j = 0; j <= N; j++){
        if (j<=99){
            shift0[j] = 1;
            shift1[j] = 1;
        }
        else if (j <= 900){
            shift0[j] = 1 - alpha*pow((double)log(j), pi);
            shift1[j] = shift0[j];
        }
        else{
            shift0[j] = 1 - alpha*pow((double)log(j), 2*pi);
            shift1[j] = 1 - alpha*pow((double)log(j), 3*pi);
        }
    }
	   
    
    /** do last stage **/
    for (j = 0; j <= N; j++){
        if (j <= 99){
            optimal0[(T-1)*(N+1)+j] = j;
            optimal1[(T-1)*(N+1)+j] = j;
            execution0[(T-1)*(N+1)+j] = j;
            execution1[(T-1)*(N+1)+j] = j;
            
        }
        else if (j <= 900){
            bestone0 = 99;
            bestk0 = 99;
            bestk1 = 99;
            for (int k=1;k<=std::min(9,int(floor(j/100.0)));++k){
                candidate = shift0[100*k]*100*k;
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = 100*k;
                    bestk1 = 100*k;
                }
            }
            optimal0[(T-1)*(N+1)+j] = bestone0;
            optimal1[(T-1)*(N+1)+j] = bestone0;
            execution0[(T-1)*(N+1)+j] = bestk0;
            execution1[(T-1)*(N+1)+j] = bestk1;
            
        }
        else{
            bestone0 = 99;
            bestone1 = 99;
            bestk0 = 99;
            bestk1 = 99;
            for (int k=1;k<=9;++k){
                candidate = shift0[100*k]*100*k;
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = 100*k;
                }
                if (candidate > bestone1){
                    bestone1 = candidate;
                    bestk1 = 100*k;
                }
            }
            for (int k=1; k<=int(floor(j/1000.0));++k){
                candidate = shift0[1000*k]*1000*k;
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = 1000*k;
                }
                candidate = shift1[1000*k]*1000*k;
                if (candidate > bestone1){
                    bestone1 = candidate;
                    bestk1 = 1000*k;
                    
                }
            }
            optimal0[(T-1)*(N+1)+j] = bestone0;
            optimal1[(T-1)*(N+1)+j] = bestone1;
            execution0[(T-1)*(N+1)+j] = bestk0;
            execution1[(T-1)*(N+1)+j] = bestk1;
        }
    }
    
    /* in the middle */
    for (t= T - 2; t>= 0; t--){
        for (j = 0; j <= N; j++){
            
            bestone0 = 0;
            bestone1 = 0;
            bestk0 = 0;
            bestk1 = 0;
            /** enumerate possibilities **/
            for (int k=0; k<=std::min(99,j);++k){
                candidate = shift0[k]*(k+optimal0[(t+1)*(N+1)+j-k]);
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = k;
                }
                if (candidate > bestone1){
                    bestone1 = candidate;
                    bestk1 = k;
                }
            }
            for (int k=1;k<=std::min(9,int(floor(j/100.0)));++k){
                candidate = shift0[100*k]*(100*k+optimal1[(t+1)*(N+1)+j-100*k]);
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = 100*k;
                }
                if (candidate > bestone1){
                    bestone1 = candidate;
                    bestk1 = 100*k;
                    
                }
            }
            for (int k=1; k<=int(floor(j/1000.0));++k){
                candidate = shift0[1000*k]*(1000*k+optimal1[(t+1)*(N+1)+j-1000*k]);
                if (candidate > bestone0){
                    bestone0 = candidate;
                    bestk0 = 1000*k;
                    
                }
                candidate = shift1[1000*k]*(1000*k+optimal1[(t+1)*(N+1)+j-1000*k]);
                if (candidate > bestone1){
                    bestone1 = candidate;
                    bestk1 = 1000*k;
                    
                }
            }
            
            
            optimal0[t*(N + 1) + j] = bestone0;
            optimal1[t*(N + 1) + j] = bestone1;
            execution0[t*(N + 1) + j] = bestk0;
            execution1[t*(N + 1) + j] = bestk1;
        }
        printf("done with stage t = %d\n", t);
    }
    
    printf("optimal value for trade sequencing = %4.9f\n", optimal0[N]);
    double P = 1.0;
    double sum = 0.0;
    int N_remaining = N;
    int large = 0;
    int bestk;
    for (int t = 0; t < T; ++t){
        if (large == 0){
            bestk = execution0[t*(N+1)+N_remaining];
            P *= shift0[bestk];
        }
        else{
            bestk = execution1[t*(N+1)+N_remaining];
            P *= shift1[bestk];
        }
        if (bestk >= 100){
            large = 1;
        }
        else{
            large = 0;
        }
        N_remaining -= bestk;
        best_executions[t]=bestk;
        prices[t] = P;
        printf("When t = %d, we sell %d, the price now is %g\n",t, bestk, P);
        sum += P * bestk;
    }
    printf("check by dot product: %4.9f\n",sum);
    return retcode;

}
