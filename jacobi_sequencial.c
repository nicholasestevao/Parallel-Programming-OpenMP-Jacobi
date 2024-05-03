// to compile: make seq || make all
// to execute: ./jacobi_sequencial <ordem_matriz> <seed>

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>
#include <time.h>

#define MAX_ITERACOES 50000
#define MAX_MATRIX_VALUE 100
#define PRECISAO_JACOBI 0.0001

// Inicializa a matriz A e o vetor B com valores aleatorios
void init_matrix(double *matrix, double *vet_b, double *vet_diag, int N) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j< N; j++) {
            matrix[i*N + j] = rand() % MAX_MATRIX_VALUE;
        }

        double soma_linha = 0;
        for(int j = 0; j < N; j++) {
            soma_linha += fabs(matrix[i*N + j]);
        }

        if(fabs(matrix[i*N + i]) < soma_linha - fabs(matrix[i*N + i])) {
            matrix[i*N +i] = soma_linha + 1;
        }

        vet_b[i] = rand()%100;
    }
}

// Normaliza a matriz A e o vetor B e armazena a diagonal original da matriz A
void normalize_matrix(double *matrix, double *vet_b, double *vet_diag, int N) {
    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            if(i != j){
                matrix[i*N + j] = matrix[i*N + j] / matrix[i*N + i];
            }
        }
        vet_b[i] = vet_b[i] / matrix[i*N +i];
        vet_diag[i] = matrix[i*N + i];
        matrix[i*N + i] = 0;
    }
}

// Calculo do novo vetor X
void calculate_new_x(double *matrix, double *vet_b, double *vet_x, double *vet_new_x, int N) {
    for(int i = 0; i < N; i++) {
        vet_x[i] = vet_new_x[i];
        vet_new_x[i] = vet_b[i];
    }

    for(int i = 0; i < N; i++) {
        for(int j = 0; j < N; j++){
            vet_new_x[i] -= matrix[i*N + j] * vet_x[j];
        }
    }
}

// Calculo do erro (criterio de parada)
void calculate_error(double *vet_x, double *vet_new_x, double *diff, double *error, int N) {
    double max_diff = 0;
    double max_new_x = fabs(vet_new_x[0]);
    
    for(int i = 0; i < N; i++) {
        diff[i] = fabs(vet_new_x[i] - vet_x[i]);
        if(diff[i] > max_diff) {
            max_diff = diff[i];
        }

        if(fabs(vet_new_x[i]) > max_new_x) {
            max_new_x = fabs(vet_new_x[i]);
        }
    }

    *error = max_diff / max_new_x;
}

int main(int argc,char **argv) {
    if ( argc  != 3) {
	    printf("Wrong arguments. Please use main <ordem_matriz> <seed> <option_debug>\n");
	    exit(0);
	}

    int N = atoi(argv[1]);
    int seed = atoi(argv[2]);

    double * matrix = (double *) malloc(sizeof(double *) * N * N);
    double * vet_b = (double *) malloc(sizeof(double) * N);
    double * vet_diag = (double *) malloc(sizeof(double) * N);
    if (matrix == NULL || vet_b == NULL || vet_diag == NULL) {
        printf("Erro de alocação de memória\n");
        exit(1);
    }
    
    double wtime;
    double cpu_time_used;
    clock_t start, end;

    // Record the starting time
    start = clock();

    srand(seed);

    init_matrix(matrix, vet_b, vet_diag, N);
    
    wtime = omp_get_wtime();

    normalize_matrix(matrix, vet_b, vet_diag, N);

    double * vet_x = (double *) malloc(sizeof(double) * N);
    double * vet_new_x = (double *) malloc(sizeof(double) * N);
    double * diff = (double *) malloc(sizeof(double) * N);
    if (vet_x == NULL || vet_new_x == NULL || diff == NULL) {
        printf("Erro de alocação de memória\n");
        exit(1);
    }

    double error = 1;

    // Inicializacao dos vetores X e novo X
    for(int i = 0; i < N; i++) {
        vet_x[i] = vet_b[i]; 
        vet_new_x[i] = vet_x[i];
    }

    int cont = 0;


    while(error > PRECISAO_JACOBI && cont < MAX_ITERACOES) {
        calculate_new_x(matrix, vet_b, vet_x, vet_new_x, N);
        calculate_error(vet_x, vet_new_x, diff, &error, N);

        cont++;
    }

    wtime = omp_get_wtime() - wtime;
    printf("User time = %f ms \n", wtime*1000);

    printf("\nDigite o indice da equacao que deseja substituir: ");
    int linha;
    scanf("%d", &linha);
    double result = 0;
    if (linha >= 0 && linha < N) {
        for(int i = 0; i < N; i++) {
            if(i != linha){
                matrix[linha*N + i] *= vet_diag[linha];
            } else {
                matrix[linha*N + i] = vet_diag[linha];
            }
            result += matrix[linha*N +i] * vet_x[i];
        }

        printf("Resultado da atribuicao na linha %d (%d iteracoes): %.6f\n", linha, cont, result);
        printf("Valor esperado: %f\n", vet_b[linha]*vet_diag[linha]);
        printf("Erro: %.6f\n", error);
    }

    // Record the ending time
    end = clock();

    // Calculate the CPU time used
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // Output the CPU time and real time
    printf("CPU time: %f seconds\n", cpu_time_used);

    free(matrix);
    free(vet_b);
    free(vet_diag);
    free(vet_x);
    free(vet_new_x);
    free(diff);
    
    return 0;
}
