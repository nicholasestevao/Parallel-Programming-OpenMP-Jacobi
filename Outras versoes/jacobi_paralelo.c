// to compile: make seq || make all
// to execute: ./jacobi_paralelo <ordem_matriz> <seed>

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
        // Gera uma linha da matriz A
        for(int j = 0; j < N; j++) {
            matrix[i*N + j] = rand() % MAX_MATRIX_VALUE;
        }

        // Soma a linha atual da matriz A
        double soma_linha = 0;
        for(int j = 0; j < N; j++) {
            soma_linha += fabs(matrix[i*N + j]);
        }
        // Verifica se a matriz eh diagonalmente dominante
        if(fabs(matrix[i*N + i]) < soma_linha - fabs(matrix[i*N + i])) { // Diagonal deve ser maior que a soma dos outros elementos da linha
            matrix[i*N +i] = soma_linha + 1; // corrige a diagonal para ser maior que a soma dos outros elementos da linha
        }

        // Gera elemento do vetor B
        vet_b[i] = rand()%100;
    }
}

// Normaliza a matriz A e o vetor B e armazena a diagonal original da matriz A
void normalize_matrix(double *matrix, double *vet_b, double *vet_diag, int N, int T) {
    int i, j = 0;

    #pragma omp parallel num_threads(T) shared(matrix, vet_b, vet_diag, N) private (i, j)
    {
        #pragma omp single
        {
            for(i = 0; i < N; i++){
                #pragma omp task depend(in: matrix[i*N + i])
                {
                    vet_b[i] = vet_b[i] / matrix[i*N +i];
                    vet_diag[i] = matrix[i*N + i];
                } 
                #pragma omp task depend(out: matrix[i*N + i])
                for(j = 0; j < i; j++){
                    matrix[i*N + j] = matrix[i*N + j] / matrix[i*N + i]; // normaliza cada linha em relacao ao elemento da diagonal
                }
                #pragma omp task depend(out: matrix[i*N + i])
                for(j = i+1; j < N; j++){
                    matrix[i*N + j] = matrix[i*N + j] / matrix[i*N + i]; // normaliza cada linha em relacao ao elemento da diagonal
                }
                #pragma omp taskwait
                {
                    matrix[i*N + i] = 0;   // zera a diagonal da matriz A
                }
            }
        }
    }
}

// Calculo do novo vetor X
void calculate_new_x(double *matrix, double *vet_b, double *vet_x, double *vet_new_x, int N, int T) {
    int i, j = 0;
    // Atualiza o vetor X para a proxima iteracao
    #pragma omp parallel for simd num_threads(T) shared(vet_new_x, matrix, vet_x, vet_b, N) private(i)
    for(i = 0; i < N; i++){
        vet_x[i] = vet_new_x[i]; // vetor X recebe o novo vetor X (proximo chute)
        vet_new_x[i] = vet_b[i]; // vetor novo X sempre comeca com B
    }

    #pragma omp parallel for collapse(2) num_threads(T) private(i, j) reduction(+:vet_new_x[:N])
    for(i = 0; i < N; i++){
        for(j = 0; j< N; j++){
            vet_new_x[i] -= matrix[i*N + j] * vet_x[j];
        }
    }
}

// Calculo do erro (criterio de parada)
void calculate_error(double *vet_x, double *vet_new_x, double *diff, double *error, int N, int T) {
    double max_diff = 0;
    double max_new_x = fabs(vet_new_x[0]);
    int i = 0;
    
    #pragma omp parallel for simd num_threads(T) firstprivate(diff, vet_x, vet_new_x, N) private(i) reduction(max:max_diff, max_new_x)
    for(i = 0; i < N; i++){
        diff[i] = fabs(vet_new_x[i] - vet_x[i]); // calcula diferenca entre o novo vetor X e o vetor X
        if(diff[i] > max_diff){
            max_diff = diff[i]; // calcula o maior valor da diferenca
        }

        if(fabs(vet_new_x[i]) > max_new_x){
            max_new_x = fabs(vet_new_x[i]); // calcula o maior valor do novo vetor X
        }
    }

    *error = max_diff / max_new_x;
}

int main(int argc,char **argv) {
    // Argumentos de entrada
    if ( argc  != 4) {
	    printf("Wrong arguments. Please use main <ordem_matriz> <seed> <option_debug>\n");
	    exit(0);
	}

    int N = atoi(argv[1]);
    int seed = atoi(argv[2]);
    int T = atoi(argv[3]);

    // Alocacao de memoria para matriz A e vetor B
    double * matrix = (double *) malloc(sizeof(double *) * N * N);
    double * vet_b = (double *) malloc(sizeof(double) * N);
    double * vet_diag = (double *) malloc(sizeof(double) * N); // Vetor que armazena a diagonal original da matriz A para posterior substituicao na equacao
    if (matrix == NULL || vet_b == NULL || vet_diag == NULL) {
        printf("Erro de alocação de memória\n");
        exit(1);
    }
    
    double wtime;
    double cpu_time_used;
    clock_t start, end;

    // Grava o tempo inicial
    start = clock();

    // Define a semente para geracao de numeros aleatorios
    srand(seed);

    init_matrix(matrix, vet_b, vet_diag, N);

    wtime = omp_get_wtime();

    normalize_matrix(matrix, vet_b, vet_diag, N, T);

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
        // Calculo do novo vetor X  -> x[i]k+1 = B*[i] - (A*[i j].x[j]k), para i <> j e 0 >= j < n
        calculate_new_x(matrix, vet_b, vet_x, vet_new_x, N, T);
        // Calculo do erro (criterio de parada)
        calculate_error(vet_x, vet_new_x, diff, &error, N, T);

        cont++;
    }

    // Fim do calculo do tempo de execucao
    wtime = omp_get_wtime() - wtime;
    printf("User time = %f ms \n", wtime*1000);

    printf("\nDigite o indice da equacao que deseja substituir: ");
    int linha;
    scanf("%d", &linha);
    double result = 0;
    if (linha >= 0 && linha < N) {
        for(int i = 0; i < N; i++) {
            // Reconstroi a linha original da matriz A (sem normalizacao)
            if(i != linha) {
                matrix[linha*N + i] *= vet_diag[linha];
            } else {
                matrix[linha*N + i] = vet_diag[linha];
            }
            // Avalia equacao com o valor do vetor X
            result += matrix[linha*N +i] * vet_x[i];
        }

        printf("Resultado da atribuicao na linha %d (%d iteracoes): %.6f\n", linha, cont, result);
        printf("Valor esperado: %f\n", vet_b[linha]*vet_diag[linha]);
        printf("Erro: %.6f\n", error);
    }

    // Grava o tempo final
    end = clock();

    // Calcula o tempo de CPU utilizado
    cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
    // Saída o tempo de CPU e tempo real
    printf("CPU time: %f seconds\n", cpu_time_used);

    free(matrix);
    free(vet_b);
    free(vet_diag);
    free(vet_x);
    free(vet_new_x);
    free(diff);
    
    return 0;
}
