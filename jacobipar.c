// to compile: make seq || make all
// to execute: ./jacobiseq <ordem_matriz> <seed> <option_debug>
//

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define MAX_ITERACOES 300
#define MAX_MATRIX_VALUE 100

//#pragma omp declare simd linear(matrix:1) linear(vet_x:1) linear(vet_new_x:1)
void calculate_new_x (float * matrix, float * vet_x, float * vet_new_x, float * vet_b, int N, int T)
{
    // Atualiza o vetor X para a proxima iteracao
    #pragma omp parallel for simd num_threads(T) shared(vet_new_x) firstprivate(matrix, vet_x, vet_b, N)
    for(int i = 0; i< N; i++){
        vet_x[i] = vet_new_x[i]; // vetor X recebe o novo vetor X (proximo chute)
        vet_new_x[i] = vet_b[i]; // vetor novo X sempre comeca com B
        for(int j = 0; j< N; j++){
            vet_new_x[i] -= matrix[i*N + j] * vet_x[j];
        }
    }
}

void calculate_error(float * vet_x, float * vet_new_x, float * diff, int N, int debug, float * error, int T)
{
    float max_diff = 0;
    float max_new_x = fabs(vet_new_x[0]);
    #pragma omp parallel for simd num_threads(T) firstprivate(diff, vet_x, vet_new_x, N) reduction(max:max_diff, max_new_x)
    for(int i = 0; i< N; i++){
        diff[i] = fabs(vet_new_x[i] - vet_x[i]); // calcula diferenca entre o novo vetor X e o vetor X
       if(diff[i] > max_diff){
            max_diff = diff[i]; // calcula o maior valor da diferenca
        }
        if(fabs(vet_new_x[i]) > max_new_x){
            max_new_x = fabs(vet_new_x[i]); // calcula o maior valor do novo vetor X
        }
    }
    *error = max_diff / max_new_x;

    if(debug == 3){
        //printf(" Max new X: %f %f", max_diff, max_new_x); // --- DEBUG
        printf("\nError: %f\n", *error);            
    }
}

int main(int argc,char **argv){

    double wtime;
    
    // Argumentos de entrada
    if ( argc  != 5)
	{
	    printf("Wrong arguments. Please use main <ordem_matriz> <seed> <num_threads> <option_debug>\n");
	    exit(0);
	}

    int N = atoi(argv[1]);
    float seed = *argv[2];
    int T = atoi(argv[3]);
    int debug = atoi(argv[4]);

    // Alocacao de memoria para matriz A e vetor B
    float * matrix = (float *) malloc(sizeof(float *) * N * N); // matriz linearizada
    float * vet_b = (float *) malloc(sizeof(float) * N);

    // Vetor que armazena a diagonal original da matriz A para posterior substituicao na equacao
    float * vet_diag = (float *) malloc(sizeof(float) * N);

    // Define a semente para geracao de numeros aleatorios
    srand(seed);

    // Inicializa a matriz A e o vetor B com valores aleatorios
    for(int i = 0; i< N; i++){
        // Gera umal linha da matriz A
        for(int j = 0; j< N; j++){
            matrix[i*N + j] = rand()%MAX_MATRIX_VALUE;
        }

        // Soma a linha atual da matriz A
        float soma_linha = 0;
        for(int j = 0; j< N; j++){
            soma_linha += fabs(matrix[i*N + j]);
        }

        // Verifica se a matriz eh diagonalmente dominante
        if(fabs(matrix[i*N + i]) < soma_linha - fabs(matrix[i*N + i])){ // Diagonal deve ser maior que a soma dos outros elementos da linha
            //printf("Corrigindo diagonal\n");
            matrix[i*N +i] = soma_linha; // corrige a diagonal para ser maior que a soma dos outros elementos da linha
        }
        // Imprime a linha da matriz A --- DEBUG
        if(debug == 1){
            for(int j = 0; j< N; j++){
                printf("%.1f\t", matrix[i*N + j]);
            }
            printf("\n");
        }

        // Gera elemento do vetor B
        vet_b[i] = rand()%100;
    }


    // Trecho para Debug do programa  --- DEBUG
    if(debug == 2){
        // Imprime a matriz A e o vetor B na tela --- DEBUG
        for(int i = 0; i< N; i++){
            for(int j = 0; j< N; j++){
                printf("%.1f\t", matrix[i*N + j]);
            }
            printf("=\t%.1f", vet_b[i]);
            printf("\n");
        }
        printf("\n");
    }
    
    wtime = omp_get_wtime();

    // Normaliza a matriz A e o vetor B e armazena a diagonal original da matriz A
    #pragma omp parallel num_threads(T) shared(matrix, vet_b, vet_diag, N)
    {
        #pragma omp single
        {
            for(int i = 0; i< N; i++){
                #pragma omp task depend(in: matrix[i*N + i])
                {
                    vet_b[i] = vet_b[i] / matrix[i*N +i];
                    vet_diag[i] = matrix[i*N + i];
                } 
                #pragma omp task depend(out: matrix[i*N + i])
                for(int j = 0; j< i; j++){
                    matrix[i*N + j] = matrix[i*N + j] / matrix[i*N + i]; // normaliza cada linha em relacao ao elemento da diagonal
                }
                #pragma omp task depend(out: matrix[i*N + i])
                for(int j = i+1; j< N; j++){
                    matrix[i*N + j] = matrix[i*N + j] / matrix[i*N + i]; // normaliza cada linha em relacao ao elemento da diagonal
                }
                #pragma omp taskwait
                {
                    matrix[i*N + i] = 0;   // zera a diagonal da matriz A
                }
            }
        }
    }
    

    // Imprime a matriz A e o vetor B normalizados --- DEBUG
    if(debug == 1){
        for(int i = 0; i< N; i++){
            for(int j = 0; j< N; j++){
                printf("%f ", matrix[i*N + j]);
            }
            printf("= %f", vet_b[i]);
            printf("\n");
        }
    }

    // Alocacao de memoria para o vetor solucao X, para o novo vetor X
    float * vet_x = (float *) malloc(sizeof(float) * N);
    float * vet_new_x = (float *) malloc(sizeof(float) * N);

    float error = 1; // erro inicial
    float * diff = (float *) malloc(sizeof(float) * N); // vetor para armazenar a diferenca entre o vetor X e o novo vetor X
    int flag_cancel = 0; // flag para cancelar o loop 

    // Inicializacao dos vetores X e novo X
    //#pragma omp parallel for num_threads(T) shared(vet_x, vet_b, vet_new_x, N)
    for(int i = 0; i< N; i++){
        vet_x[i] = vet_b[i]; // chute inicial é o vetor B
        vet_new_x[i] = vet_x[i]; // novo vetor X sempre comeca igual o vetor B pois x[i]k+1 = B*[i] - (A*[i j].x[j]k), para i <> j e 0 >= j < n
    }

    int cont = 0; // contador de iteracoes --- DEBUG
    //#pragma omp parallel for num_threads(T) shared(matrix, vet_b, vet_x, vet_new_x, flag_cancel) firstprivate(error, N, debug)
    //for(cont = 0; cont < MAX_ITERACOES; cont++){ // loop para realizar iteracoes ate satisfazer o criterio de parada
    while(error > 0.001 /*&& cont < MAX_ITERACOES*/){
        
        // Calculo do novo vetor X  -> x[i]k+1 = B*[i] - (A*[i j].x[j]k), para i <> j e 0 >= j < n
        calculate_new_x(matrix, vet_x, vet_new_x, vet_b, N, T);

         if(debug  == 4){
            printf("\nNew x:");
            for(int i = 0; i< N; i++){
                printf("%f ", vet_new_x[i]);
            }
        }

        // Imprime o vetor X e o novo vetor X --- DEBUG
        if(debug == 1){
            printf("\nOld x:");
            for(int i = 0; i< N; i++){
                printf("%f ", vet_x[i]);
            }
            printf("\nNew x:");
            for(int i = 0; i< N; i++){
                printf("%f ", vet_new_x[i]);
            }
            printf("\n");
        }

        // Calculo do erro (criterio de parada)
        calculate_error(vet_x, vet_new_x, diff, N, debug, &error, T);
        
        if(debug == 1){
            //printf(" Max new X: %f", max_new_x); // --- DEBUG
            printf("\nError: %f\n", error);            
        }

        if(error < 0.001){ // sai do loop while quando satisfaz o erro minimo
            break;
        }

    
        cont++;
    }

    // Fim do calculo do tempo de execucao
    wtime = omp_get_wtime() - wtime;
    printf("Elapsed wall clock time = %f  \n", wtime);

    // Imprime o vetor solucao X
    if(debug == 1){
        //printf("\nVetor solucao: ");
        for(int i = 0; i< N; i++){
            //printf("%.3f ", vet_x[i]);
        }
    }

    printf("\nDigite o indice da equacao que deseja substituir: ");
    int linha;
    scanf("%d", &linha);
    float result = 0;
    if (linha >=0 && linha < N)
    {
        for(int i = 0; i< N; i++){
            // Reconstroi a linha original da matriz A (sem normalizacao)
            if(i != linha){
                matrix[linha*N + i] *= vet_diag[linha];
            }else{
                matrix[linha*N + i] = vet_diag[linha];
            }
            if(debug == 1)
                printf("%.2fx_%d +\t", matrix[linha*N + i], i);
            // Avalia equacao com o valor do vetor X
            result += matrix[linha*N +i] * vet_x[i];
        }
        if(debug == 1)
            printf("= %.2f - Error: %f\n", vet_b[linha]*vet_diag[linha], error);
        printf("Resultado da atribuicao na linha %d (%d iteracoes): %.6f\n", linha, cont, result);
        printf("Erro: %.6f\n", error);
    }
    

    // Desalocacao de Memoria
    free(matrix);
    free(vet_b);
    free(vet_diag);
    free(vet_x);
    free(vet_new_x);
    free(diff);
    
    return 0;
}
