// to compile: make seq || make all
// to execute: ./jacobiseq <ordem_matriz> <seed> <option_debug>
//

#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <math.h>

#define MAX_ITERACOES 200

int main(int argc,char **argv){
    
    // Argumentos de entrada
    if ( argc  != 4)
	{
	    printf("Wrong arguments. Please use main <ordem_matriz> <num_threads> <seed> <option_debug>\n");
	    exit(0);
	}

    int N = atoi(argv[1]);
    float seed = *argv[2];
    int debug = atoi(argv[3]);

    // Alocacao de memoria para matriz A e vetor B
    float ** matrix = (float **) malloc(sizeof(float *) * N);
    float * vet_b = (float *) malloc(sizeof(float) * N);
    for(int i = 0; i< N; i++){
        matrix[i] = (float *) malloc(sizeof(float) * N);
    }

    // Vetor que armazena a diagonal original da matriz A para posterior substituicao na equacao
    float * vet_diag = (float *) malloc(sizeof(float) * N);

    // Define a semente para geracao de numeros aleatorios
    srand(seed);

    // Inicializa a matriz A e o vetor B com valores aleatorios
    for(int i = 0; i< N; i++){
        // Gera umal linha da matriz A
        for(int j = 0; j< N; j++){
            matrix[i][j] = rand()%100;
        }

        // Soma a linha atual da matriz A
        float soma_linha = 0;
        for(int j = 0; j< N; j++){
            soma_linha += fabs(matrix[i][j]);
        }

        // Verifica se a matriz eh diagonalmente dominante
        if(fabs(matrix[i][i]) < soma_linha - fabs(matrix[i][i])){ // Diagonal deve ser maior que a soma dos outros elementos da linha
            matrix[i][i] = soma_linha; // corrige a diagonal para ser maior que a soma dos outros elementos da linha
        }
        // Imprime a linha da matriz A --- DEBUG
        for(int j = 0; j< N; j++){
            printf("%.1f\t", matrix[i][j]);
        }
        printf("\n");

        // Gera elemento do vetor B
        vet_b[i] = rand()%100;
    }


    // Trecho para Debug do programa  --- DEBUG
    if(debug == 1){
        // Le a matriz A e o vetor B
        for(int i = 0; i< N; i++){
            for(int j = 0; j< N; j++){
                scanf("%f", &(matrix[i][j]));
            }
        }
        for(int i = 0; i< N; i++){
            scanf("%f", &(vet_b[i]));
        }
    }

    // Imprime a matriz A e o vetor B na tela --- DEBUG
    for(int i = 0; i< N; i++){
        for(int j = 0; j< N; j++){
            printf("%.1f\t", matrix[i][j]);
        }
        printf("=\t%.1f", vet_b[i]);
        printf("\n");
    }
    printf("\n");

    // Normaliza a matriz A e o vetor B e armazena a diagonal original da matriz A
    for(int i = 0; i< N; i++){
        for(int j = 0; j< N; j++){
            if(i != j){
                matrix[i][j] = matrix[i][j] / matrix[i][i]; // normaliza cada linha em relacao ao elemento da diagonal
            }
        }
        vet_b[i] = vet_b[i] / matrix[i][i];
        vet_diag[i] = matrix[i][i];
        matrix[i][i] = 0;   // zera a diagonal da matriz A
    }

    // Imprime a matriz A e o vetor B normalizados --- DEBUG
    for(int i = 0; i< N; i++){
        for(int j = 0; j< N; j++){
            printf("%f ", matrix[i][j]);
        }
        printf("= %f", vet_b[i]);
        printf("\n");
    }

    // Alocacao de memoria para o vetor solucao X, para o novo vetor X
    float * vet_x = (float *) malloc(sizeof(float) * N);
    float * vet_new_x = (float *) malloc(sizeof(float) * N);

    // Alocacao de memoria para as variaveis de calculo do erro (criterio de parada)
    float * diff = (float *) malloc(sizeof(float) * N);
    float max_new_x = fabs(vet_new_x[0]);
    float max_diff = 0;
    float error = 1; // erro inicial para entrar no loop while

    // Inicializacao dos vetores X e novo X
    for(int i = 0; i< N; i++){
        vet_x[i] = vet_b[i]; // chute inicial é o vetor B
        vet_new_x[i] = vet_x[i]; // novo vetor X sempre comeca igual o vetor B pois x[i]k+1 = B*[i] - (A*[i j].x[j]k), para i <> j e 0 >= j < n
    }

    int cont = 0; // contador de iteracoes --- DEBUG

    while(error > 0.001 && cont < MAX_ITERACOES){ // loop para realizar iteracoes ate satisfazer o criterio de parada
        // Inicio da iteracao
        printf("\nIteracao %d\n", cont);

        // Calculo do novo vetor X  -> x[i]k+1 = B*[i] - (A*[i j].x[j]k), para i <> j e 0 >= j < n
        for(int i = 0; i< N; i++){
            for(int j = 0; j< N; j++){
                vet_new_x[i] -= matrix[i][j] * vet_x[j];
            }
        }

        // Imprime o vetor X e o novo vetor X --- DEBUG
        printf("\nOld x:");
        for(int i = 0; i< N; i++){
            printf("%f ", vet_x[i]);
        }
        printf("\nNew x:");
        for(int i = 0; i< N; i++){
            printf("%f ", vet_new_x[i]);
        }
        printf("\nDiff: ");
        for(int i = 0; i< N; i++){
            printf("%f ", diff[i]);
        }
        printf("\n");

        // Calculo do erro (criterio de parada)
        max_diff = 0;
        max_new_x = fabs(vet_new_x[0]);
        for(int i = 0; i< N; i++){
            diff[i] = fabs(vet_new_x[i] - vet_x[i]); // calcula diferenca entre o novo vetor X e o vetor X
            if(diff[i] > max_diff){
                max_diff = diff[i]; // calcula o maior valor da diferenca
            }
            if(fabs(vet_new_x[i]) > max_new_x){
                max_new_x = fabs(vet_new_x[i]); // calcula o maior valor do novo vetor X
            }
        }
        printf(" Max new X: %f", max_new_x); // --- DEBUG

        error = max_diff / max_new_x;
        printf("\nError: %f\n", error);

        if(error < 0.001){ // sai do loop while quando satisfaz o erro minimo
            break;
        }

        // Atualiza o vetor X para a proxima iteracao
        for(int i = 0; i< N; i++){
            vet_x[i] = vet_new_x[i]; // vetor X recebe o novo vetor X (proximo chute)
            vet_new_x[i] = vet_b[i]; // vetor novo X sempre comeca com B
        }

        cont++; // contagem de iteracoes --- DEBUG
    }

    // Imprime o vetor solucao X
    printf("\nVetor solucao: ");
    for(int i = 0; i< N; i++){
        printf("%.3f ", vet_x[i]);
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
                matrix[linha][i] *= vet_diag[linha];
            }else{
                matrix[linha][i] = vet_diag[linha];
            }
            printf("%.2fx_%d +\t", matrix[linha][i], i);
            // Avalia equacao com o valor do vetor X
            result += matrix[linha][i] * vet_x[i];
        }
        printf("= %.2f\n", vet_b[linha]*vet_diag[linha]);
        printf("Resultado da atribuicao na linha %d: %.2f\n", linha, result);
    }
    
    // Desalocacao de Memoria
    for(int i = 0; i< N; i++){
        free(matrix[i]);
    }
    free(matrix);
    free(vet_b);
    free(vet_diag);
    free(vet_x);
    free(vet_new_x);
    free(diff);
    
    return 0;
}
