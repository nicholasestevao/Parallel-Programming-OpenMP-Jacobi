#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define TIME_OF_EXECUTION 30
#define ORDEM_MATRIZ 1000


int main() {
    char comando[100];
    int retorno = 0;


    // Teste Sequencial
    for (int i = 0; i < TIME_OF_EXECUTION; i++) {
        sprintf(comando, "time -p -o testeSeq.txt -a ./jacobiseq.out %d %ld", ORDEM_MATRIZ, clock());
        retorno = system(comando);

        if (retorno == -1) {
            printf("Ocorreu um erro ao tentar executar o comando.\n");
        }
    }

    // Teste Paralelo
    int max_threads = omp_get_max_threads(); // roda com o número máximmo de theads lógicas disponíveis
    for(int i = 0; i < TIME_OF_EXECUTION; i++) {
        sprintf(comando, "time -p -o testePar.txt -a ./jacobipar.out %d %ld %d", ORDEM_MATRIZ, clock(), max_threads);

        retorno = system(comando);

        if (retorno == -1) {
            printf("Ocorreu um erro ao tentar executar o comando.\n");
        }
    }
    
    return 0;
}
