#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <time.h>

#define TIME_OF_EXECUTION 30
#define ORDEM_MATRIZ 1500
#define LINE_FOR_VERIFICATION 0

void txtToCSV(char *nome_txt, char *nome_csv);

int main()
{
    char *comando = (char*) malloc(100 * sizeof(char));
    char *nome_txt = (char*) malloc(50 * sizeof(char));
    char *nome_csv = (char*) malloc(50 * sizeof(char));
    int retorno = 0;

    sprintf(nome_txt, "testeSeq_%d.txt", ORDEM_MATRIZ);
    sprintf(nome_csv, "testeSeq_%d.csv", ORDEM_MATRIZ);
    // Teste Sequencial
    for (int i = 0; i < TIME_OF_EXECUTION; i++)
    {
        sprintf(comando, "time -p -o %s -a ./jacobiseq.out %d %d %d", nome_txt, ORDEM_MATRIZ, (int)clock(), LINE_FOR_VERIFICATION);
        retorno = system(comando);

        if (retorno == -1)
        {
            printf("Ocorreu um erro ao tentar executar o comando.\n");
        }
    }
    // transforma em csv
    txtToCSV(nome_txt, nome_csv);

    // Teste Paralelo
    int max_threads = omp_get_max_threads(); // roda com o número máximmo de theads lógicas disponíveis
    for (int j = 2; j <= max_threads; j+=2) {
        sprintf(nome_txt, "testePar_%d_threads%d.txt", ORDEM_MATRIZ, j);
        sprintf(nome_csv, "testePar_%d_threads%d.csv", ORDEM_MATRIZ, j);
        for (int i = 0; i < TIME_OF_EXECUTION; i++) {
            sprintf(comando, "time -p -o %s -a ./jacobipar.out %d %d %d %d", nome_txt, ORDEM_MATRIZ, (int)clock(), j, LINE_FOR_VERIFICATION);

            retorno = system(comando);

            if (retorno == -1)
            {
                printf("Ocorreu um erro ao tentar executar o comando.\n");
            }
        }
        // transforma em csv
        txtToCSV(nome_txt, nome_csv);
    }
    return 0;
}

void txtToCSV(char *nome_txt, char *nome_csv)
{
    FILE *arq_txt = fopen(nome_txt, "r");
    if (arq_txt == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nome_txt);
        return;
    }

    FILE *arq_csv = fopen(nome_csv, "w");
    if (arq_csv == NULL)
    {
        printf("Erro ao abrir o arquivo %s\n", nome_csv);
        return;
    }
    fprintf(arq_csv, "Real,User,Sys\n");

    char *linha = (char *)malloc(100 * sizeof(char));
    float *real = (float *)malloc(TIME_OF_EXECUTION * sizeof(float));
    float *user = (float *)malloc(TIME_OF_EXECUTION * sizeof(float));
    float *sys = (float *)malloc(TIME_OF_EXECUTION * sizeof(float));

    int i = 0;
    int j = 0;
    while (fgets(linha, 100, arq_txt) != NULL)
    {
        if (linha[0] == 'r')
        {
            sscanf(linha, "real %f", &real[i]);
            j++;
        }
        else if (linha[0] == 'u')
        {
            sscanf(linha, "user %f", &user[i]);
            j++;
        }
        else if (linha[0] == 's')
        {
            sscanf(linha, "sys %f", &sys[i]);
            j++;
        }

        if (j == 3)
        {
            fprintf(arq_csv, "%f,%f,%f\n", real[i], user[i], sys[i]);
            i++;
            j = 0;
        }
    }

    fclose(arq_txt);
    fclose(arq_csv);
    free(linha);
    free(real);
    free(user);
    free(sys);

    char comando[100];
    sprintf(comando, "rm %s", nome_txt);
    if (system(comando) == -1)
    {
        printf("Erro ao tentar remover %s.\n", nome_txt);
    }

    return;
}