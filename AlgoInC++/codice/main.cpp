#include <iostream>
#include "Test.h"

int main(int argc, char *argv[])
{
    // Controllo se è stato fornito almeno un parametro
    if (argc != 2) {
        printf("\nErrore: manca parametro x");
        printf("\nUso:");
        printf("\n      0 per eseguire test sul tempo con k fissato");
        printf("\n      1 per eseguire test sul tempo con n fissato");
        printf("\n      2 per eseguire tutti i test sul tempo");
        printf("\n      3 per eseguire test sulla qualità");
        printf("\n      4 per eseguire tutti i test\n\n");
        return 1;  // esco con errore
    }
    // Converto il parametro in intero
    int x = atoi(argv[1]);
    
    // Imposta il modulo dell'universo come il primo numero primo maggiore o uguale a 2^30 (1073741824)
    size_t m = 1073741824;

    //imposta ripetizioni esperimenti
    int rep = 100;




    
    // ------------ TEST TEMPO k fissato ------------
    auto start = std::chrono::high_resolution_clock::now(); // CRONOMETRO
    std::vector<int> vettore_k_fissato = { 16, 256, 4096 }; // k fissato a salti di 2^4
    if (x == 0 || x == 2 || x == 4)
    {
        std::vector<int> n_range;
        for (int i = (1 << 10); i <= (1 << 20); i *= 2) {
            n_range.push_back(i); // n range || Ciclo per generare le potenze da 2^10 fino a 2^20
        }
        for(int k_fissato : vettore_k_fissato){
            Test::test_time_vs_n(k_fissato, n_range, rep, m);

            // Esegui lo script Python per generare il grafico dei tempi in funzione di n
            system( ("python3 graficoTempo.py k " + std::to_string(k_fissato)).c_str() );

            // system( ("open grafico_tempo_k=" + std::to_string(k_fissato) + ".png").c_str() );
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> duration = end - start;
    std::cout << "\nTempo impiegato [TEST TEMPO k fissato]: " << duration.count() << " secondi\n\n" << std::endl;




    
    // ------------ TEST TEMPO n fissato ------------
    start = std::chrono::high_resolution_clock::now(); // CRONOMETRO
    std::vector<int> vettore_n_fissato = { 1000,10000 , 100000, 1000000 }; // n fissato potenze di 10 da 10^2
    if (x == 1 || x == 2 || x == 4)
    {
        std::vector<int> k_range;
        for (int i = 1; i <= 8192; i *= 2) k_range.push_back(i); // k range || Ciclo per generare le potenze di 2 fino a 2^16

        for(int n_fissato : vettore_n_fissato){
            Test::test_time_vs_k(k_range, n_fissato, rep, m);

            // // Esegui lo script Python per generare il grafico dei tempi in funzione di k
            system( ("python3 graficoTempo.py n " + std::to_string(n_fissato)).c_str() );

            // system( ("open grafico_tempo_n=" + std::to_string(n_fissato) + ".png").c_str() );
        }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "\nTempo impiegato [TEST TEMPO n fissato]: " << duration.count() << " secondi\n\n" << std::endl;



    
    
    // ------------ TEST per analizzare la qualità ------------
    start = std::chrono::high_resolution_clock::now(); // CRONOMETRO
    if (x == 3 || x == 4)
    {
        int n = 50000;
        for(int k_fissato : vettore_k_fissato){
                Test::test_quality(k_fissato, n, rep, m);

                // Esegui lo script Python per generare il grafico dei tempi in funzione di n
                system(("python3 graficoQualita.py " + std::to_string(k_fissato)).c_str());
            }
    }
    end = std::chrono::high_resolution_clock::now();
    duration = end - start;
    std::cout << "\nTempo impiegato [TEST QUALITA']: " << duration.count() << " secondi\n\n" << std::endl;
}