#ifndef HASHFUNCTION_H // Controlla se la macro non è già definita
#define HASHFUNCTION_H // Definisci la macro per evitare inclusioni successive

#include "Prime.h"
#include <cstddef>
#include <iostream>

class HashFunction
{
private:
    size_t a; ///< Coefficiente moltiplicativo
    size_t b; ///< Coefficiente additivo
    size_t m; ///< Modulo Universo / Dimensione Firma del FSS
    size_t p; ///< Numero primo maggiore di m
    
    size_t t;
    size_t p_t;
public:
    /**
     * Costruttore.
     * @param m: Modulo per la funzione hash.
     * @param seed: Semina per il generatore casuale.
     */
    HashFunction(size_t m, size_t seed, size_t t = 0);

    /**
     * Mappa un elemento `x` utilizzando la funzione hash.
     * @param x: Valore da mappare.
     * @return Valore hashato.
     */
    size_t map(size_t x);

    /**
     * Mappa un elemento `x` utilizzando la funzione hash.
     * @param x: Valore da mappare.
     * @return Valore hashato.
     */
    std::pair<size_t, double> map(size_t x, size_t i);
};

#endif // Fine della protezione contro l'inclusione multipla