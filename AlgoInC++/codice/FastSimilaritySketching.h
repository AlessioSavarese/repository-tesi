#include <vector>
#include "HashFunction.h"
#include <iostream>

class FastSimilaritySketching
{
private:
    size_t t;                   ///< Numero di bin.
    // HashFunction *funzioniHash; ///< Funzioni hash.
    std::vector<HashFunction *> vettoreFirme; ///< Funzioni hash.

public:
    /**
     * Costruttore.
     * @param t: Numero di bin.
     * @param seed: Semina per il generatore casuale.
     */
    FastSimilaritySketching(size_t t, size_t m, size_t seed);

    /**
     * Calcola la firma di un vettore.
     * @param s: Vettore di input.
     * @return Firma come vettore di reali (double).
     */
    std::vector<double> computeSignature(std::vector<uint64_t> s);

    bool separationProcedure(size_t r, double gamma, std::vector<double> S_A, std::vector<double> S_B);
};
