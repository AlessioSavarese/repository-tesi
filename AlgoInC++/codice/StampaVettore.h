#include <iostream>
#include <vector>
/**
 * Classe di utilità per stampare un vettore.
 * @param array: Vettore da stampare.
 */
class StampaVettore
{
private:
public:
    void printVector(const std::vector<uint64_t> &vettore);

    void printVector(const std::vector<double> &vettore);
};