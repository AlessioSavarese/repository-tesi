#include <sstream>
class LettoreFile
{
public:
    static void read(std::string &filename, std::vector<std::pair<std::vector<uint64_t>, std::vector<uint64_t>>> &coppie);
};