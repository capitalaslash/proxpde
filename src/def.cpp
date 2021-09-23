#include "def.hpp"

#ifdef NDEBUG
std::ofstream Utils::debug = std::ofstream{"dev/null"};
#endif

std::ofstream Utils::filelog = std::ofstream{"minifem.log"};

