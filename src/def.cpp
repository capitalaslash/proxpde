#include "def.hpp"

namespace proxpde
{

#ifdef NDEBUG
std::ofstream Utils::debug = std::ofstream{"dev/null"};
#endif

} // namespace proxpde
