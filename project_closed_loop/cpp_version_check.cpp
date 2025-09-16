#include <iostream>
int main() {
    std::cout << "C++ version: " << __cplusplus << std::endl;
    #ifdef __cpp_lib_math_special_functions
    std::cout << "Special math functions available: " << __cpp_lib_math_special_functions << std::endl;
    #else
    std::cout << "Special math functions NOT available" << std::endl;
    #endif
    return 0;
}
