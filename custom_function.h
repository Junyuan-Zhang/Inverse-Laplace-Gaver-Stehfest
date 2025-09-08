#ifndef CUSTOM_FUNCTION_H
#define CUSTOM_FUNCTION_H

#include <gmpxx.h>

// User-defined custom functions
extern const bool HAS_ANALYTICAL_SOLUTION;
extern mpf_class custom_target_function(const mpf_class& s);
extern mpf_class custom_analytical_function(const mpf_class& t);

#endif // CUSTOM_FUNCTION_H
