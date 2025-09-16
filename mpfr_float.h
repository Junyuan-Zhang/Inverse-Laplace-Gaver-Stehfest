#ifndef MPFR_FLOAT_H
#define MPFR_FLOAT_H

#include <mpfr.h>
#include <string>

// MPFR Wrapper for easier integration
class MPFRFloat {
private:
    mpfr_t value;
    static int default_precision;
    
public:
    MPFRFloat(double d = 0.0) {
        mpfr_init2(value, default_precision);
        mpfr_set_d(value, d, MPFR_RNDN);
    }
    
    MPFRFloat(const char* str) {
        mpfr_init2(value, default_precision);
        mpfr_set_str(value, str, 10, MPFR_RNDN);
    }
    
    MPFRFloat(const MPFRFloat& other) {
        mpfr_init2(value, mpfr_get_prec(other.value));
        mpfr_set(value, other.value, MPFR_RNDN);
    }
    
    MPFRFloat& operator=(const MPFRFloat& other) {
        if (this != &other) {
            mpfr_set_prec(value, mpfr_get_prec(other.value));
            mpfr_set(value, other.value, MPFR_RNDN);
        }
        return *this;
    }
    
    ~MPFRFloat() {
        mpfr_clear(value);
    }
    
    // Arithmetic operations
    MPFRFloat operator+(const MPFRFloat& other) const {
        MPFRFloat result;
        mpfr_add(result.value, value, other.value, MPFR_RNDN);
        return result;
    }
    
    MPFRFloat operator-(const MPFRFloat& other) const {
        MPFRFloat result;
        mpfr_sub(result.value, value, other.value, MPFR_RNDN);
        return result;
    }
    
    MPFRFloat operator*(const MPFRFloat& other) const {
        MPFRFloat result;
        mpfr_mul(result.value, value, other.value, MPFR_RNDN);
        return result;
    }
    
    MPFRFloat operator/(const MPFRFloat& other) const {
        MPFRFloat result;
        mpfr_div(result.value, value, other.value, MPFR_RNDN);
        return result;
    }
    
    MPFRFloat& operator+=(const MPFRFloat& other) {
        mpfr_add(value, value, other.value, MPFR_RNDN);
        return *this;
    }
    
    MPFRFloat& operator-=(const MPFRFloat& other) {
        mpfr_sub(value, value, other.value, MPFR_RNDN);
        return *this;
    }
    
    // Conversion functions
    double to_double() const {
        return mpfr_get_d(value, MPFR_RNDN);
    }
    
    std::string to_string(int digits = 30) const {
        char* str = nullptr;
        mpfr_asprintf(&str, "%.*Rf", digits, value);
        std::string result(str);
        mpfr_free_str(str);
        return result;
    }
    
    mpfr_t& get_mpfr_t() { return value; }
    const mpfr_t& get_mpfr_t() const { return value; }
    
    static void set_default_precision(int bits) {
        default_precision = bits;
        mpfr_set_default_prec(bits);
    }
};

#endif // MPFR_FLOAT_H