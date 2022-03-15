//
// This file is part of myInterpolator - image void filling
// Copyright (C) 2022  Dirk 'jtk' Frommholz, DLR OS-SEC
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.
//

// Basic floating-point comparison functions

#ifndef COMMON_FPCOMP_BASE_H
#define COMMON_FPCOMP_BASE_H


namespace common {


    /************************************************************************/


    /**
     * Fp_equal compares two floating point numbers with 
     * respect to a relative error (a factor from 0 .. 1.0 regarding 
     * the greater of the two numbers). The number type (typically
     * float or double) is specified using the template parameter.
     *
     * from:
     *
     * D.E. Knuth: The Art Of Computer Programming Vol. 2
     * Addison-Wesley 1997 
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error 
     * acceptable for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if the numbers are equal, false otherwise
     */
    template<typename T> 
    inline bool fp_equal(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {

        // absolute difference test
        const T diff=std::abs(f_a-f_b);
        if (diff<=f_absoluteError) {
            return true;
        }
 
        // absolute difference is small, do the relative difference test
        const T fabs_a=std::abs(f_a);
        const T fabs_b=std::abs(f_b);
        const T largestNumber=(fabs_b>fabs_a) ? fabs_b : fabs_a;
 
        if (diff<=largestNumber*f_relativeError) {
            return true;
        }

        // unequal
        return false;
    }


    /************************************************************************/


    /**
     * Fp_unequal works the same way as fp_equal except that it returns
     * true if the two numbers are not considered to be equal.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error acceptable
     * for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if the numbers are not equal, false otherwise
     */
    template<typename T> 
    inline bool fp_unequal(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {
        return !fp_equal(f_a, f_b, f_relativeError, f_absoluteError);
    }


    /************************************************************************/


    /**
     * Fp_less returns true if a floating-point number f_a is less
     * than and not equal to f_b.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error 
     * acceptable for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a<f_b
     */
    template<typename T> 
    inline bool fp_less(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {

        return (fp_unequal(f_a, f_b, f_relativeError, f_absoluteError) &&
        (f_a<f_b));    
    }


    /************************************************************************/


    /**
     * Fp_lessOrEqual returns true if a floating-point number f_a is less
     * or equal than f_b.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error 
     * acceptable for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a<=f_b
     */
    template<typename T> 
    inline bool fp_lessOrEqual(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {

        return (fp_equal(f_a, f_b, f_relativeError, f_absoluteError) ||
        (f_a<f_b));    
    }


    /************************************************************************/


    /**
     * Fp_greater returns true if a floating-point number f_a is greater
     * than and not equal to f_b.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error 
     * acceptable for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a>f_b
     */
    template<typename T> 
    inline bool fp_greater(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {

        return (fp_unequal(f_a, f_b, f_relativeError, f_absoluteError) &&
        (f_a>f_b));    
    }


    /************************************************************************/


    /**
     * Fp_greaterOrEqual returns true if a floating-point number f_a is 
     * greater or equal than f_b.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_relativeError the maximum relative (percentage) error 
     * acceptable for considering the numbers equal
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a>=f_b
     */
    template<typename T> 
    inline bool fp_greaterOrEqual(T f_a, T f_b, T f_relativeError, 
    T f_absoluteError) {

        return (fp_equal(f_a, f_b, f_relativeError, f_absoluteError) ||
        (f_a>f_b));    
    }


    /************************************************************************/


    /**
     * Fp_equalAbsolute compares two floating point numbers using an
     * absolute error criterion only. If the absolute value of the
     * difference of the two numbers to compare is less than this value, 
     * the numbers are considered equal. 
     *
     * This function is faster then fp_equal but does not take the 
     * amplitude of the numbers into account (that is, two large numbers
     * may not be considered equal though you expect them to be).
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if the numbers are equal, false otherwise
     */
    template<typename T> 
    inline bool fp_equalAbsolute(T f_a, T f_b, T f_absoluteError) {
        return (std::abs(f_a-f_b)<f_absoluteError);
    }


    /************************************************************************/


    /**
     * Fp_unequalAbsolute returns the opposite of fp_equalAbsolute.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if the numbers are unequal, false otherwise
     */
    template<typename T> 
    inline bool fp_unequalAbsolute(T f_a, T f_b, T f_absoluteError) {
        return fp_equalAbsolute(f_a, f_b, f_absoluteError);
    }


    /************************************************************************/


    /**
     * Fp_lessAbsolute returns true if a floating-point number f_a is less
     * than and not equal to f_b using an absolute threshold only.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a<f_b
     */
    template<typename T> 
    inline bool fp_lessAbsolute(T f_a, T f_b, T f_absoluteError) {

        return (fp_unequalAbsolute(f_a, f_b, f_absoluteError) &&
        (f_a<f_b));    
    }


    /************************************************************************/


    /**
     * Fp_lessOrEqualAbsolute returns true if a floating-point 
     * number f_a is less or equal to f_b using an absolute threshold only.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a<=f_b
     */
    template<typename T> 
    inline bool fp_lessOrEqualAbsolute(T f_a, T f_b, T f_absoluteError) {

        return (fp_equalAbsolute(f_a, f_b, f_absoluteError) ||
        (f_a<f_b));    
    }


    /************************************************************************/


    /**
     * Fp_greaterAbsolute returns true if a floating-point number f_a is 
     * greater than and not equal to f_b using an absolute threshold only.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a>f_b
     */
    template<typename T> 
    inline bool fp_greaterAbsolute(T f_a, T f_b, T f_absoluteError) {

        return (fp_unequalAbsolute(f_a, f_b, f_absoluteError) &&
        (f_a>f_b));    
    }


    /************************************************************************/


    /**
     * Fp_greaterOrEqualAbsolute returns true if a floating-point 
     * number f_a is greater or equal to f_b using an absolute threshold only.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_absoluteError the maximum absolute error acceptable for 
     * considering the number equal
     * @return true if f_a>=f_b
     */
    template<typename T> 
    inline bool fp_greaterOrEqualAbsolute(T f_a, T f_b, T f_absoluteError) {

        return (fp_equalAbsolute(f_a, f_b, f_absoluteError) ||
        (f_a>f_b));    
    }


    /************************************************************************/


} // namespace common

#endif

