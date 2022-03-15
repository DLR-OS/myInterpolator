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

/*
 * coolmath.h - cmath routines with extensions
 *
 */

#ifndef COMMON_COOLMATH_H
#define COMMON_COOLMATH_H

#include <cmath>
#include <cstddef>
#include <algorithm>

#include "fpcomp.h"

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

namespace common {


    //
    // MIN/MAX FUNCTIONS
    //

    /**
     * Min2 returns the minimum of 2 numbers. It is designed to
     * work with integer types.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @return the maximum of f_a and f_b
     */
    template <typename T> inline T min2(T f_a, T f_b) {
        return (f_a<f_b) ? f_a : f_b;
    }


    /************************************************************************/


    /**
     * Max2 returns the maximum of 2 numbers. It is designed to
     * work with integer types.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @return the maximum of f_a and f_b
     */
    template <typename T> inline T max2(T f_a, T f_b) {
        return (f_a>f_b) ? f_a : f_b;
    }


    /************************************************************************/


    /**
     * Min2Index returns the number (0 or 1) of the minimum valued parameter.
     * If both values are equal, it returns 0.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @return  0 if f_a<=f_b, 1 otherwise
     */
    template <typename T>
    inline int min2Index(T f_a, T f_b) {
        return (f_a<=f_b) ? 0 : 1;
    }


    /************************************************************************/


    /**
     * Max2Index returns the number of the parameter that has the
     * maximum value.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @return 0 is f_a>=f_b, 1 otherwise
     */
    template <typename T> inline int max2Index(T f_a, T f_b) {
        return (f_a>=f_b) ? 0 : 1;
    }


    /************************************************************************/


    /**
     * Min3 returns the minimum of 3 numbers. It is designed to
     * work with integer types.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @return the minimum of f_a, f_b and f_c
     */
    template <typename T> inline T min3(T f_a, T f_b, T f_c) {
        return ((f_a < f_b) ? (( f_a < f_c) ? f_a : f_c) :((f_b < f_c) ? 
        f_b : f_c));
    }


    /************************************************************************/


    /**
     * Min3Index returns the number of the parameter that is the smallest
     * where counting starts at zero. 
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @return 0 if f_a, 1 if f_b and 2 if f_c is the smallest number
     */
    template <typename T> inline unsigned int min3Index(T f_a, T f_b, T f_c) {
        return ((f_a <= f_b) ? (( f_a <= f_c) ? 0 : 2) :((f_b <= f_c) ? 
        1 : 2));
    }


    /************************************************************************/


    /**
     * Min4 returns the minimum of 4 numbers. It is designed to
     * work with integer types, i.e. normal comparisons are used.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @param f_d the fourth number
     * @return the minimum of f_a, f_b, f_c and f_d
     */
    template <typename T> inline T min4(T f_a, T f_b, T f_c, T f_d) {
        return min2(f_a, min3(f_b, f_c, f_d));
    }


    /************************************************************************/


    /**
     * Min4Index returns the number of the parameter that is the smallest
     * where counting starts at zero. 
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @param f_d the fourth number
     * @return 0 if f_a, 1 if f_b, 2 if f_c or 3 if f_d is the smallest 
     * number
     */
    template <typename T> inline unsigned int min4Index(T f_a, T f_b, T f_c, 
    T f_d) {

        switch (min3Index(f_a, f_b, f_c)) {
            case 0:
                return f_d<f_a ? 3 : 0;
            case 1:
                return f_d<f_b ? 3 : 1;
            default:
                return f_d<f_c ? 3 : 2;
        }
    }


    /************************************************************************/


    /**
     * Max3 returns the maximum of 3 numbers. It is designed to
     * work with integer types.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @return the maximum of f_a, f_b and f_c
     */
    template <typename T> inline T max3(T f_a, T f_b, T f_c) {
        return ((f_a > f_b) ? (( f_a > f_c) ? f_a : f_c) :((f_b > f_c) ? 
        f_b : f_c));
    }


    /************************************************************************/


    /**
     * Max3Index returns the index of the parameter that has the maximum 
     * value. On equality, the index of the first equal parameter set is 
     * returned. 
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @return 0 if f_a, 1 if f_b or 2 if f_c is the maximum parameter
     */
    template <typename T> inline int max3Index(T f_a, T f_b, T f_c) {
        return ((f_a >= f_b) ? (( f_a >= f_c) ? 0 : 2) :((f_b >= f_c) ? 
        1 : 2));
    }


    /************************************************************************/


    /**
     * Max4 returns the maximum of 4 numbers. It is designed to
     * work with integer types, i.e. normal comparisons are used.
     *
     * @param f_a the first number
     * @param f_b the second number
     * @param f_c the third number
     * @param f_d the fourth number
     * @return the maximum of f_a, f_b, f_c and f_d
     */
    template <typename T> inline T max4(T f_a, T f_b, T f_c, T f_d) {
        return max2(f_a, max3(f_b, f_c, f_d));
    }


    /************************************************************************/


    /**
     * Sort3 sorts three parameters which must support the < operation.
     *
     * @param f_a the first datum
     * @param f_b the second datum
     * @param f_c the third datum 
     */
    template <typename T> inline void sort3(T& f_a, T& f_b, T& f_c) {

        // three-element bubble sort
        if (!(f_a<f_b)) {
            std::swap(f_a, f_b);
        }

        if (!(f_b<f_c)) {
            std::swap(f_b, f_c);
        }

        if (!(f_a<f_b)) {
            std::swap(f_a, f_b);
        }
    }


    /************************************************************************/


    //
    // ANGLE FUNCTIONS
    //


    /**
     * Rad2Deg converts radians to degrees.
     *
     * @param f_radAngle the angle in radians
     * @return the angle in degrees
     */
    double rad2deg(double f_radAngle);    


    /************************************************************************/


    /**
     * Deg2Rad converts degrees to radians.
     *
     * @param f_degAngle the angle in degrees
     * @return the angle in radians
     */
    double deg2rad(double f_degAngle);


    /************************************************************************/


    //
    // MISC MATH FUNCTIONS
    //


    /**
     * RoundNat implements natural rounding rounding, i.e. 45.5 -> 46 
     * and -45.5 -> -46.
     *
     * @param f_a the number to be rounded, must be of a built-in floating 
     * point type (float or double)
     */
    template <typename T> T roundNat(T f_a) {
        T truncPart;

        // get fractional part of f_a with sign of f_a
        T fracPart=std::modf(f_a, &truncPart);

        // double fractional part to see where to round towards,
        // and round it (i.e. 2*0.5 -> 1, rounded to 1 and 2*-0.5=-1
        // rounded to -1);
        T roundFracPart;
        std::modf(fracPart+fracPart, &roundFracPart);

        // addition works losslessy on IEEE754 fp types
        return truncPart+roundFracPart;
    }


    /************************************************************************/


    /**
     * Sgn returns 1 if the number is positive, -1 if it is negative and
     * 0 if it is zero (the so called signum). The implementation is
     * branch-less.
     *
     * @param f_a the number to compute the signum for
     * @return the signum of the number
     */
    template <typename T>
    inline int sgn(T f_a) {        
        return (T(0)<f_a)-(f_a<T(0));
    }


    /************************************************************************/


    /**
     * Square returns the square of a number.
     *
     * @return the square of a number
     */
    template<typename T>
    T square(T f_x) {
        return f_x*f_x;
    }


    /************************************************************************/


    /**
     * Cubic returns the power of 3 of a number.
     *
     * @return the square of a number
     */
    template<typename T>
    T cubic(T f_x) {
        return f_x*f_x*f_x;
    }


    /************************************************************************/


    /**
     * GemanMcClure computes the robust error function named after their
     * inventors (x^2/(s^2+x^2)).
     *
     * @param f_s determines where the function begins to descend smoothly
     * @param f_x the sample
     * @return the GemanMcClure value
     */
    double gemanMcClure(double f_s, double f_x);


    /************************************************************************/


    /**
     * GaussianDistribution computes the Gaussian normal distribution at f_x with
     * the expectation value f_mu and the variance f_sigma.
     *
     * @param f_x the function parameter
     * @param f_mu the expectation value
     * @param f_sigma the variance
     * @return the Gaussian normal distribution f(x, mu, sigma)
     */
    double gaussianDistribution(double f_x, double f_mu, double f_sigma);


    /************************************************************************/


    /**
     * GaussianCDF computes the cumulative distribution function of the
     * Gaussian normal distributiuon N(0,1).
     *
     * @param f_x the function parameter, an arbitrary number
     * @return the value of the CDF for f_x
     */
    double gaussianCDF(double f_x);


    /************************************************************************/


    /**
     * PeriodIncluding maps an integer value x to the periodic range 
     * [0..max], i.e. x=-1 will become 1, x=max+1 will become max-1 and so 
     * on, effectively implementing a shifted periodic discrete triangle 
     * function without repetition at the borders.\n\n
     *
     * T must be a signed integer type.
     *
     * @param f_x the value to be periodized
     * @param f_max the upper bound of the period, including f_max
     */
    template<typename T> T periodIncluding(const T& f_x, const T& f_max) {
        return f_max-abs(abs(f_x%(f_max<<1))-f_max);
    }


    /************************************************************************/


    /**
     * Period maps an integer value x to the periodic range [0..max), i.e.
     * x=-1 will become 1, x=max+1 will become max-1 and so on, effectively
     * implementing a shifted periodic discrete triangle function without
     * repetition at the borders.\n\n
     *
     * T must be a signed integer type.
     *
     * @param f_x the value to be periodized
     * @param f_max the upper bound of the period, excluding, i.e. the
     * maximum value of the sequence will be f_max-1
     */
    template<typename T> T period(const T& f_x, const T& f_max) {
        return periodIncluding(f_x, f_max-1);
    }


    /************************************************************************/


    /**
     * Bounce maps an integer value x to the periodic range [0..max) with
     * repetition at the borders, i.e. the sequenece -2 -1 0 1 2 3 4 5 6 7 8
     * will become -1 0 0 1 2 3 3 2 1 0 0 with max=4.\n\n
     *
     * T must be a signed integer type.
     *
     * @param f_x the value to be periodized
     * @param f_max the upper bound of the period (excluding, i.e. the
     * maximum value in the bounced sequence will be f_max-1)
     */
    template<typename T> T bounce(const T& f_x, const T& f_max) {
        if (f_x>=0) { 
            return (f_x/f_max) & 1 ? 
                f_max-1-(f_x%f_max) :
                (f_x%f_max);
        } else {
            return ((f_x+1)/f_max) & 1 ? 
                f_max-1+((f_x+1)%f_max) :
                -((f_x+1)%f_max);
        }
    }


    /************************************************************************/


    /**
     * FixZeroSign turns a negative zero into an unsigned zero. This function
     * works if for the data type T the value -0 equals 0 like for IEEE754
     * floating-point values. For typical integer encodings this function
     * does nothing.
     *
     * @tparam T the data type of the value to be converted
     * @param f_a the value to be converted 
     * @return zero if f_a is minus zero, any other values remain unchanged
     */
    template<typename T> T fixZeroSign(const T& f_a) {
        return (f_a==static_cast<T>(-0.0)) ? static_cast<T>(0) : f_a;
    }


    /************************************************************************/


    /**
     * Clip ensures that the value f_a is inside the range f_low...f_high
     * including the bounds.
     *
     * @tparam T the type of the value and the bounds
     * @param f_a the value to be adjusted if necessary
     * @param f_low the lower bound
     * @param f_high the higher bound
     */
    template<typename T> T clip(const T& f_a, const T& f_low, const 
    T& f_high) {
        return f_a<f_high ? (f_a>f_low ? f_a : f_low) : f_high;
    }


    /************************************************************************/


    

} // namespace common

#endif


