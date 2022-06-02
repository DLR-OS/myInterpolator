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

// fpcomp.h - Convenient floating-point comparison functions with type-
//            dependent absolute/relative error thresholds


#ifndef COMMON_FP_COMP_H
#define COMMON_FP_COMP_H

#include "fpcomp_base.h"


namespace common {


    /**
     * CFpCompTrait is a class template that provides the relative
     * and absolute error thresholds that are needed for floating-
     * point comparisons. This class canm be specialized, and will
     * form the default for the built-in type float and user-defined
     * floating point data types.
     *
     * @author Dirk "jtk" Frommholz
     * @date May 08, 2013
     */
    template <typename T> 
    class CFpCompTrait {
        public:


            /**
             * GetFpCompAbsoluteError returns the absolute error to be used
             * for the comparison of floating-point numbers.
             *
             * @return the absolute error
             */
            static const T getFpCompAbsoluteError() {
                return static_cast<T>(0.000001);
            }


            /**
             * GetFpCompRelativeError returns the relative error to be used
             * for the comparison of floating-point numbers.
             *
             * @return the relative error
             */
            static const T getFpCompRelativeError() {
                return static_cast<T>(0.000001);
            }
    };


    /************************************************************************/


    /**
     * CFpCompTrait<double> is a more accurate specialization of
     * CFpCompTrait for double-precision floating-point numbers.
     */
    template <> 
    class CFpCompTrait<double> {
        public:


            /**
             * GetFpCompAbsoluteError returns the absolute error to be used
             * for the comparison of double-precision floating-point numbers.
             *
             * @return the absolute error
             */
            static const double getFpCompAbsoluteError() {
                return 0.00000001;
            }


            /**
             * GetFpCompRelativeError returns the relative error to be used
             * for the comparison of double-precision floating-point numbers.
             *
             * @return the relative error
             */
            static const double getFpCompRelativeError() {
                return 0.00000001;
            }
    };


    /************************************************************************/

    //
    // accurate fp comparisons (absolute and relative threshold)
    //


    /**
     * Fp_equal returns true if the floating point numbers f_a and f_b
     * are considered to be equal within the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a==f_b within the error bounds
     */
    template<typename T>
    bool fp_equal(const T& f_a, const T& f_b) {
        return fp_equal<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompRelativeError(), 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_unequal returns true if the floating point numbers f_a and f_b
     * are considered to be not equal within the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a!=f_b within the error bounds
     */
    template<typename T>
    bool fp_unequal(const T& f_a, const T& f_b) {
        return !fp_equal<T>(f_a, f_b);
    }


    /************************************************************************/


    /**
     * Fp_less returns true if the floating point numbers f_a and f_b
     * are considered to be less wrt the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a<f_b within the error bounds
     */
    template<typename T>
    bool fp_less(const T& f_a, const T& f_b) {
        return fp_less<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompRelativeError(), 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_lessOrEqual returns true if the floating point numbers f_a and f_b
     * are considered to be less or equal wrt the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a<=f_b within the error bounds
     */
    template<typename T>
    bool fp_lessOrEqual(const T& f_a, const T& f_b) {
        return fp_lessOrEqual<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompRelativeError(), 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_greater returns true if the floating point numbers f_a and f_b
     * are considered to be greater wrt the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a>f_b within the error bounds
     */
    template<typename T>
    bool fp_greater(const T& f_a, const T& f_b) {
        return fp_greater<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompRelativeError(), 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_greaterOrEqual returns true if the floating point numbers f_a and f_b
     * are considered to be greater or equal wrt the absolute and relative error
     * bounds as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a>=f_b within the error bounds
     */
    template<typename T>
    bool fp_greaterOrEqual(const T& f_a, const T& f_b) {
        return fp_greaterOrEqual<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompRelativeError(), 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/
    /************************************************************************/
    /************************************************************************/


    //
    // more inaccurate but faster fp comparisons (absolute threshold only)
    //


    /**
     * Fp_equalAbsolute returns true if the floating point numbers f_a and 
     * f_b are considered to be equal within the absolute error bound only
     * as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a==f_b within the absoloute error bound only
     */
    template<typename T>
    bool fp_equalAbsolute(const T& f_a, const T& f_b) {
        return fp_equalAbsolute<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompAbsoluteError()        
        );
    }


    /************************************************************************/


    /**
     * Fp_lessAbsolute returns true if the floating point numbers f_a and f_b
     * are considered to be less wrt the absolute error bound only
     * as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a<f_b within the absolute error bound only
     */
    template<typename T>
    bool fp_lessAbsolute(const T& f_a, const T& f_b) {
        return fp_less<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_lessOrEqualAbsoloute returns true if the floating point numbers 
     * f_a and f_b are considered to be less or equal wrt the absolute error 
     * bound only as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a<=f_b within the absolute error bound
     */
    template<typename T>
    bool fp_lessOrEqualAbsolute(const T& f_a, const T& f_b) {
        return fp_lessOrEqual<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_greaterAbsoloute returns true if the floating point numbers f_a and
     * f_b are considered to be greater wrt the absolute error bound
     * as defined by the CFpCompTrait class belonging to the
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a>f_b within the error bounds
     */
    template<typename T>
    bool fp_greaterAbsolute(const T& f_a, const T& f_b) {
        return fp_greater<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /************************************************************************/


    /**
     * Fp_greaterOrEqualAbsolute returns true if the floating point numbers 
     * f_a and f_b are considered to be greater or equal wrt the absolute 
     * error bound as defined by the CFpCompTrait class belonging to the 
     * floating point type of the parameters.
     * 
     * @param f_a the first (left) number to be compared
     * @param f_b the second (right) number to be compared
     * @return true if f_a>=f_b within the error bounds
     */
    template<typename T>
    bool fp_greaterOrEqualAbsolute(const T& f_a, const T& f_b) {
        return fp_greaterOrEqual<T>(
            f_a, 
            f_b, 
            CFpCompTrait<T>::getFpCompAbsoluteError()
        );
    }


    /********************************************************************/


} // namespace common


#endif
