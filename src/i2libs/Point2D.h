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



#ifndef DATA_POINT2D_H
#define DATA_POINT2D_H

// have C++11 ?
#if COMMON_CPLUSPLUS_VERSION > 199711L
    #include <initializer_list>
#endif

#include <cmath>
#include "basictypes.h"
#include "coolmath.h"


namespace data {


    /**
     * _CPoint2D stores 2D coordinates related to 3D objects, for instance
     * texture coordinates.
     * 
     * This is a configurable helper class and not intended for direct
     * use. Use CPoint2D below instead.
     *
     * @author Dirk "jtk" Frommholz
     * @version 1.0a
     * @date Feb 28, 2007
     *
     * Revisions:
     *
     * 1.0a
     * - initial revision
     *
     */
    template <typename _CoordType>
    struct _CPoint2D {


        /**
         * Alias of the coordinate type for easy access using the class name.
         */
        typedef _CoordType CoordType;


        /**
         * The horizontal coordinate of the point.
         */
        CoordType m_i;


        /**
         * The vertical coordinate of the point.
         */
        CoordType m_j;


        /**
         * The default constructor.
         */
        _CPoint2D():m_i(0), m_j(0) { } 


        /**
         * The convenience constructor.
         *
         * @param f_i the horizontal coordinate
         * @param f_j the vertical coordinate
         */ 
        _CPoint2D(CoordType f_i, CoordType f_j):
        m_i(f_i), m_j(f_j) { }


        #if COMMON_CPLUSPLUS_VERSION > 199711L

            /**
             * The constructor using C++11 initializer lists. Not efficient
             * but convenient.
             *
             * @param f_il the initializer list to set up the 2D point with;
             * if less than three elements are passed those coordinates not
             * provided will be set to zero
             */ 
            _CPoint2D(std::initializer_list<CoordType> f_il) {    
                m_i=(f_il.size()>0) ? *(f_il.begin()) : 0;
                m_j=(f_il.size()>1) ? *(f_il.begin()+1) : 0;
            }

        #endif


        /**
         * SetCoords sets the two coordinates.
         *           
         * @param f_i the horizontal coordinate
         * @param f_j the vertical coordinate
         */
        void setCoords(CoordType f_i, CoordType f_j) {
            m_i=f_i; 
            m_j=f_j;
        }


        /**
         * AtIndex returns the first or second read-only coordinate 
         * according to the index f_index which must be 0 or 1 respectively.
         * Indices that are out of range are wrap around using the modulus
         * function.
         *
         * @param f_index the index that will be translated to one of the 
         * point's coordinates
         * @return the read-only coordinate belonging to the index
         */
        const CoordType& atIndex(unsigned f_index) const {
            if (f_index & 1) {
                return m_j;
            } else {
                return m_i;
            }
        }


        /**
         * AtIndex returns the first or second coordinate according to the
         * index f_index which must be 0 or 1 respectively. Indices
         * that are out of range are wrap around using the modulus
         * function.
         *
         * @param f_index the index that will be translated to one of the 
         * point's coordinates
         * @return the coordinate belonging to the index
         */
        CoordType& atIndex(unsigned f_index) {
            if (f_index & 1) {
                return m_j;
            } else {
                return m_i;
            }
        }


        /**
         * Operator-= subtract f_subtrahend from this point.
         * 
         * @param f_subtrahend the right hand side point to be subtracted 
         * from this one
         * @return this point after subtraction
         */
        _CPoint2D operator-=(const _CPoint2D& f_subtrahend) {
            m_i-=f_subtrahend.m_i;
            m_j-=f_subtrahend.m_j;
            return *this;
        }


        /**
         * Operator- subtracts a 2-space point from this one by subtracting 
         * its components.
         *
         * @param f_subtrahend the 2d point to be subtracted
         * @return the difference of the two points     
         */
        _CPoint2D operator-(const _CPoint2D& f_subtrahend) const {
            return _CPoint2D(
                m_i-f_subtrahend.m_i,
                m_j-f_subtrahend.m_j
            );
        }


        /**
         * Operator- (unary) inverts the sign of each component of the point.     
         *
         * @return the negated 2-space element
         */
        _CPoint2D operator-() const {
            return _CPoint2D(-m_i, -m_j);
        }


        /**
         * Operator+= adds f_summand from this point.
         * 
         * @param f_summand the right hand side point to be added to
         * this one
         * @return this point after addition
         */
        _CPoint2D operator+=(const _CPoint2D& f_summand) {
            m_i+=f_summand.m_i;
            m_j+=f_summand.m_j;
            return *this;
        }


        /**
         * Operator+ adds a 2-space point from this one by adding 
         * its components.
         *
         * @param f_summand the 2d point to be added
         * @return the sum of the two points     
         */
        _CPoint2D operator+(const _CPoint2D& f_summand) const {
            return _CPoint2D(
                m_i+f_summand.m_i,
                m_j+f_summand.m_j
            );
        }


        /**
         * Operator* multiplies the 2-space vector with the
         * factor element-by-element.
         *
         * @param f_factor the factor to be multiplied
         * @return the scaled vector
         */
        _CPoint2D operator*(const _CoordType f_factor) const {
            return _CPoint2D(
                m_i*f_factor,                
                m_j*f_factor
            );        
        }


        /**
         * Operator* computes the scalar product (dot product)
         * of two points.
         *
         * @param f_p the right-hand side point
         * @return the scalar product
         */
        CoordType operator*(const _CPoint2D f_p) const {
            return (m_i*f_p.m_i)+(m_j*f_p.m_j);              
        }


        /**
         * Multiply multiplies the 2-space vector with the
         * given vector element-by-element.
         *
         * @param f_factor the factor to be multiplied
         * @return the scaled vector
         */
        _CPoint2D multiply(const _CPoint2D& f_rhs) const {
            return _CPoint2D(
                m_i*f_rhs.m_i,                
                m_j*f_rhs.m_j
            );        
        }


        /**
         * Operator/ divides the 2-space vector with the
         * divisor element-by-element.
         *
         * @param f_divisor the divisor to be used
         * @return the divided vector
         */
        _CPoint2D operator/(const _CoordType f_divisor) const {
            return _CPoint2D(
                m_i/f_divisor,                
                m_j/f_divisor
            );        
        }


        /**
         * Operator/ performs element-wise division of the vector elements.
         *
         * @param f_divisor the 2D point for element-wise division
         * @return a new 2D point with its coordinates divided by the
         * respective coordinates of f_divisor
         */
         _CPoint2D operator/(const _CPoint2D& f_divisor) const {

            return _CPoint2D(
                m_i/f_divisor.m_i,
                m_j/f_divisor.m_j
            );
        }


        /**
         * SqLength returns the squared length of the vector.
         * Useful for comparison purposes as it avoids computing the
         * square root.
         *
         * @return the length of the vector
         */
        CoordType sqLength() const {
            return m_i*m_i+m_j*m_j;
        }


        /**
         * Length returns the Euclidean length of the vector.
         *
         * @return the length of the vector
         */
        CoordType length() const {
            return ::sqrt(m_i*m_i+m_j*m_j);
        }


        /**
         * Normalize normalizes this point to unit length.
         *
         * @return the normalized point
         */
        _CPoint2D normalize() const {
            return *this/length();
        }
        

        /**
         * The <-operator returns true if this point is less than the
         * point passed. Less means that the y coordinate needs to be less
         * first, and if it is equal, the x coordinate decides.
         *
         * @param f_p the point to be hashed
         */
        bool operator<(const _CPoint2D& f_rhs) const {

            if (std::is_integral<_CoordType>::value) {
                if (m_j<f_rhs.m_j) {        
                    return true;
                };
    
                if (m_j==f_rhs.m_j) {                
                    return m_i<f_rhs.m_i;
                }
                        
                return false;

            } else {
                if (common::fp_less(static_cast<double>(m_j), 
                static_cast<double>(f_rhs.m_j))) {            
                    return true;
                };
    
                if (common::fp_equal(static_cast<double>(m_j), 
                static_cast<double>(f_rhs.m_j))) {
                    return common::fp_less(static_cast<double>(m_i), 
                    static_cast<double>(f_rhs.m_i));
                }

            
                return false;
            }
        }


        /**
         * Operator== returns true if and only if the 3-space point
         * passed equals this point. The three coordinates are
         * compared bit by bit for integral coord types and using fp_equal
         * for non-integral coord types.
         *
         * @param f_rhs the 3-space point to compare this one to
         * @return true on equality, false else
         */
        bool operator==(const _CPoint2D& f_rhs) const {
            return std::is_integral<_CoordType>::value ?
                ((m_i==f_rhs.m_i) && (m_j==f_rhs.m_j)) :
                (common::fp_equal(static_cast<double>(m_i), 
                static_cast<double>(f_rhs.m_i)) && 
                common::fp_equal(static_cast<double>(m_j), 
                static_cast<double>(f_rhs.m_j)));        

        }


        /**
         * Operator== returns true if and only if the 3-space point
         * passed equals this point. The three coordinates are
         * compared bit by bit.
         *
         * @param f_rhs the 2-space point to compare this one to
         * @return true on equality, false else
         */
        bool operator!=(const _CPoint2D& f_rhs) const {
            return !(*this==f_rhs);
        }


        /**
         * Operator() implements a hash function object for unordered_map.
         *
         * @param f_p the point to be hashed
         */
        size_t operator()(const _CPoint2D& f_p) const {
            return std::is_integral<_CoordType>::value ?
                static_cast<size_t>(f_p.m_i^f_p.m_j) :
                static_cast<size_t>(f_p.m_i*f_p.m_j);
        }


        /**
         * Abs removes the sign of both coordinates.
         *
         * @return a copy of the point with non-negative coordinates only
         */
        _CPoint2D abs() const {
            return std::is_integral<_CoordType>::value ?
                _CPoint2D(::abs(m_i), ::abs(m_j)) :
                _CPoint2D(fabs(static_cast<double>(m_i)), 
                fabs(static_cast<double>(m_j)));
        }

    };


    /**
     * CPoint2D stores 2D coordinates of for instance textures that are 
     * going to be mapped to 3D objects. For explanations, see _CPoint2D, 
     * but use this definition instead of _CPoint2D.
     */
    typedef _CPoint2D<sample_type> CPoint2D;


    /**
     * CPoint2DInt stores integer 2D coordinates of for 
     * instance texture dimensions. For explanations, see _CPoint2D, 
     * but use this definition instead of _CPoint2D. Note that some
     * functions of _CPoint2D might not make sense with integer points.
     */
    typedef _CPoint2D<pos_type> CPoint2DInt;


} //namespace data




#endif
