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

#ifndef DATA_POINT3D_H
#define DATA_POINT3D_H

#include <cmath>
#include "basictypes.h"
#include "coolmath.h"

// have C++11 ?
#if COMMON_CPLUSPLUS_VERSION > 199711L
    #include <initializer_list>
#endif


namespace data {


    /**
     * _CPoint3D stores the coordinates of a point in 3-space.
     * This is a configurable helper class and not intended for direct
     * use. Use CPoint3D below instead.
     *
     * @author Dirk "jtk" Frommholz
     * @version 1.0b
     * @date Jan 29, 2007
     *
     * Revisions:
     *
     * 1.1
     * - removed color values as they are not needed anymore
     *
     * 1.0a
     * - initial revision
     *
     */
    template <typename _CoordType>
    struct _CPoint3D {


        /**
         * Alias of the coordinate type for easy access using the class name.
         */
        typedef _CoordType CoordType;


        /**
         * The x coordinate of the point.
         */
        CoordType m_x;


        /**
         * The y coordinate of the point.
         */
        CoordType m_y;


        /**
         * The z coordinate of the point.
         */
        CoordType m_z;    


        /**
         * The constructor.
         *
         * @param f_x the x coordinate
         * @param f_y the y coordinate
         * @param f_z the z coordinate
         */ 
        _CPoint3D(CoordType f_x=0, CoordType f_y=0, CoordType f_z=0) {    
            m_x=f_x; 
            m_y=f_y;
            m_z=f_z;    
        }


        #if COMMON_CPLUSPLUS_VERSION > 199711L

            /**
             * The constructor using C++11 initializer lists. Not efficient
             * but convenient.
             *
             * @param f_il the initializer list to set up the 3D point with;
             * if less than three elements are passed those coordinates not
             * provided will be set to zero
             */ 
            _CPoint3D(std::initializer_list<CoordType> f_il) {    
                m_x=(f_il.size()>0) ? *(f_il.begin()) : 0;
                m_y=(f_il.size()>1) ? *(f_il.begin()+1) : 0;
                m_z=(f_il.size()>2) ? *(f_il.begin()+2) : 0;
            }

        #endif


        /**
         * SetCoords sets all three coordinates at once.
         *     
         * @param f_x the x coordinate
         * @param f_y the y coordinate
         * @param f_z the z coordinate
         */
        void setCoords(CoordType f_x, CoordType f_y, CoordType f_z) {
            m_x=f_x; 
            m_y=f_y;
            m_z=f_z;
        }


        /**
         * AtIndex returns the x, y or z coordinate according to the
         * index f_index which must be 0, 1 or 2 respectively. Indices
         * that are out of range are wrap around using the modulus
         * function.
         *
         * @param f_index the index that will be translated to one of the 
         * point's coordinates
         * @return the coordinate belonging to the index
         */
        CoordType& atIndex(unsigned f_index) {
            switch (f_index % 3) {
                case 0:
                    return m_x;
                break;
                case 1:
                    return m_y;
                break;
                default:
                    return m_z;
            }
        }


        /**
         * AtIndex (const version) is the read-only version of the above.
         *
         * @param f_index the index that will be translated to one of the 
         * point's coordinates
         * @return the read-only coordinate belonging to the index
         */
        const CoordType& atIndex(unsigned f_index) const {
            switch (f_index % 3) {
                case 0:
                    return m_x;
                break;
                case 1:
                    return m_y;
                break;
                default:
                    return m_z;
            }
        }


        /**
         * Operator+ adds a 3-space point to this one by adding its components.
         *
         * @param f_summand the 3d point to be added
         * @return the sum of the two points     
         */
        _CPoint3D operator+(const _CPoint3D& f_summand) const {
            return _CPoint3D(
                m_x+f_summand.m_x,                
                m_y+f_summand.m_y,
                m_z+f_summand.m_z
            );        
        }


        /**
         * Operator- subtracts a 3-space point from this one by subtracting 
         * its components.
         *
         * @param f_subtrahend the 3d point to be subtracted
         * @return the difference of the two points     
         */
        _CPoint3D operator-(const _CPoint3D& f_subtrahend) const {
            return _CPoint3D(
                m_x-f_subtrahend.m_x,                
                m_y-f_subtrahend.m_y,
                m_z-f_subtrahend.m_z
            );        
        }


        /**
         * Operator- (unary) inverts the sign of each component of the point.     
         *
         * @return the negated 3-space element
         */
        _CPoint3D operator-() const {
            return _CPoint3D(
                -m_x,                
                -m_y,
                -m_z
            );        
        }


        /**
         * Sin computes the sine for each element.
         *
         * @return the sine'd point
         */
        _CPoint3D sin() const {
            return _CPoint3D(
                ::sin(m_x),
                ::sin(m_y),
                ::sin(m_z)
            );
        }


        /**
         * Cos computes the cosine for each element.
         *
         * @return the cosine'd point
         */
        _CPoint3D cos() const {
            return _CPoint3D(
                ::cos(m_x),
                ::cos(m_y),
                ::cos(m_z)
            );
        }


        /**
         * Operator+= adds another point to this one and replaces its
         * elements by the sum.
         *
         * @param f_p the point to be added
         * @return this point, replaced by the sum
         */
        _CPoint3D& operator+=(const _CPoint3D& f_p) {
            m_x+=f_p.m_x;
            m_y+=f_p.m_y;
            m_z+=f_p.m_z;
            return *this;
        }


        /**
         * Operator-= subtracts another point from this one and replaces the
         * elements of this pointer by the difference.
         *
         * @param f_p the point to be subtracted
         * @return this point, replaced by the difference
         */
        _CPoint3D& operator-=(const _CPoint3D& f_p) {
            m_x-=f_p.m_x;
            m_y-=f_p.m_y;
            m_z-=f_p.m_z;
            return *this;
        }


        /**
         * RotX rotates the vector around the x axis using the specified
         * angle in radians.
         *
         * @param f_angle the angle in radians
         * @return the vector rotated around the x axis
         */
        _CPoint3D rotX(double f_angle) const {
            return _CPoint3D(
                m_x,
                m_y*::cos(f_angle)-m_z*::sin(f_angle),
                m_y*::sin(f_angle)+m_z*::cos(f_angle)           
            );
        }


        /**
         * RotY rotates the vector around the y axis using the specified
         * angle in radians.
         *
         * @param f_angle the angle in radians
         * @return the vector rotated around the y axis
         */
        _CPoint3D rotY(double f_angle) const {
            return _CPoint3D(
                m_x*::cos(f_angle)+m_z*::sin(f_angle),
                m_y,
               -m_x*::sin(f_angle)+m_z*::cos(f_angle)
            );
        }


        /**
         * RotZ rotates the vector around the z axis using the specified
         * angle in radians.
         *
         * @param f_angle the angle in radians
         * @return the vector rotated around the z axis
         */
        _CPoint3D rotZ(double f_angle) const {
            return _CPoint3D(
                m_x*::cos(f_angle)-m_y*::sin(f_angle),
                m_x*::sin(f_angle)+m_y*::cos(f_angle),
                m_z
            );
        }


        /**
         * Operator*= performs element-wise multiplication of the vector elements
         * with a real number and replaces this vector by the result.
         *
         * @param f_factor the real-valued factor     
         * @return this vector after scaling
         */
         _CPoint3D& operator*=(CoordType f_factor) {        
            m_x*=f_factor;
            m_y*=f_factor;
            m_z*=f_factor;

            return *this;        
        }


        /**
         * Operator* performs element-wise multiplication of all vector elements
         * with a real number.     
         *
         * @param f_factor the real-valued factor
         * @return the scaled vector
         */
         _CPoint3D operator*(CoordType f_factor) const {
            return _CPoint3D(
                m_x*f_factor,
                m_y*f_factor,
                m_z*f_factor
            );
        }


        /**
         * Operator/ performs element-wise division of all vector elements
         * with a real number.     
         *
         * @param f_divisor the real-valued non-zero divisor
         * @return the downscaled vector
         */
         _CPoint3D operator/(CoordType f_divisor) const {
            return _CPoint3D(
                m_x/f_divisor,
                m_y/f_divisor,
                m_z/f_divisor
            );
        }


        /**
         * Operator* returns the scalar product of this vector and the specified
         * right side.
         *
         * @param f_right the real-valued right-side vector
         * @return the scalar product of this vector and f_right
         */
         CoordType operator*(const _CPoint3D& f_right) const {
            return (m_x*f_right.m_x)+(m_y*f_right.m_y)+(m_z*f_right.m_z);
        }


        /**
         * Multiply returns the vector that results from element-wise 
         * multiplication with another vector.
         *
         * @param f_right the real-valued right-side vector
         * @return this vector elementwisely multiplied with f_right
         */
         _CPoint3D multiply(const _CPoint3D& f_right) const {
            return _CPoint3D(
                m_x*f_right.m_x,
                m_y*f_right.m_y,
                m_z*f_right.m_z
            );        
        }


        /**
         * Divide returns the vector that results from element-wise 
         * division with another vector.
         *
         * @param f_right the real-valued right-side vector
         * @return this vector elementwisely divided by f_right
         */
         _CPoint3D divide(const _CPoint3D& f_right) const {
            return _CPoint3D(
                m_x/f_right.m_x,
                m_y/f_right.m_y,
                m_z/f_right.m_z
            );        
        }


        /**
         * Abs makes all point coordinates positive.
         *
         * @return this vector elementwisely divided by f_right
         */
         _CPoint3D abs() const {
            return _CPoint3D(
                fabs(m_x),
                fabs(m_y),
                fabs(m_z)
            );        
        }


        /**
         * Length returns the length of the vector.
         *
         * @return the length of the vector
         */
        CoordType length() const {
            return sqrt(m_x*m_x+m_y*m_y+m_z*m_z);
        }


        /**
         * Normalize scales the length of the vector to f_length
         * without changing its direction. However, for non-floating point
         * vectors the length may get truncated.
         *
         * @param f_length the new length, defaults to 1.0
         * @return the length of the vector
         */
        _CPoint3D normalize(double f_length=1) const {
            return (*this)*(f_length/(*this).length());        
        }


        /**
         * SqLength returns the squared length of the vector.
         * Useful for comparison purposes as it avoids computing the
         * square root.
         *
         * @return the length of the vector
         */
        CoordType sqLength() const {
            return m_x*m_x+m_y*m_y+m_z*m_z;
        }


        /**
         * Cross computes the cross product of this vector with
         * another one.
         *     
         * @param f_b the vector to compute the cross product with
         * @return the cross product of this vector and f_vec (*this x f_vec)
         */
        _CPoint3D cross(const _CPoint3D& f_b) const {
            return _CPoint3D(
                m_y*f_b.m_z-m_z*f_b.m_y,
                m_z*f_b.m_x-m_x*f_b.m_z,
                m_x*f_b.m_y-m_y*f_b.m_x
            );
        }


        /**
         * Operator() implements a hash function object for unordered_map
         * will work sufficiently for integers only
         *
         * @param f_p the point to be hashed
         */
        std::size_t operator()(const _CPoint3D& f_p) const {
            return static_cast<std::size_t>(1187*f_p.m_x+571*f_p.m_y+151*
            f_p.m_z);
        }


        /**
         * Operator== returns true if and only if the 3-space point
         * passed equals this point. Depending on the coordinate type
         * either the exact or floating-point comparisons are used.
         *
         * @param f_rhs the 3-space point to compare this one to
         * @return true on equality, false else
         */
        bool operator==(const _CPoint3D& f_rhs) const {
            return std::is_integral<_CoordType>::value ?
                ((m_x==f_rhs.m_x) && (m_y==f_rhs.m_y) && (m_z==f_rhs.m_z)) :

                // casts are necessary to suppress annoying compile-time 
                // warnings in case of int's
                (common::fp_equal<_CoordType>(static_cast<_CoordType>(m_x), 
                static_cast<_CoordType>(f_rhs.m_x)) && 
                common::fp_equal<_CoordType>(static_cast<_CoordType>(m_y), 
                static_cast<_CoordType>(f_rhs.m_y)) && 
                common::fp_equal<_CoordType>(static_cast<_CoordType>(m_z), 
                static_cast<_CoordType>(f_rhs.m_z)));        
        }


        /**
         * Operator!= returns true if and only if the 3-space point
         * passed does not equal this point. This is operator== 
         * negated.
         *
         * @param f_rhs the 3-space point to compare this one to
         * @return true on inequality, false else
         */
        bool operator!=(const _CPoint3D& f_rhs) const {
            return !(*this==f_rhs);
        }


        /**
         * Operator< returns true if and only if this 3-space point
         * is lexicographically strictly less than the right hand side point.
         *
         * @param f_rhs the 3-space point to compare this one to
         * @return true if this point is lexicographically less then f_rhs
         */
        bool operator<(const _CPoint3D& f_rhs) const {
            if (std::is_integral<_CoordType>::value) {
                if (m_x!=f_rhs.m_x ) {
                    return m_x<f_rhs.m_x;
                }

                if (m_y!=f_rhs.m_y) {
                    return m_y<f_rhs.m_y;
                }

                return m_z<f_rhs.m_z;
            }

            // have fp coordinates, use consistent fp_equal family of
            // comparisons
            if (common::fp_unequal<_CoordType>(static_cast<_CoordType>(m_x),
            static_cast<_CoordType>(f_rhs.m_x))) {
                return m_x<f_rhs.m_x;
            }

            if (common::fp_unequal<_CoordType>(static_cast<_CoordType>(m_y),
            static_cast<_CoordType>(f_rhs.m_y))) {
                return m_y<f_rhs.m_y;
            }

            return common::fp_less<_CoordType>(static_cast<_CoordType>(m_z),
            static_cast<_CoordType>(f_rhs.m_z));
        }


        /**
         * FixZeroSign fixes the negative zero and turns it into a
         * "positive" zero. This function only makes sense if the
         * CoordType follows an encoding scheme where -0 == 0, i.e. 
         * IEE754 and is meaningless on integer numbers.
         *
         * @return this 3D point with any negative zero turned into a
         * "positive" one, i.e. -0 becomes 0
         */
        _CPoint3D fixZeroSign() const {
            return _CPoint3D(
                common::fixZeroSign(m_x),
                common::fixZeroSign(m_y),   
                common::fixZeroSign(m_z)
            );
        }

    };


    /**
     * CPoint3D stores the floating point coordinates of a point in 3-space.
     * For explanations, see _CPoint3D, but use this definition instead of
     * _CPoint3D.
     */
    typedef _CPoint3D<sample_type> CPoint3D;



    /**
     * CPoint3DInt stores the unsigned int coordinates of a point in
     * a discrete signed integer 3-space. For explanations, see _CPoint3D,
     but use this definition instead of _CPoint3D.
     */
    typedef _CPoint3D<pos_type> CPoint3DInt;


} // namespace data


#endif
