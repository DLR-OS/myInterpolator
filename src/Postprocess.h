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

#ifndef IMAGE_POSTPROCESS_H
#define IMAGE_POSTPROCESS_H

#include "Image.h"
#include "Progress.h"

namespace image {


    /**
     * CPostprocess smoothens the interpolated pixels.
     *
     * @author Dirk 'jtk' Frommholz
     * @date March 03, 2022
     */
    class CPostprocess {
        public:


            /**
             * MakeBinaryMask creates a binary mask image for a floating-
             * point mask image as a more compact representation (typically 
             * 1/32 of original space allocated). The binary masks can be 
             * used with the filter functions below.
             *
             * @param f_fpMask the floating-point source mask image with
             * interpolated pixels being indicated by positive values
             * @param f_bitMask the binary single-channel representation 
             * (1 bit per pixel) of f_fpMask where set bits indicate pixels 
             * with a positive value in f_fpMask's first channel; 
             * will be entirely set up by this function and have 
             * f_fpMask.width/8 by f_fpMask.height pixels
             */
            void makeBinaryBitMask(const CFloatImage& f_fpMask, CU8Image& 
            f_bitMask);


            /**
             * FilterGauss applies a 2D Gaussian to the source image
             * and writes the result to the destination image. The pixels
             * to be filtered are those with a non-zero bit in the compact
             * binary f_bitMasklImage. This function is implemented as a
             * multithreaded straightforward 2D convolution.
             *
             * @param f_srcDestImage the image to be filtered, also holds the
             * result afterwards
             * @param f_bitMaskImage the single-channel compact binary image 
             * describing the pixels to be filtered (set bits)
             * @param f_tempImage a temporary image of the same 
             * characteristics as the source/destination image holding
             * the intermediate output of horizontal filtering
             * @param f_filterSize the horizontal and vertical filter
             * size, will be adjusted to positive odd integers if necessary
             * @param f_progress_p the optional progress indicator
             */
            void filterGauss(CFloatImage& f_srcDestImage, const CU8Image& 
            f_bitMaskImage, image::CFloatImage& f_tempImage, const data::
            CPoint2DInt& f_filterSize, common::CProgress* const 
            f_progress_p=0);


            /**
             * FilterGaussSep applies a 2D Gaussian to the source image
             * which will be replaced by the result. The pixels to be
             * filtered are those marked with a set bit in the compact 
             * binary f_bitMask. Since this function splits up multithreaded 
             * filtering into a horizontal and vertical 1D pass over the 
             * images, it will be faster than a naive 2D convolution.
             * However, a temporary image of equal characteristics as the 
             * source/destination image must be passed.
             *
             * @param f_srcDestImage the image to be filtered, also holds the
             * result afterwards
             * @param f_bitMaskImage the single-channel compact binary image 
             * describing the pixels to be filtered (set bits)
             * @param f_tempImage a temporary image of the same 
             * characteristics as the source/destination image holding
             * the intermediate output of horizontal filtering
             * @param f_filterSize the horizontal and vertical filter
             * size, will be adjusted to positive odd integers if necessary
             * @param f_progress_p the optional progress indicator
             */
            void filterGaussSep(CFloatImage& f_srcDestImage, const CU8Image& 
            f_bitMaskImage, image::CFloatImage& f_tempImage, const data::
            CPoint2DInt& f_filterSize, common::CProgress* const 
            f_progress_p=0);
    };
} // namespace image

#endif
