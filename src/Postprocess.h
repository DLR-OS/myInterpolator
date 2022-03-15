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
             * FilterGauss applies a 2D Gaussian to the source image
             * and writes the result to the destination image. The pixels
             * to be interpolated are those mask image pixel with a 
             * positive intensity. This function is implemented as a
             * multithreaded straightforward 2D convolution.
             *
             * @param f_srcImage the image to be filtered
             * @param f_maskImage a single-channel image describing the 
             * pixels to be filtered (positive intensities)
             * @param f_destImage the output image, must have the same
             * dimensions and channel count as the source image
             * @param f_filterSize the horizontal and vertical filter
             * size, will be adjusted to positive odd integers if necessary
             * @param f_progress_p the optional progress indicator
             */
            void filterGauss(const CFloatImage& f_srcImage, CFloatImage& 
            f_maskImage, CFloatImage& f_destImage, const data::CPoint2DInt& 
            f_filterSize, common::CProgress* const f_progress_p=0);


            /**
             * FilterGaussSep applies a 2D Gaussian to the source image
             * which will be replaced by the result. The pixels
             * to be interpolated are those mask image pixel with a 
             * positive intensity. This function splits up multithreaded 
             * filtering into a horizontal and vertical 1D pass over the 
             * images and hence will be faster than a naive 2D convolution.
             * However, a temporary image of the same dimensions and channel
             * count as the source/destination image must be passed.
             *
             * @param f_srcDestImage the image to be filtered, also holds the
             * result
             * @param f_maskImage a single-channel image describing the 
             * pixels to be filtered (positive intensities)
             * @param f_tempImage a temporary image of the same 
             * characteristics as the source/destination image holding the
             * intermediate output of horizontal filtering
             * @param f_filterSize the horizontal and vertical filter
             * size, will be adjusted to positive odd integers if necessary
             * @param f_progress_p the optional progress indicator
             */
            void filterGaussSep(CFloatImage& f_srcDestImage, CFloatImage& 
            f_maskImage, image::CFloatImage& f_tempImage, const data::
            CPoint2DInt& f_filterSize, common::CProgress* const 
            f_progress_p=0);
    };
} // namespace image

#endif
