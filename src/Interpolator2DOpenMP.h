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

#ifndef IMAGE_INTERPOLATOR_2D_OPENMP_H
#define IMAGE_INTERPOLATOR_2D_OPENMP_H

#ifdef _OPENMP
    #include <omp.h>
#endif

#include <vector>
#include "Image.h"
#include "Progress.h"

namespace image {


    /**
     * CInterpolator2DOpenMP interpolates images using gap-less 
     * Bresenham sweeps. Mutithreading is via OpenMP, however, code
     * is organized so an OpenCL version can be derived straighforwardly.
     *
     * @author Dirk 'jtk' Frommholz
     * @date Dec 28, 2020
     */
    class CInterpolator2DOpenMP {
        public:


            /**
             * The constructor.
             * 
             * @param f_inImage the input image to be interpolated
             * @param f_outImage the interpolated output image
             * @param f_weightImage the image containing the aggregated
             * weights per hole pixel (must have been initialized to the
             * minimum value supported by the pixel data type)
             * @param f_bgnd the color of the holes per input color channel 
             * @param f_dirCount the number of directions
             * @param f_angularOffset the angle in degrees to start the
             * sweeps with 
             * @param f_idwSmoothness the smoothness exponent
             * @param f_disableOC true to disable oversampling compensation
             * and use the plain IDW formula
             * @param f_progress_p the progress indicator
             */
            CInterpolator2DOpenMP( CFloatImage& f_inImage, CFloatImage&
            f_outImage, CFloatImage& f_weightImage, const std::vector<
            float>& f_bgnd, const data::pos_type f_dirCount, const float 
            f_idwSmoothness, const float f_angularOffset, const bool 
            f_disableOC, common::CProgress* f_progress_p=0);


            /**
             * Interpolate interpolates holes using the parameters passed
             * to the constructor.
             *        
             * @param f_ppFilterSize the kernel size for postprocessing
             * (negative coordinates=disable)
             * @param f_pointStats the number of foreground, copied 
             * background and  interpolated background pixels in 
             * .m_x, .m_y and .m_z respectively 
             * @param f_timingStats the runtime of initialization, interpolation,
             * normalization, and filtering respectively
             */
            void interpolate(const data::CPoint2DInt& f_ppFilterSize, 
            data::CPoint3DInt* f_pointStats_p=0, std::vector<double>* 
            f_timingStats_p=0);


            /**
             * The destructor.
             */
            ~CInterpolator2DOpenMP() { }


        private:


            /**
             * SweepMode defines how to traverse the image pixels in
             * order to interpolate the holes from neighboring foreground
             * pixels. One or two of these sweeps must be performed to touch 
             * each pixel at least once. The sweep mode to be chosen depends
             * on the path direction.
             */
            enum SweepMode {


                /**
                 * Don't execute the sweep (no-op). Default.
                 */
                SWEEP_NONE                  = 0,


                /**
                 * Interpolate using paths starting at the top image border.
                 */
                SWEEP_STARTS_TOP_ROW        = 1,


                /**
                 * Interpolate using paths starting at the bottom image 
                 * border.
                 */
                SWEEP_STARTS_BOTTOM_ROW     = 2,


                /**
                 * Interpolate using paths starting at the left image border.
                 */
                SWEEP_STARTS_LEFT_COL       = 4,


                /**
                 * Interpolate using paths starting at the right image 
                 * border.
                 */
                SWEEP_STARTS_RIGHT_COL      = 8
            };


            /**
             * The input image to be interpolated.
             */
            CFloatImage& m_inImage;


            /**
             * The interpolated output image.
             */
            CFloatImage& m_outImage;


            /**
             * The image storing the aggregated pixel weights, for 
             * normalization. 
             */
            CFloatImage& m_weightImage;


            /**
             * The smoothness to be used during inverse distance weighting,
             * i.e. the exponent of the distance. 
             */
            float m_idwSmoothness;


            /**
             * The background color.
             */
            std::vector<float> m_bgnd;


            /**
             * The number of path directions to be used for the interpolation.
             */
            data::pos_type m_dirCount;


            /**
             * The angle in degrees of the first path to be constructed.
             */
            float m_angularOffset;


            /**
             * True if oversampling compensation shall not be used but
             * the plain IDW formula.
             */
            bool m_disableOC;


        private:


            /**
             * The optional progress indicator.
             */
            common::CProgress* const m_progress_p;


            /**
             * PixelIsBgnd returns true if the input image pixel is
             * of the background color and false otherwise.                                                     
             *
             * @param f_pos the pixel position to be tested
             * @return true if the pixel at f_pos is background, false
             * if not
             */
            bool pixelIsBgnd(const data::CPoint2DInt& f_pos) const;


            /**
             * ComputeLineRuns precalculates the run sequence of a discrete
             * line using Bresenhams 4-connected component line tracing 
             * algorithm. The runs describe horizontal or vertical line
             * segments arbitrarily oriented paths are approximated with. 
             *
             * @param f_start the start position of the line
             * @param f_dir the direction vector of the line
             * @param f_length the line length
             * @param f_runs the resulting run sequence where each line run 
             * is described by its offset to the line start and the length
             * of the run
             */
            data::pos_type computeLineRuns(const data::CPoint2DInt& f_start, 
            const data::CPoint2DInt& f_dir, data::CPoint2DInt::CoordType 
            f_length, std::vector<data::CPoint2DInt>& f_runs);

            /**
             * SweepArbDirectionReplay performs inverse distance waiting
             * along a line run sequence.
             *
             * @param f_dirCount the direction count
             * @param f_spIndex the start point index
             * @param f_f_startPixels_p the image positions of the paths along
             * which IDW is to be performed along
             * @param f_prefAndNonPrefDirStep one discrete line steps in the
             * preferred (major) and non-preferred (minor) line direction
             * for the horizontal and vertical direction (must be 1 or -1)
             * @param  f_prefDirOnlyStep one discrete line step in the
             * preferred (major) line directio
             * @param f_distanceInc the distance increment from 
             * one line pixel to the next one on the continuous line
             * (i.e. the Euclidean distance between two neighboring
             * line pixels)
             * @param f_idwSmoothness the IDW smoothness (>0)
             * @param f_pathRuns_p the run sequence of the path of direction 
             * f_dir to be traced; must have enough runs to cover the image
             * (i.e. >max(w,h)) since there is no length checking here!
             * @param f_skipRuns 0=do not initially skip f_spIndex line runs 
             * (primary sweeps=, 1=initially skip f_spIndex line runs 
             * (secondary sweeps) before performing IDW along the paths
             * @param f_inImage_p the input image
             * @param f_weigthImage_p the image storing the IDWs
             * @param f_outImage_p the partially interpolated output
             * image
             * @param f_width the common width of the images
             * @param f_height the common height of the images
             * @param f_nrOfChannels the common channel count of the images
             */
            void sweepArbDirectionReplay(const uint64_t f_dirCount,
            const uint64_t f_spIndex,  const data::CPoint2DInt* 
            f_startPixels_p, const data::CPoint2DInt& f_prefAndNonPrefDirStep,
            const data::CPoint2DInt& f_prefDirOnlyStep, const float 
            f_distanceInc, const float f_idwSmoothness, const data::
            CPoint2DInt* f_pathRuns_p, const data::pos_type f_skipRuns,
            const float* f_inImage_p, float* f_weightImage_p, float* 
            f_outImage_p, const uint64_t f_width, const uint64_t 
            f_height, const uint64_t f_nrOfChannels);


            /**
             * SweepPrimary performs the primary sweep over the image for
             * to given path direction.
             *
             * @param f_sweepMode indicates at which image border the primary
             * sweep shall start (see sweep scheme figure) 
             * @param f_direction the direction of the set of parallel paths
             * in pixel space as a integer 2D vector           
             * @param f_prefAndNonPrefDirStep precomputed pixel increment/
             * decrement in the preferred/major and non-preferred (minor) 
             * line direction, i.e., for horizontal directions this will be 
             * (+/-1;0)
             * @param f_prefDirOnlyStep precomputed pixel increment/
             * decrement in the preferred/major line direction only
             * @param f_distanceInc the Euclidean distance increment between
             * two pixels on the paths of the current direction
             * @param f_idwSmoothness the exponent controlling the 
             * interpolation smoothness
             * @param f_pathRuns the line run sequence for the current 
             * direction, will be computed by this function to be passed to
             * the secondary sweep; the vector must have been already resized 
             * to more elements than the maximum of the image dimensions
             * @param f_forward true if the paths are followed in forward 
             * direction or not, this determines against which image border
             * we have to test whether this image frame is left
             */ 
            data::pos_type sweepPrimary(const SweepMode 
            f_sweepMode, const data::CPoint2DInt& f_direction, const data::
            CPoint2DInt& f_prefAndNonPrefDirStep,
            const data::CPoint2DInt& f_prefDirOnlyStep, const float 
            f_distanceInc, const float f_idwSmoothness, std::vector<data::
            CPoint2DInt>& f_pathRuns, std::vector<data::CPoint2DInt>& 
            f_startPixels, bool f_forward);


            /**
             * SweepSecondary performs the aecondary sweep over the image for
             * to given path direction.
             *
             * @param f_sweepMode indicates at which image border the primary
             * sweep shall start (see sweep scheme figure) 
             * @param f_direction the direction of the set of parallel paths
             * in pixel space as a integer 2D vector           
             * @param f_prefAndNonPrefDirStep precomputed pixel increment/
             * decrement in the preferred/major and non-preferred (minor) 
             * line direction, i.e., for horizontal directions this will be 
             * (+/-1;0)
             * @param f_prefDirOnlyStep precomputed pixel increment/
             * decrement in the preferred/major line direction only
             * @param f_distanceInc the Euclidean distance increment between
             * two pixels on the paths of the current direction
             * @param f_idwSmoothness the exponent controlling the 
             * interpolation smoothness
             * @param f_pathRuns the line run sequence for the current 
             * direction, from the primary sweep
             * @param f_forward true if the paths are followed in forward 
             * direction or not, this determines against which image border
             * we have to test whether this image frame is left
             */ 
            void sweepSecondary(const SweepMode f_sweepMode, const 
            data::CPoint2DInt& f_prefAndNonPrefDirStep, const data::
            CPoint2DInt& f_prefDirOnlyStep, const float f_distanceInc,
            const float f_idwSmoothness, const std::vector<data::
            CPoint2DInt>& f_pathRuns, const data::pos_type f_runCount,
            std::vector<data::CPoint2DInt>& f_startPixels, bool f_forward);

            
            /**
             * InitWeightImage initializes the weight bitmap which is used
             * to normalized the aggregated intensities after interpolation
             * to preserve the average image intensity. Valid (existing) 
             * pixels will bet set to a negative value.
             */
            void initWeightImage();


            /**
             * InterpolateOptimizedBresenham interpolates the image. The
             * function runs the primary and secondary sweeps for the
             * desired number of path directions.
             */
            void interpolateOptimizedBresenham();


            /**
             * Normalize normalizes the interpolated pixels by the 
             * aggregated weights preserving the average intensity.
             */
            data::CPoint3DInt normalize();
    };


} // namespace image

#endif // IMAGE_INTERPOLATOR_2D_H
