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

#include <vector>
#include "Postprocess.h"


/****************************************************************************/


void image::CPostprocess::makeBinaryBitMask(const CFloatImage& f_fpMask, 
CU8Image& f_bitMask) {

    // create bit mask image
    f_bitMask.create(common::max2<data::pos_type>(1, f_fpMask.getWidth()/8), 
    f_fpMask.getHeight(), 1, 0);

    // fill bits    
    for (data::pos_type y=0; y<f_fpMask.getHeight(); ++y) {
        for (data::pos_type x=0; x<f_fpMask.getWidth(); ++x) {

            const uint8_t pixelMask=*(f_fpMask.getIterator(x, y, 0))>0.0f ?
            (1<<(x%8)): 0;
            
            *(f_bitMask.getIterator(x/8, y, 0))|=pixelMask;
        }
    }
}


/****************************************************************************/


void image::CPostprocess::filterGauss(CFloatImage& f_srcDestImage, 
const CU8Image& f_bitMaskImage, CFloatImage& f_tempImage, const data::
CPoint2DInt& f_filterSize, common::CProgress* const f_progress_p) {

    // guarantee odd filters
    const data::CPoint2DInt oddFilterSize(
        f_filterSize.m_i % 2 ? f_filterSize.m_i : f_filterSize.m_i+1,
        f_filterSize.m_j % 2 ? f_filterSize.m_j : f_filterSize.m_j+1
    );

    // half filter size
    const data::CPoint2DInt halfKs(
        common::max2<data::CPoint2DInt::CoordType>(1, oddFilterSize.m_i/2),
        common::max2<data::CPoint2DInt::CoordType>(1, oddFilterSize.m_j/2)
    );

    #pragma omp parallel for schedule (dynamic)
    for (data::pos_type y=0; y<f_srcDestImage.getHeight(); ++y) {

        if (f_progress_p) {
            
            #pragma omp critical
            {
                const double percentage=static_cast<double>(y+1)/
                f_srcDestImage.getHeight();

                f_progress_p->setPercentage(percentage);
                f_progress_p->showProgress();
            }
        }

        for (data::pos_type x=0; x<f_srcDestImage.getWidth(); ++x) {    

            // pixel to be filtered, kernel center
            const data::CPoint2DInt centerPos(x, y);

            // filter interpolated pixels only from fully h-filtered
            // temporary image - implictly keep non-interpolated pixels
            const uint8_t mustFilter=*(f_bitMaskImage.getIterator(
            x/8, y, 0)) & (1 << (x%8));

            if (mustFilter) {

                data::CPoint2DInt kernelPos;
                float outSample=0;
                float weightSum=0;

                // init output pixel for accumulation
                float* tempIt=f_tempImage.getIterator(centerPos, 0);
                for (data::pos_type ch=0; ch<f_tempImage.getNrOfChannels(); 
                ++ch, ++tempIt) {
                    *tempIt=0;
                }
                tempIt-=f_tempImage.getNrOfChannels();

                for (data::pos_type ky=-halfKs.m_j; ky<=halfKs.m_j; ++ky) {
                    
                    kernelPos.m_j=common::period(ky+y, f_srcDestImage.
                    getHeight());

                    for (data::pos_type kx=-halfKs.m_i; kx<=halfKs.m_i; 
                    ++kx) {

                        kernelPos.m_i=common::period(kx+x, f_srcDestImage.
                        getWidth());

                        float* srcIt=f_srcDestImage.getIterator(kernelPos, 0);
                        const float weight=common::gaussianDistribution(kx, 0, 
                        oddFilterSize.m_i/3)*common::gaussianDistribution(ky, 0, 
                        oddFilterSize.m_j/3);
                        weightSum+=weight;                        
                        
                        for (data::pos_type ch=0; ch<f_srcDestImage.
                        getNrOfChannels(); ++ch, ++tempIt, ++srcIt) {
                            const float* srcIt=f_srcDestImage.getIterator(
                            kernelPos, ch);
                            *tempIt+=*srcIt*weight;
                        }   

                        srcIt-=f_srcDestImage.getNrOfChannels();                     
                        tempIt-=f_tempImage.getNrOfChannels();                     
                    } // kx
                } // ky

                // normalize                
                for (data::pos_type ch=0; ch<f_tempImage.
                getNrOfChannels(); ++ch, ++tempIt) {                            
                    *tempIt/=(fabs(weightSum)>0) ? weightSum : 1;
                }
                tempIt-=f_tempImage.getNrOfChannels();
            } else {

                // do not filter, plain copy
                for (data::pos_type ch=0; ch<f_srcDestImage.getNrOfChannels(); 
                ++ch) {
                    *(f_tempImage.getIterator(centerPos, ch))=
                    *(f_srcDestImage.getIterator(centerPos, ch));
                } 

            } // mustFilter
        } // x
    } // y

    // copy result back to input image
    f_srcDestImage=f_tempImage;
}


/****************************************************************************/

//
// Note that we must
//
// 1. filter the input entirely into the temp image in the horizontal pass
// 2. filter the temp image vertically keeping the non-interpolated points
//
// to get the equivalent of a genuine 2D convolution on the interpolated
// pixels. Due to a separate binary bitmask, filtering can take placed in 
// an interleaved fashion and does not have to be run band-wise.
//
void image::CPostprocess::filterGaussSep(image::CFloatImage& f_srcDestImage, 
const image::CU8Image& f_bitMaskImage, image::CFloatImage& f_tempImage, const
data::CPoint2DInt& f_filterSize, common::CProgress* const f_progress_p) {

    const data::CPoint2DInt oddFilterSize(
        f_filterSize.m_i % 2 ? f_filterSize.m_i : f_filterSize.m_i+1,
        f_filterSize.m_j % 2 ? f_filterSize.m_j : f_filterSize.m_j+1
    );
    const data::CPoint2DInt halfFilterSize=oddFilterSize/2;

    //
    // horizontal pass (srcDestImage -> tempImage)
    // filter everything (!), also the existing pixels to achieve the
    // equivalent of a naive 2D convolution
    //

    // precompute horizontal kernel and its normalization factor
    // assum filter size is 3*sigma (99% of values), so sigma=filterSize/3
    std::vector<float> hKernel(common::max2<data::CPoint2DInt::CoordType>(1, 
    oddFilterSize.m_i));
    float hNormFactor=0;

    for (auto kx=0; kx<hKernel.size(); ++kx) {
        hKernel[kx]=common::gaussianDistribution(kx-halfFilterSize.m_i, 0, 
        oddFilterSize.m_i/3);    
        hNormFactor+=hKernel[kx];
    }
    
    hNormFactor=1.0/hNormFactor;

    #pragma omp parallel for schedule (dynamic)
    for (data::pos_type y=0; y<f_srcDestImage.getHeight(); ++y) {    

        if (f_progress_p) {
            
            #pragma omp critical
            {
                // horizontal filter is only half of the cake
                const double percentage=0.5*(y+1)/f_srcDestImage.getHeight();

                if (f_progress_p->getPercentage()<percentage) {
                    f_progress_p->setPercentage(percentage);
                }

                f_progress_p->showProgress();
            }
        }

        for (data::pos_type x=0; x<f_srcDestImage.getWidth(); ++x) {    

            float* destIt=f_tempImage.getIterator(x, y, 0);

            // initialize target pixel for aggregation
            for (data::pos_type ch=0; ch<f_tempImage.getNrOfChannels(); 
            ++ch, ++destIt) {
                *destIt=0;
            }
            destIt-=f_tempImage.getNrOfChannels();

            // yes, apply low-pass kernel horizontally
            for (auto kx=-halfFilterSize.m_i; kx<=halfFilterSize.m_i; 
            ++kx) {
                float* srcIt=f_srcDestImage.getIterator(common::period(
                x+kx, f_srcDestImage.getWidth()), y, 0);

                for (data::pos_type ch=0; ch<f_srcDestImage.
                getNrOfChannels(); ++ch, ++srcIt, ++destIt) {
                    *destIt+=*srcIt*hKernel[kx+halfFilterSize.m_i]*
                    hNormFactor;
                }

                destIt-=f_tempImage.getNrOfChannels();
            }

        } // x
    } // y

    //
    // vertical pass (tempImage -> srcDestImage)
    // now filter selectively the pixels only that have a set bit inside
    // trhe binary bit mask
    //    

    // precompute vertical kernel and its normalization factor
    std::vector<float> vKernel(common::max2<data::CPoint2DInt::CoordType>(1, 
    oddFilterSize.m_j));
    float vNormFactor=0;

    for (auto ky=0; ky<vKernel.size(); ++ky) {
        vKernel[ky]=common::gaussianDistribution(ky-halfFilterSize.m_j, 0, 
        oddFilterSize.m_j/3);    
        vNormFactor+=vKernel[ky];
    }
    
    vNormFactor=1.0/vNormFactor;

    #pragma omp parallel for schedule (dynamic)
    for (data::pos_type x=0; x<f_srcDestImage.getWidth(); ++x) {    

        if (f_progress_p) {
            
            #pragma omp critical
            {
                // vertical filter is other half of the cake
                const double percentage=0.5+(0.5*(x+1)/f_srcDestImage.
                getWidth());

                if (f_progress_p->getPercentage()<percentage) {
                    f_progress_p->setPercentage(percentage);
                }
                f_progress_p->showProgress();
            }
        }

        for (data::pos_type y=0; y<f_srcDestImage.getHeight(); ++y) {    

            float* destIt=f_srcDestImage.getIterator(x, y, 0);

            // filter interpolated pixels only from fully h-filtered
            // temporary image - implictly keep non-interpolated pixels
            const uint8_t mustFilter=*(f_bitMaskImage.getIterator(
            x/8, y, 0)) & (1 << (x%8));

            if (mustFilter) {

                // initialize target pixel for aggregation
                for (data::pos_type ch=0; ch<f_tempImage.getNrOfChannels(); 
                ++ch, ++destIt) {
                    *destIt=0;
                }
                destIt-=f_srcDestImage.getNrOfChannels();

                // yes, apply low-pass kernel horizontally
                for (auto ky=-halfFilterSize.m_j; ky<=halfFilterSize.m_j; 
                ++ky) {
                    float* srcIt=f_tempImage.getIterator(x, common::period(
                    y+ky, f_srcDestImage.getHeight()), 0);

                    for (data::pos_type ch=0; ch<f_srcDestImage.
                    getNrOfChannels(); ++ch, ++srcIt, ++destIt) {
                        *destIt+=*srcIt*vKernel[ky+halfFilterSize.m_j]*
                        vNormFactor;
                    }

                    destIt-=f_srcDestImage.getNrOfChannels();
                }
            } // filter
        } // x
    } // y
}


/****************************************************************************/



