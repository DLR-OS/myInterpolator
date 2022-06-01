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

#include <chrono>

#include "coolmath.h"
#include "Interpolator2DOpenMP.h"
#include "Postprocess.h"


/****************************************************************************/


image::CInterpolator2DOpenMP::CInterpolator2DOpenMP(CFloatImage& 
f_inImage, CFloatImage& f_outImage, CFloatImage& f_weightImage, const 
std::vector<float>& f_bgnd, const data::pos_type f_dirCount, const float 
f_angularOffset, const float f_idwSmoothness, const bool f_disableOC, 
common::CProgress* const f_progress_p):
m_inImage(f_inImage), m_outImage(f_outImage), m_idwSmoothness(
f_idwSmoothness), m_weightImage(f_weightImage), m_dirCount(f_dirCount), 
m_angularOffset(f_angularOffset), m_bgnd(f_bgnd), m_disableOC(f_disableOC), 
m_progress_p(f_progress_p) {
}


/****************************************************************************/


data::pos_type image::CInterpolator2DOpenMP::computeLineRuns(const data::
CPoint2DInt& f_start, const data::CPoint2DInt& f_dir, data::CPoint2DInt::
CoordType f_length, std::vector<data::CPoint2DInt>& f_runs) {
       
    // init Bresenham variables
    data::CPoint2DInt linePos(f_start);  
 
    const data::pos_type dx=abs(f_dir.m_i);
    const data::pos_type sx=(f_dir.m_i<0) ? -1 : 1;
    const data::pos_type dy=-abs(f_dir.m_j);
    const data::pos_type sy=(f_dir.m_j<0) ? -1 : 1;
    data::pos_type bresErr=dx+dy;
    data::pos_type bresErr2;    

    // run pointer
    std::vector<data::CPoint2DInt>::iterator runIt=f_runs.begin();   
    data::CPoint2DInt::CoordType startOfCurrentRun=0;
    data::CPoint2DInt::CoordType nrOfRuns=0;   

    // path loop (Bresenham)
    while (f_length>0) {

        // prepare initial run
        if (nrOfRuns==0) {
            nrOfRuns=1;
            runIt->setCoords(startOfCurrentRun, 1);
        } else {

            // increment current run length else
            ++(runIt->m_j);
        }

        // update path length
        --f_length;
        
        bresErr2=2*bresErr;
        bool curRunIsFinished=false;

        if (bresErr2>dy) { 
            bresErr+=dy; 
            linePos.m_i+=sx; 

            // making a horizontal step on a vertically dominating line
            // means the current run is finished            
            curRunIsFinished=(-dy>=dx);
        } 

        if (bresErr2<dx) { 
            bresErr+=dx; 
            linePos.m_j+=sy;

            // making a vertical step on a horizontally dominating line
            // means the current run is finished
            curRunIsFinished=(curRunIsFinished==false) ? (dx>-dy) : true;
        } 

        // must start new run ?
        if (curRunIsFinished) {
            startOfCurrentRun+=runIt->m_j;

            ++runIt;
            runIt->m_i=startOfCurrentRun;
            runIt->m_j=0;
            ++nrOfRuns;
        }

    } // while

    // return number of runs
    return nrOfRuns;
} 


/****************************************************************************/


bool image::CInterpolator2DOpenMP::pixelIsBgnd(const data::CPoint2DInt& 
f_pos) const {

    const float* inIt=m_inImage.getIterator(f_pos);

    for (const auto& b: m_bgnd) {
        if (common::fp_unequal(b, *inIt++)) {
            return false;
        }
    }

    return true;
}


/****************************************************************************/


void image::CInterpolator2DOpenMP::sweepArbDirectionReplay(
    const uint64_t f_dirCount,
    const uint64_t f_spIndex, 
    const data::CPoint2DInt* f_startPixels_p,
    const data::CPoint2DInt& f_prefAndNonPrefDirStep, 
    const data::CPoint2DInt& f_prefDirOnlyStep,
    const float f_distanceInc, 
    const float f_idwSmoothness,
    const data::CPoint2DInt* f_pathRuns_p,
    const data::pos_type f_skipRuns,
    // imagery
    const float* f_inImage_p,
    float* f_weightImage_p,
    float* f_outImage_p,
    const uint64_t f_width,
    const uint64_t f_height,
    const uint64_t f_nrOfChannels
) {

    // actual start pixel index, must start at f_skipRuns=1 on 
    // secondary sweeps to avoid overlaps
    const uint64_t spIndex=f_spIndex+f_skipRuns;

    // the distance to the current pixel (avoid sqrt calculation)
    // negative if we havn't seen a fgnd pixel yet
    float distanceToCurPixel=-1;

    //
    // init Bresenham variables
    //

    // f_startPixels must be outside the image frame when
    // f_skipRuns=1 (secondary sweeps) for the sucessor of the maximum
    // valid run index from the primary sweep !!!                          
    data::CPoint2DInt linePos(f_startPixels_p[spIndex]);

    // current run
    // we skip the first spIndex runs (only if f_skipRuns is 1
    // for the secondary sweeps)
    // this may end up outside the range of valid runs from the
    // primary sweep by 1 on secondary sweeps, but since in this
    // case - as stated above - we'll be outside the image frame
    // no harm will be done
    const data::CPoint2DInt* runIt=f_pathRuns_p+f_skipRuns*spIndex;

    // the length of the current run
    data::pos_type runLength=runIt->m_j;

    // initialize image pointers
    const data::pos_type pixelOffset=(linePos.m_j*f_width+linePos.m_i)*
    f_nrOfChannels;

    // current input image position (channel #0)
    const float* inImagePtr_p=f_inImage_p+pixelOffset;
    // current weight image position (channel #0)
    float* weightImagePtr_p=f_weightImage_p+pixelOffset;
    // current output image position (channel #0)
    float* outImagePtr_p=f_outImage_p+pixelOffset;
    // last foreground pixel (channel #0)
    const float* lastFgndPtr_p=0;

    // switch to next pixel on major+minor step
    const int64_t prefAndNonPrefDirPixelUpdate=(
    f_prefAndNonPrefDirStep.m_j*f_width+f_prefAndNonPrefDirStep.
    m_i)*f_nrOfChannels;

    // switch to next pixel on major step only
    const int64_t prefOnlyPixelUpdate=(f_prefDirOnlyStep.m_j*f_width+
    f_prefDirOnlyStep.m_i)*f_nrOfChannels;

    // path loop (modified Bresenham)
    while ((linePos.m_j>=0) && (linePos.m_i>=0) && (linePos.m_i<
    f_width) && (linePos.m_j<f_height)) {

        // check if input pixel is fgnd or bgnd by looking at the first
        // color channel of the weight image (this comparison is more 
        // robust to numerical issues than ">=0" when the weight image
        // has been initialized to something very less then -1.0)
        if (*weightImagePtr_p>-1.0f) { 

            // bgnd ? check if the fgnd pixel pointer has been initialized
            // if not, we cannot interpolate this bgnd pixel and skip it
            if (lastFgndPtr_p) {
                
                // distance to this bgnd pixel
                distanceToCurPixel+=f_distanceInc;

                // compensate for oversampling close pixels, i.e.
                // weigh them lower by 1/oversampling aka 1/(dirCount/U)
                // =(U/dirCount)=(8*distance)/dirCount
                const float directionalWeight=(f_dirCount>0) ? 
                (8*distanceToCurPixel)/f_dirCount : 1;

                // yes, perform IDW
                // compute distance-dependent weight                
                const float weight=directionalWeight/pow(distanceToCurPixel, 
                f_idwSmoothness);

                // interpolate
                for (uint64_t ch=0; ch<f_nrOfChannels; ++ch, 
                ++inImagePtr_p, ++weightImagePtr_p, ++outImagePtr_p, 
                ++lastFgndPtr_p) {
                    *outImagePtr_p+=*lastFgndPtr_p*weight;
                    *weightImagePtr_p+=weight;
                }                

                // reset to channel #0
                inImagePtr_p-=f_nrOfChannels;
                weightImagePtr_p-=f_nrOfChannels;
                outImagePtr_p-=f_nrOfChannels;
                lastFgndPtr_p-=f_nrOfChannels;
            }
        } else {

            // copy input to output without changes
            for (uint64_t ch=0; ch<f_nrOfChannels; ++ch, ++inImagePtr_p, 
            ++outImagePtr_p) {
                *outImagePtr_p=*inImagePtr_p;
            }

            // reset to channel #0
            inImagePtr_p-=f_nrOfChannels;
            outImagePtr_p-=f_nrOfChannels;

            // fgnd -> remember as the last neighbor seen 
            lastFgndPtr_p=inImagePtr_p;
            distanceToCurPixel=0.0f;
        }

        // switch to next pixel of the current run
        // then check if current run got entirely processed
        if (--runLength<=0) {

            // make one step in both the non-preferred direction
            linePos+=f_prefAndNonPrefDirStep;

            // update image pointers
            inImagePtr_p+=prefAndNonPrefDirPixelUpdate;
            weightImagePtr_p+=prefAndNonPrefDirPixelUpdate;
            outImagePtr_p+=prefAndNonPrefDirPixelUpdate;

            // re-initialize run length from next run
            ++runIt;
            runLength=runIt->m_j;
        }
        else {

            // run not complete yet, just step in preferred direction
            linePos+=f_prefDirOnlyStep;

            // update image pointers
            inImagePtr_p+=prefOnlyPixelUpdate;
            weightImagePtr_p+=prefOnlyPixelUpdate;
            outImagePtr_p+=prefOnlyPixelUpdate;
        }
    } // whileInImageAndHaveRunsLeft
}


/****************************************************************************/


// use full Bresenham paths across the image in two non-overlapping sweeps
//
// OpenMP version, parallel per row/column
// Use dynamic scheduling since the length of the Bresenham paths will vary
// a lot

//
// PRIMARY SWEEP
//
data::pos_type image::CInterpolator2DOpenMP::sweepPrimary(const SweepMode 
f_sweepMode, const data::CPoint2DInt& f_direction, const data::CPoint2DInt& 
f_prefAndNonPrefDirStep, const data::CPoint2DInt& f_prefDirOnlyStep, 
const float f_distanceInc, const float f_idwSmoothness, std::vector<data::
CPoint2DInt>& f_pathRuns, std::vector<data::CPoint2DInt>& f_startPixels, 
bool f_forward) {

    // start path cost aggregation from top/bottom border ?
    if ((f_sweepMode & SWEEP_STARTS_TOP_ROW) || (f_sweepMode & 
    SWEEP_STARTS_BOTTOM_ROW)) {
        
        if (f_forward) {

            // first path of primary sweep is the longest, it always 
            // starts in an image corner, increment is +/-1
            const data::CPoint2DInt sweepStart(
                0, 
                (f_sweepMode & SWEEP_STARTS_TOP_ROW) ? 0 : 
                m_inImage.getHeight()-1
            );

            // we must continue tracing the Bresenham line until we
            // hit the opposite image border or we'll miss some pixels
            // during the secondary trace (i.e. in the 45° case, on
            // rectangular images)
            //
            // thus, since Bresenham always makes one step in the major
            // direction, the path length is the image height (since if
            // the image is wider than high, in case of the worst-case
            // path angle of 45 degrees, we'll hit the top/bottom row
            // before the left/right column, and otherwise, if the 
            // image is higher than wide we'll hit the l/r col first
            // and then the top/bottom, so we must take the height
            // as well)
            const data::pos_type runCount=computeLineRuns(sweepStart, 
            f_direction, m_inImage.getHeight(), f_pathRuns);

            // precompute Bresenham start positions for each run
            #pragma omp parallel for
            for (data::CPoint2DInt::CoordType x=0; x<m_inImage.
            getWidth(); ++x) { 
                f_startPixels[x].setCoords(x, sweepStart.m_j);
            }

            #pragma omp parallel for schedule(dynamic)
            for (data::CPoint2DInt::CoordType x=sweepStart.m_i; x<m_inImage.
            getWidth(); ++x) { 

                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, 
                x, f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 0, m_inImage.getIterator(), 
                m_weightImage.getIterator(), m_outImage.getIterator(), 
                m_inImage.getWidth(), m_inImage.getHeight(), 
                m_inImage.getNrOfChannels());
            }

            return runCount;

        } else {

            // same properties but reverse
            const data::CPoint2DInt sweepStart(
                m_inImage.getWidth()-1, 
                (f_sweepMode & SWEEP_STARTS_TOP_ROW) ? 0 : 
                m_inImage.getHeight()-1
            );

            // path length will always be the image height in reverse mode
            // as well
            const data::pos_type runCount=computeLineRuns(sweepStart, 
            f_direction, m_inImage.getHeight(), f_pathRuns);

            // precompute Bresenham start positions for each run
            #pragma omp parallel for
            for (data::CPoint2DInt::CoordType x=sweepStart.m_i; x>=0; --x) {
                f_startPixels[x].setCoords(x, sweepStart.m_j);
            }

            // parallel IDW calculation
            // doesn't matter if we go from sweepStart->0 or reverse
            #pragma omp parallel for schedule(dynamic)            
            for (data::CPoint2DInt::CoordType x=0; x<=sweepStart.m_i; ++x) {
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, x, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness,  
                f_pathRuns.data(), 0, m_inImage.getIterator(), 
                m_weightImage.getIterator(), m_outImage.getIterator(),
                m_inImage.getWidth(), m_inImage.getHeight(), 
                m_inImage.getNrOfChannels());                
            }

            return runCount;
        }
    } // top/bottom border

    // start path cost aggregation from left/right border ?
    if ((f_sweepMode & SWEEP_STARTS_LEFT_COL) || (f_sweepMode & 
    SWEEP_STARTS_RIGHT_COL)) {
        if (f_forward) {

            // first path of primary sweep is the longest, it always 
            // starts in an image corner, increment is +/-1
            const data::CPoint2DInt sweepStart(
                (f_sweepMode & SWEEP_STARTS_LEFT_COL) ? 0 : m_inImage.
                getWidth()-1,
                0
            );

            // path length will always be the image width
            const data::pos_type runCount=computeLineRuns(sweepStart, 
            f_direction, m_inImage.getWidth(), f_pathRuns);

            // precompute Bresenham start positions for each run
            #pragma omp parallel for
            for (data::CPoint2DInt::CoordType y=sweepStart.m_j; y<m_inImage.
            getHeight(); ++y) {
                f_startPixels[y].setCoords(sweepStart.m_i, y);
            }

            // forward sweep direction
            #pragma omp parallel for schedule(dynamic)
            for (data::CPoint2DInt::CoordType y=sweepStart.m_j; y<m_inImage.
            getHeight(); ++y) {
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, y, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 0, m_inImage.getIterator(), m_weightImage.
                getIterator(), m_outImage.getIterator(), m_inImage.
                getWidth(), m_inImage.getHeight(), m_inImage.
                getNrOfChannels());
            }

            return runCount;

        } else {

            // reverse sweep direction
            const data::CPoint2DInt sweepStart(
                (f_sweepMode & SWEEP_STARTS_LEFT_COL) ? 0 : m_inImage.
                getWidth()-1,
                m_inImage.getHeight()-1
            );

            // path length will always be the image width in reverse
            const data::pos_type runCount=computeLineRuns(sweepStart, 
            f_direction, m_inImage.getWidth(), f_pathRuns);

            // precompute Bresenham start positions for each run
            #pragma omp parallel for
            for (data::CPoint2DInt::CoordType y=sweepStart.m_j; y>=0; --y) {
                f_startPixels[y].setCoords(sweepStart.m_i, y);
            }

            // will again do 0->sweepStart instead of sweepStart->0
            #pragma omp parallel for schedule(dynamic)
            for (data::CPoint2DInt::CoordType y=0; y<=sweepStart.m_j; ++y) {
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, y, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 0, m_inImage.getIterator(), m_weightImage.
                getIterator(), m_outImage.getIterator(), m_inImage.
                getWidth(), m_inImage.getHeight(), m_inImage.
                getNrOfChannels());
            }

            return runCount;
        } 
    } // left/right border

    // invalid sweep mode
    if (f_sweepMode>SWEEP_STARTS_RIGHT_COL) {
        throw common::CException("Invalid primary sweep mode");
    }

    throw common::CException("No primary sweep executed");
}


/****************************************************************************/


void image::CInterpolator2DOpenMP::sweepSecondary(const SweepMode 
f_sweepMode, const data::CPoint2DInt& f_prefAndNonPrefDirStep, const data::
CPoint2DInt& f_prefDirOnlyStep, const float f_distanceInc, const float 
f_idwSmoothness, const std::vector<data::CPoint2DInt>& f_pathRuns, const 
data::pos_type f_runCount, std::vector<data::CPoint2DInt>& f_startPixels, 
bool f_forward) {

    // set the start position at the runCount index to something 
    // outside the image as a stop mark or the run pointer inside
    // sweepArbDirectionReplay() may become invalid            
    f_startPixels[f_runCount].setCoords(-1, -1);

    // start path cost aggregation from top/bottom border ?
    if ((f_sweepMode & SWEEP_STARTS_TOP_ROW) || (f_sweepMode & 
    SWEEP_STARTS_BOTTOM_ROW)) {

        if (f_forward) {            

            // precompute start pixels
            // start at run #0 since we offset the actual start run for the
            // secondary sweep inside sweepArbDirectionReplay()
            #pragma omp parallel for
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                f_startPixels[runNo].setCoords(f_pathRuns[runNo].m_i, 
                (f_sweepMode & SWEEP_STARTS_TOP_ROW) ? 0 : m_inImage.
                getHeight()-1);
            }

            // start with zero-th line run since the actual start run
            // will be computed inside sweepArbDirectionReplay() to skip
            // those pixels touched by the primary sweep, enable run skipping
            #pragma omp parallel for schedule(dynamic)
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, runNo, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 1, m_inImage.getIterator(), m_weightImage.
                getIterator(), m_outImage.getIterator(), m_inImage.
                getWidth(), m_inImage.getHeight(), m_inImage.
                getNrOfChannels());
            }

        } else {

            // reverse order, r->l

            // precompute start pixels
            #pragma omp parallel for
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 

                f_startPixels[runNo].setCoords(m_inImage.getWidth()-1-
                f_pathRuns[runNo].m_i, (f_sweepMode & 
                SWEEP_STARTS_TOP_ROW) ? 0 : m_inImage.getHeight()-1);
            }

            #pragma omp parallel for schedule(dynamic)
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, runNo, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 1, m_inImage.getIterator(), m_weightImage.
                getIterator(), m_outImage.getIterator(), m_inImage.
                getWidth(), m_inImage.getHeight(), m_inImage.
                getNrOfChannels());
            }
        }
    } // top/bottom border


    // start path cost aggregation from left/right border ?
    if ((f_sweepMode & SWEEP_STARTS_LEFT_COL) || (f_sweepMode & 
    SWEEP_STARTS_RIGHT_COL)) {
        if (f_forward) {

            // precompute start pixels
            #pragma omp parallel for
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 

                f_startPixels[runNo].setCoords(f_sweepMode & 
                SWEEP_STARTS_LEFT_COL ? 0 : m_inImage.
                getWidth()-1, f_pathRuns[runNo].m_i);
            }

            #pragma omp parallel for schedule(dynamic)
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, runNo, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, 
                f_pathRuns.data(), 1, m_inImage.getIterator(), 
                m_weightImage.getIterator(), m_outImage.getIterator(), 
                m_inImage.getWidth(), m_inImage.getHeight(), m_inImage.
                getNrOfChannels());      
            }

        } else {

            // reverse sweep direction

            // precompute start pixels
            #pragma omp parallel for
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                f_startPixels[runNo].setCoords(f_sweepMode & 
                SWEEP_STARTS_LEFT_COL ? 0 : m_inImage.
                getWidth()-1, m_inImage.getHeight()-1-f_pathRuns[
                runNo].m_i);
            }

            #pragma omp parallel for schedule(dynamic)
            for (data::pos_type runNo=0; runNo<f_runCount; ++runNo) { 
                sweepArbDirectionReplay(m_disableOC ? 0 : m_dirCount, runNo, 
                f_startPixels.data(), f_prefAndNonPrefDirStep, 
                f_prefDirOnlyStep, f_distanceInc, f_idwSmoothness, f_pathRuns.
                data(), 1, m_inImage.getIterator(), m_weightImage.getIterator(), 
                m_outImage.getIterator(), m_inImage.getWidth(), m_inImage.
                getHeight(), m_inImage.getNrOfChannels());      
            }
        } 
    } // left/right border

    // invalid sweep mode
    if (f_sweepMode>SWEEP_STARTS_RIGHT_COL) {
        throw common::CException("Invalid secondary sweep mode");
    }
}
 

/****************************************************************************/


void image::CInterpolator2DOpenMP::interpolateOptimizedBresenham() {        

    // an array storing the Bresenham runs of a path along which will be
    // interpolated
    // maximally spans the entire image, however, we add 1 to eliminate 
    // testing whether we have run out of runs
    std::vector<data::CPoint2DInt> pathRuns(common::max2(m_inImage.
    getWidth(), m_inImage.getHeight())+1);
    
    // the start positions of each Bresenham trace inside the input image
    // just preallocated, will be filled before primary/secondary sweeps
    std::vector<data::CPoint2DInt> startPixels(pathRuns.size());

    // show some progress
    if (m_progress_p) {
        m_progress_p->init("Image interpolation");
    }

    // process all directions
    for (data::pos_type dirNo=0; dirNo<m_dirCount; ++dirNo) {

        if (m_progress_p) {
            m_progress_p->setPercentage(static_cast<double>(dirNo)/
            m_dirCount);
            m_progress_p->showProgress();
        }

        // compute direction in fp, round to pixels
        const double alpha=fmod(dirNo*(360.0/
        m_dirCount)+m_angularOffset, 360.0);

        const data::CPoint2DInt pathDir(
            static_cast<data::CPoint2DInt::CoordType>(common::roundNat(
            cos(common::deg2rad(alpha))*2*common::max2(m_inImage.getWidth(),
            m_inImage.getHeight()))),
            static_cast<data::CPoint2DInt::CoordType>(common::roundNat(
            sin(common::deg2rad(alpha))*2*common::max2(m_inImage.getWidth(),
            m_inImage.getHeight())))
        );

        //
        // compute the discrete single line steps for this direction
        //

        // preferred direction is horizontal ?      
        const bool prefDirIsX=(abs(pathDir.m_i)>=abs(pathDir.m_j));

        // h/v steps to make in case we walk in the non-preferred
        // direction; this also includes a step in the preferred
        // direction in 4cc Bresenham
        const data::CPoint2DInt prefAndNonPrefDirStep(
            (pathDir.m_i>0) ? 1 : -1,
            (pathDir.m_j>0) ? 1 : -1
        );

        const data::CPoint2DInt prefDirOnlyStep(
            prefDirIsX ? prefAndNonPrefDirStep.m_i : 0,
            prefDirIsX ? 0 : prefAndNonPrefDirStep.m_j
        );

        // the length increment per pixel, for simple weight calculation
        // will be >=1 during the discrete trace
        const float distanceIncrement=sqrt(pathDir.m_i*pathDir.m_i+
        pathDir.m_j*pathDir.m_j)/common::max2(abs(pathDir.m_i), 
        abs(pathDir.m_j)); 
    
        //
        // traverse paths for all directions
        //

        //
        // do primary/secondary sweep according to path direction
        //
        // We must make sure all pixels are touched once per path direction;
        // this is guaranteed if the following things are obeyed to:
        //
        // 1. Performing a primary sweep for each per l/r border column 
        // or top/down border row depending on the major path direction
        // 
        // 2. Performing a secondary sweep offset by the length of the
        // k-th Bresenham run starting with k=1 (pixels for k=0 touched 
        // by primary sweep)
        // 
        // We implement this by first tracing a "prototype" line during
        // the primary sweep where tracing stops when we hit the opposite
        // image border (i.e. start trace at left -> prototype line ends
        // at right border) even if we leave the image.
        // The prototype line is stored as a sequence of Bresenham runs
        // which are then passed to the secondary sweep. The number
        // of the contained runs then exactly fits to touch the pixels
        // not touched during the primary run, and no more than that.
        //
        
        // evaluate vertical direction
        switch (common::sgn(pathDir.m_j)) {
            // purely horizontal paths, one sweep only needed
            case 0:

                if (pathDir.m_i>0) {
                    // left-to-right, cols up->down
                    sweepPrimary(SWEEP_STARTS_LEFT_COL, pathDir, 
                    prefAndNonPrefDirStep, prefDirOnlyStep, 
                    distanceIncrement, m_idwSmoothness, pathRuns, 
                    startPixels, true);
                } else {
                    // right-to-left, cols up->down
                    sweepPrimary(SWEEP_STARTS_RIGHT_COL, pathDir, 
                    prefAndNonPrefDirStep, prefDirOnlyStep, 
                    distanceIncrement, m_idwSmoothness, pathRuns, 
                    startPixels, true);
                }
            break; // purely horizontal sweeps

            // positive vertical slope, i.e. downward paths in image space
            case 1:
                switch (common::sgn(pathDir.m_i)) {
                    case 0:
                        // purely vertical, 1 sweep only, cols up->down
                        sweepPrimary(SWEEP_STARTS_TOP_ROW, pathDir, 
                        prefAndNonPrefDirStep, prefDirOnlyStep, 
                        distanceIncrement, m_idwSmoothness, pathRuns, 
                        startPixels, true);
                    break;
                    case 1:      
                        if (pathDir.m_i>=pathDir.m_j) {
                            // flat positive slope in x
                            // sweep #1 - start from left col, rows up->down
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_LEFT_COL, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, true);

                            // sweep #2 - start from top row, cols l->r,
                            // offset/increment by initial Bresenham run from sweep #1
                            sweepSecondary(SWEEP_STARTS_TOP_ROW, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns,
                            runCount, startPixels, true);
                        } else {
                            // steep positive slope in x            
                            // sweep #1 - start from top row, cols l->r
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_TOP_ROW, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, true);

                            // sweep #2 - start from left col, rows up->down
                            // offset/increment by initial Bresenham run from sweep #1
                            sweepSecondary(SWEEP_STARTS_LEFT_COL, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            runCount, startPixels, true);
                        }
                    break;
                    case -1:                        
                        if (-pathDir.m_i>=pathDir.m_j) {
                            // flat negative slope in x          
                            // sweep #1 - start from right col, rows up->down
                            // (cols up->down will guarantuee a non-truncated
                            // initial Bresenham run for sweep #2)
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_RIGHT_COL, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, true);

                            // sweep #2 - start from top row, cols right->left
                            // (right->left will guarantee gapless continuation
                            // of traced paths) 
                            sweepSecondary(SWEEP_STARTS_TOP_ROW, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns,
                            runCount, startPixels, false);
                        } else {
                            // steep negative slope in x
                            // sweep #1 - start from top row, cols right->left
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_TOP_ROW, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, false);

                            // sweep #2 - start from right col, rows up->down
                            sweepSecondary(SWEEP_STARTS_RIGHT_COL, 
                            prefAndNonPrefDirStep, prefDirOnlyStep,     
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            runCount, startPixels, true);
                        }
                    break;
                }
            break; // positive vertical slope

            // negative vertical slope in image space, i.e. upward paths
            // in image space
            case -1:
                switch (common::sgn(pathDir.m_i)) {
                    case 0:
                        // purely vertical, 1 sweep only, upwards, cols l->r
                        sweepPrimary(SWEEP_STARTS_BOTTOM_ROW, pathDir, 
                        prefAndNonPrefDirStep, prefDirOnlyStep, 
                        distanceIncrement, m_idwSmoothness, pathRuns, 
                        startPixels, true);
                    break;
                    case 1:      
                        if (pathDir.m_i>=-pathDir.m_j) {
                            // flat positive slope in x 
                            // sweep #1 - start from left col, rows down->up
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_LEFT_COL, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, false);

                            // sweep #2 - start from bottom row, cols left->right
                            sweepSecondary(SWEEP_STARTS_BOTTOM_ROW, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            runCount, startPixels, true);
                        } else {
                            // steep positive slope in x 
                            // sweep #1 - start from bottom row, cols left->right
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_BOTTOM_ROW, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, true);

                            // sweep #2 - start from left col, rows down->up
                            sweepSecondary(SWEEP_STARTS_LEFT_COL, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns,
                            runCount, startPixels, false);
                        }
                    break;
                    case -1:                              
                        if (-pathDir.m_i>=-pathDir.m_j) {
                            // flat negative slope in x 
                            // sweep #1 - start from right col, rows down->up
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_RIGHT_COL, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, false);

                            // sweep #2 - start from bottom row, cols right->left
                            sweepSecondary(SWEEP_STARTS_BOTTOM_ROW, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            runCount, startPixels, false);
                        } else {
                            // step negative slopes in x
                            // sweep #1 - start from bottom row, cols right->left
                            const data::pos_type runCount=
                            sweepPrimary(SWEEP_STARTS_BOTTOM_ROW, pathDir, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            startPixels, false);

                            // sweep #2 - start from right col, rows down->up                            
                            sweepSecondary(SWEEP_STARTS_RIGHT_COL, 
                            prefAndNonPrefDirStep, prefDirOnlyStep, 
                            distanceIncrement, m_idwSmoothness, pathRuns, 
                            runCount, startPixels, false);
                        }
                    break;
                }
            break; // negative vertical slope
        } // vertical slope

    } // pathNo
} 


/****************************************************************************/


data::CPoint3DInt image::CInterpolator2DOpenMP::normalize() {

    // normalize
    if (m_progress_p) {
        m_progress_p->init("Normalization");
    }

    // statistics
    data::CPoint3DInt globalPointStats;

    #pragma omp parallel for schedule(dynamic)
    for (data::pos_type y=0; y<m_weightImage.getHeight(); ++y) {

        // initialize per-thread iterator set 
        const float* inIt=m_inImage.getIterator(0, y, 0);
        const float* weightIt=m_weightImage.getIterator(0, y, 0);    
        float* outIt=m_outImage.getIterator(0, y, 0);
        data::CPoint3DInt  localPointStats;

        if (m_progress_p) {
            
            #pragma omp critical
            {
                const double percentage=static_cast<double>(y+1)/
                m_weightImage.getHeight();

                m_progress_p->setPercentage(percentage);
                m_progress_p->showProgress();
            }
        }

        for (data::pos_type x=0; x<m_weightImage.getWidth(); ++x) {
            for (data::pos_type ch=0; ch<m_weightImage.getNrOfChannels(); 
            ++ch, ++inIt, ++outIt, ++weightIt) {

                // we don't need to normalize non-bgnd pixels, and for
                // the fgnd/bgnd check, we abuse the weight image
                if (*weightIt>-1) {

                    // interpolated bgnd pixel ?
                    if (*weightIt==0.0f) {
                        // no, plain copy
                        *outIt=m_bgnd[ch];    
                        ++localPointStats.m_y;
                    } else {
                        // yes, normalize
                        *outIt/=(*weightIt); 
                        ++localPointStats.m_z;
                    }
                    
                } else {
                    // fgnd pixel, copy
                    *outIt=*inIt;
                    ++localPointStats.m_x;
                }
            }    
        }

        // aggregate thread-local stats into global stats
        #pragma omp critical
        globalPointStats+=localPointStats;
    }

    return globalPointStats/m_weightImage.getNrOfChannels();
}


/****************************************************************************/


void image::CInterpolator2DOpenMP::initWeightImage() {
    
    #pragma omp parallel for schedule(dynamic)
    for (data::pos_type y=0; y<m_weightImage.getHeight(); ++y) {
         for (data::pos_type x=0; x<m_weightImage.getWidth(); ++x) {         

            const data::CPoint2DInt pos(x, y);
            // must be < -1 for fgnd pixels and zero for bgnd
            // to enable some optimizations in tracing code
            const CFloatImage::SampleType initialWeight=pixelIsBgnd(pos) ? 
            0 : -10;

            CFloatImage::Iterator weightIt=m_weightImage.getIterator(pos, 0);
            for (data::pos_type ch=0; ch<m_weightImage.getNrOfChannels(); 
            ++ch, ++weightIt) {                
                *weightIt=initialWeight;
            } // ch
        } // x
    } // y
}


/****************************************************************************/


void image::CInterpolator2DOpenMP::interpolate(const data::CPoint2DInt& 
f_ppFilterSize, data::CPoint3DInt* f_pointStats_p, std::vector<double>* 
f_timingStats_p) {

    // init weight image to accelerate bgnd check
    auto startTime=std::chrono::system_clock::now();    
    initWeightImage();
    if (f_timingStats_p) {
        f_timingStats_p->clear();

        f_timingStats_p->push_back(std::chrono::duration_cast<std::chrono::
        milliseconds>(std::chrono::system_clock::now()-startTime).count());
    }

    // interpolate along paths
    startTime=std::chrono::system_clock::now();    
    interpolateOptimizedBresenham();  
    if (f_timingStats_p) {
        f_timingStats_p->push_back(std::chrono::duration_cast<std::chrono::
        milliseconds>(std::chrono::system_clock::now()-startTime).count());
    }

    // intensity normalization
    startTime=std::chrono::system_clock::now();    
    const data::CPoint3DInt pointStats=normalize();  

    if (f_pointStats_p) {
        *f_pointStats_p=pointStats;
    }

    if (f_timingStats_p) {
        f_timingStats_p->push_back(std::chrono::duration_cast<std::chrono::
        milliseconds>(std::chrono::system_clock::now()-startTime).count());
    }
    
    CPostprocess pp;

    // postprocessing (classic 2D convolution)   
    //if (f_ppFilterSize.m_i>0) {

    //    m_progress_p->init("Postprocessing");
    //    startTime=std::chrono::system_clock::now();    
    //    pp.filterGauss(m_outImage, m_weightImage, m_inImage, f_ppFilterSize, 
    //    m_progress_p);
    //    m_outImage=m_inImage;

    //    if (f_timingStats_p) {
    //        f_timingStats_p->push_back(std::chrono::duration_cast<std::
    //        chrono::milliseconds>(std::chrono::system_clock::now()-startTime).
    //        count());
    //    }
    //}

    // postprocessing (separable filter)   
    if (f_ppFilterSize.m_i>0) {
        m_progress_p->init("Postprocessing (separable filter)");
        startTime=std::chrono::system_clock::now();    
        pp.filterGaussSep(m_outImage, m_weightImage, m_inImage, f_ppFilterSize, 
        m_progress_p);

        if (f_timingStats_p) {
            f_timingStats_p->push_back(std::chrono::duration_cast<std::
            chrono::milliseconds>(std::chrono::system_clock::now()-startTime).
            count());
        }
    }
}


/****************************************************************************/



