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

#ifndef IMAGE_IP_PARAMETERS_H
#define IMAGE_IP_PARAMETERS_H


/**
 * CIpParameters stores the interpolation parameters as provided 
 * on the command line.
 */
struct CIpParameters {


    /**
     * The URI of the input image to be interpolated.
     */
    std::string m_inputImageUri;


    /**
     * The URI of the interpolated output image.
     */
    std::string m_outputImageUri;


    /**
     * The smoothness to be used during inverse distance weighting,
     * i.e. the exponent of the distance. 
     */
    double m_idwSmoothness;


    /**
     * The background color identifying holes to the interpolated.
     */
    std::vector<float> m_bgndCol;


    /**
     * The number of directions for path-wise interpolation.
     */
    data::pos_type m_dirCount;


    /**
     * The start angle of the first interpolation direction in degrees.
     */
    double m_angularOffset;


    /**
     * The filter size for postprocessing (negative=disabled).
     */
    data::CPoint2DInt m_ppFilterSize;


    /**
     * The number of OpenMP processing threads.
     */
    int m_threads;


    /**
     * True if contour oversampling compensation is disabled, i.e, when
     * true the plain IDW formula will be applied which leads to 
     * oversampling of close-by pixels by multiple paths.
     */
    bool m_disableOC;


    /**
     * The initializing constructor.
     */
    CIpParameters(): m_threads(-1), m_angularOffset(0), m_dirCount(128),
    m_idwSmoothness(2.0), m_disableOC(false), m_ppFilterSize(-1, -1) 
    { }


    /**
     * Print prints the interpolation parameters to the screen.
     */
    void print() {
        std::cout << "Interpolation parameters used" << std::endl;
        std::cout << "\tInput image URI        : " << m_inputImageUri
        << std::endl;

        std::string bgndColStr=" ";
        for (auto& col: m_bgndCol) {
            bgndColStr+=common::CConverter::xToStr<float>(col)+' ';
        }
        std::cout << "\tBackground color       :  [" << bgndColStr << "]"
        << std::endl;

        std::cout << "\tOutput image URI       : " << m_outputImageUri
        << std::endl;
        std::cout << "\tSmoothness exponent    : " << m_idwSmoothness 
        << std::endl;
        std::cout << "\tPath direction count   : " << m_dirCount
        << std::endl;
        std::cout << "\tPath angular offset    : " << m_angularOffset << 
        " degrees" << std::endl;
        std::cout << "\tCompensate oversampling: " << (m_disableOC ? 
        "no" : "yes") << std::endl;
        std::cout << "\tPostprocessing         : " << (m_ppFilterSize.m_i<0 ?
        "(disabled)" : common::CConverter::xToStr(m_ppFilterSize.m_i)+" x "+
        common::CConverter::xToStr(m_ppFilterSize.m_j)) << std::endl;

        
        #ifdef _OPENMP
        std::cout << "\tOpenMP threads         : " << (m_threads<0 ? 
        "auto" : common::CConverter::xToStr(m_threads)) << std::endl;
        #else
        std::cout << "\tThreads                      : no MT, 1" << std::endl;
        #endif

        std::cout << std::endl;
    }


    /**
     * The destructor.
     */
    ~CIpParameters() { }
};


#endif
