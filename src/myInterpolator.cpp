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
 * myInterpolator - sample program performing hole interpolation
 *
 * @author Dirk "jtk" Frommholz
 * @date December 23, 2020
 */

#include <iostream>
#include "coolmath.h"
#include "StdProgress.h"
#include "Image.h"
#include "Interpolator2DOpenMP.h"
#include "IpParameters.h"

#ifdef _OPENMP
    #include <omp.h>
#endif


/****************************************************************************/


/**
 * PrintUsage prints the usage information.
 *
 * @param f_programName the program name as of the command line
 * @param f_headerOnly print the software information header only not the 
 * full help screen
 */
void printUsage(const char* f_programName_p, bool f_headerOnly) {

    std::cout << f_programName_p << " - My interpolation utility" << 
    std::endl;
    std::cout << "Created by Dirk 'jtk' Frommholz, DLR OS-SEC Berlin-"
    "Adlershof" <<  std::endl;
    std::cout << "This is free software with absolutely no warranty, you "
    "may redistribute it under the conditions of the GPL v3 or later" << 
    std::endl;
    std::cout << std::endl;

    if (f_headerOnly) {
        return;
    }

    // usage
    std::cout << "usage: " << f_programName_p << " [options] <inImageUri> "
    "<outImageUri> " << std::endl << std::endl;

    std::cout << "options (in alphabetical order):" << std::endl << 
    std::endl;

    std::cout << "--bgnd <n> <int0> ... <intn-1>" << std::endl;
    std::cout << "\t\t" << "sets the background (hole) color "
    "for n components (defaults to 0 if not provided)" <<  std::endl << 
    std::endl;

    std::cout << "--dir <n>" << std::endl;
    std::cout << "\t\t" << "sets the number of interpolation directions" << 
    std::endl << std::endl;

    std::cout << "--disableOC" << std::endl;
    std::cout << "\t\t" << "disables oversamplig compensation, i.e., "
    "do not weigh close-by neighbor pixels less and distant pixels "
    "more in order to to simulate regular contour sampling" << 
    std::endl << std::endl;

    std::cout << "--idwSmoothness <s>" << std::endl;
    std::cout << "\t\t" << "sets the exponent to be applied to the "
    "distance to neighboring pixels in inverse distance weighting" <<
    "(1=linear, 2=quadratic decay)" << std::endl << std::endl;

    std::cout << "--offset <deg>" << std::endl;
    std::cout << "\t\t" << "sets the angular offset of the interpolation "
    "paths in degrees" << std::endl << std::endl;

    std::cout << "--postprocess sx sy" << std::endl;
    std::cout << "\t\t" << "smoothen interpolated pixels with Gauss filter " <<
    "of kernel size sx by sy " << std::endl << std::endl;
 
    std::cout << "--threads <count>" << std::endl;
    std::cout << "\t\t" << "sets the number of processing threads" << 
    std::endl << std::endl;
}


/****************************************************************************/


/**
 * ParseParameters parses the command line arguments and sets up a
 * CIpParameters record accordingly.
 * 
 * @param f_argc the argument count as of main()
 * @param f_argv_p the command line arguments as of main()
 * @param f_ipParameters the interpolation parameter object to be set up 
 * from the command line options
 */
void parseParameters(int f_argc, char* f_argv_p[], CIpParameters& 
f_ipParameters) {

    int curDefaultOption=0;

    for (int argNo=1; argNo<f_argc; ++argNo) {

        const std::string curArg=f_argv_p[argNo];

        // help
        if ((curArg=="-help") || (curArg=="--help") || (curArg=="-h") || 
        (curArg=="-?") || (curArg=="--?")) {
            printUsage(f_argv_p[0], false);
            exit(0);
        }

        // threads
        if (curArg=="--threads") {
            if (argNo+1>=f_argc) {
                throw common::CException("Missing argument in "+curArg);
            }

            ++argNo;
            f_ipParameters.m_threads=common::max2<int>(1, common::
            CConverter::strToX<int>(f_argv_p[argNo]));
            continue;            
        }

        // IDW smoothness
        if (curArg=="--idwSmoothness") {
            if (argNo+1>=f_argc) {
                throw common::CException("Missing argument in " + curArg);
            }

            ++argNo;
            f_ipParameters.m_idwSmoothness=common::CConverter::strToX<float>
            (f_argv_p[argNo]);

            continue;
        }

        // bgnd color
        if (curArg=="--bgnd") {
            if (argNo+1>=f_argc) {
                throw common::CException("Missing color channel count in " +
                curArg);
            }

            ++argNo;
            const size_t chCount=abs(common::CConverter::strToX<int32_t>
            (f_argv_p[argNo]));

            // copy bgnd col
            if (argNo+chCount>=f_argc) {
                throw common::CException("Not enough bgnd intensity "
                "values in " +curArg);
            }

            f_ipParameters.m_bgndCol.clear();
            for (size_t c=0; c<chCount; ++c) {
                ++argNo;
                f_ipParameters.m_bgndCol.push_back(common::CConverter::
                strToX<float>(f_argv_p[argNo]));
            }
            
            continue;
        }

        // direction count
        if (curArg=="--dir") {
            if (argNo+1>=f_argc) {
                throw common::CException("Missing argument in " + curArg);
            }

            ++argNo;
            f_ipParameters.m_dirCount=common::CConverter::strToX<data::
            pos_type>(f_argv_p[argNo]);
            continue;
        }

        // disable oversamplig compensation
        if (curArg=="--disableOC") {
            f_ipParameters.m_disableOC=true;
            continue;    
        }

        // angular offset
        if (curArg=="--offset") {
            if (argNo+1>=f_argc) {
                throw common::CException("Missing argument in " + curArg);
            }

            ++argNo;
            f_ipParameters.m_angularOffset=common::CConverter::strToX<float>(
            f_argv_p[argNo]);
            continue;
        }

        // postprocess
        if (curArg=="--postprocess") {
            if (argNo+2>=f_argc) {
                throw common::CException("Missing arguments in " + curArg);
            }

            ++argNo;
            f_ipParameters.m_ppFilterSize.m_i=common::CConverter::strToX<
            data::CPoint2DInt::CoordType>(f_argv_p[argNo]);
            ++argNo;
            f_ipParameters.m_ppFilterSize.m_j=common::CConverter::strToX<
            data::CPoint2DInt::CoordType>(f_argv_p[argNo]);
            continue;
        }

        // default option
        // don't forget to review the sanity check below on changes
        switch (curDefaultOption) {
            case 0:
                // input image
                f_ipParameters.m_inputImageUri=curArg;
            break;
            case 1: 
                // output image
                f_ipParameters.m_outputImageUri=curArg;
            break;
        }

        ++curDefaultOption;
    } // argNo

    //
    // sanity checks
    //

    // announce number of threads to OpenMP, auto-set if negative
    #ifdef _OPENMP
        if (f_ipParameters.m_threads<0) {
            f_ipParameters.m_threads=omp_get_max_threads();
        }

        omp_set_num_threads(f_ipParameters.m_threads); 
    #else
        f_ipParameters.m_threads=1;
    #endif

    // have all mandatory arguments ?
    if (curDefaultOption!=2) {
        throw common::CException("Missing "+common::CConverter::xToStr(
        2-curDefaultOption)+" mandatory argument(s), -h for help");
    }
}


/****************************************************************************/


/**
 * SetupImages loads the input image and sets up the output and weight 
 * images.
 *
 * @param f_ipParameters the interpolation parameter record providing the 
 * image URIs
 * @param f_inputImage the input image
 * @param f_outputImage the output image
 * @param f_weightImage the weight image
 */
bool setupImages(CIpParameters& f_ipParameters, image::CFloatImage&
f_inputImage, image::CFloatImage& f_outImage, image::CFloatImage& 
f_weightImage) {

    // load input image
    try {
        f_inputImage.readTiff(f_ipParameters.m_inputImageUri, false);
    } catch (common::CException& e) {
        std::cout << "Cannot read input image (" << e.what() << ")" <<
        std::endl;
        return false;
    }

    // create output and weight image to have the same properties as
    // the input bitmap
    f_outImage.create(f_inputImage.getWidth(), f_inputImage.getHeight(), 
    f_inputImage.getNrOfChannels());
    f_weightImage.create(f_inputImage.getWidth(), f_inputImage.getHeight(), 
    f_inputImage.getNrOfChannels());

    // check against bgnd color
    if (f_ipParameters.m_bgndCol.size()!=f_inputImage.getNrOfChannels()) {
        if (f_ipParameters.m_bgndCol.size()) {
            throw common::CException("Mismatch between input image channel "
            "count and explicity given background color ("+common::
            CConverter::xToStr(f_inputImage.getNrOfChannels())+" vs. "+
            common::CConverter::xToStr(f_ipParameters.m_bgndCol.
            size())+")");
        } else {
            f_ipParameters.m_bgndCol.resize(f_inputImage.getNrOfChannels(), 
            0);

            std::cout << "Background color will be implicitly zeroed" 
            << std::endl;
        }
    }

    // all OK
    return true;
}


/****************************************************************************/ 


void writeOutputImage(image::CFloatImage& f_outImage, const data::pos_type
f_origSampleType, data::pos_type f_origBpc, const std::string& f_outUri) {    

    std::string outType="float";

    switch (f_origSampleType) {            

        case SAMPLEFORMAT_UINT:
            switch (f_origBpc) {
                case 8: 
                    f_outImage.writeTiff<uint8_t>(f_outUri);
                    outType="uint8";
                break;
                case 16: 
                    f_outImage.writeTiff<uint16_t>(f_outUri);
                    outType="uint16";
                break;
                case 32: 
                    f_outImage.writeTiff<uint32_t>(f_outUri);
                    outType="uint32";
                break;
                case 64: 
                    f_outImage.writeTiff<uint64_t>(f_outUri);
                    outType="uint64";
                break;
                default:
                    // use sample type of output image as of its creation
                    f_outImage.writeTiff(f_outUri);                    
            }
        break;

        case SAMPLEFORMAT_INT:
            switch (f_origBpc) {
                case 8: 
                    f_outImage.writeTiff<int8_t>(f_outUri);
                    outType="int8";
                break;        
                case 16: 
                    f_outImage.writeTiff<int16_t>(f_outUri);
                    outType="int16";
                break;
                case 32: 
                    f_outImage.writeTiff<int32_t>(f_outUri);
                    outType="int32";
                break;
                case 64: 
                    f_outImage.writeTiff<int64_t>(f_outUri);
                    outType="int64";
                break;
                default:
                    // use sample type of output image as of its creation
                    f_outImage.writeTiff(f_outUri);                    
            }
        break;

        case SAMPLEFORMAT_IEEEFP:
            switch (f_origBpc) {
                case 64: 
                    f_outImage.writeTiff<double>(f_outUri);
                    outType="double";
                break;
                default:
                    // use sample type of output image as of its creation
                    f_outImage.writeTiff(f_outUri);                    
            }
        break;

        default:
            // use sample type of output image as of its creation
            f_outImage.writeTiff(f_outUri);
    }

    std::cout << "Output image will be of type " << outType << std::endl;
}


/****************************************************************************/ 


/**
 * Main is the program entry point.
 *
 * @param f_argc the argument count
 * @param f_argv_p the command line arguments
 */
int main(int f_argc, char* f_argv_p[]) {
    
    try {    

        // parse parameters
        printUsage(f_argv_p[0], true);

        CIpParameters ipParameters;
        parseParameters(f_argc, f_argv_p, ipParameters);

        // load, create images
        image::CFloatImage inImage, outImage, weightImage;
        setupImages(ipParameters, inImage, outImage, weightImage);
 
        // print interpolation parameters
        ipParameters.print();        

        //  
        // now run the interpolation
        //        
        common::CStdProgress stdProgress;
        stdProgress.setPercentageSteps(0.01);

        data::CPoint3DInt pointStats;
        std::vector<double> timingStats;
        const data::pos_type totalPixels=inImage.getWidth()*inImage.
        getHeight();

        // OpenMP mode
        image::CInterpolator2DOpenMP ip2d(inImage, outImage, weightImage, 
        ipParameters.m_bgndCol, ipParameters.m_dirCount, ipParameters.
        m_angularOffset, ipParameters.m_idwSmoothness, ipParameters.
        m_disableOC, &stdProgress);

        ip2d.interpolate(ipParameters.m_ppFilterSize, &pointStats,
        &timingStats);

        // point stats
        std::cout << std::endl;        
        std::cout << totalPixels << " input image pixels read ("<< inImage.
        getWidth() << " x " << inImage.getHeight() << ")" << std::endl;
        std::cout << pointStats.m_x+pointStats.m_y+pointStats.m_z << 
        " input image pixels processed" << std::endl;            
        std::cout << pointStats.m_x << " foreground pixels (" << 
        (100.0*pointStats.m_x/totalPixels) << "%)" << std::endl;
        std::cout << pointStats.m_y << " copied background pixels (" << 
        (100.0*pointStats.m_y/totalPixels) << "%)" << std::endl;
        std::cout << pointStats.m_z << " interpolated background pixels ("
        << (100.0*pointStats.m_z/totalPixels) << "%)" << std::endl;

        // timing stats
        std::cout << std::endl;
        if (timingStats.size()>0) {
            std::cout << "weight image initialization took " << 
            timingStats[0] << " ms (" << timingStats[0]/1000 << " s)" << 
            std::endl;
        }

        if (timingStats.size()>0) {
            std::cout << "IDW interpolation took " << timingStats[1] <<
            " ms (" << timingStats[1]/1000 << " s)" << std::endl;
        }

        if (timingStats.size()>2) {
            std::cout << "intensity normalization took " << timingStats[2] <<
            " ms (" << timingStats[2]/1000 << " s)" << std::endl;
        }

        if (timingStats.size()>3) {
            std::cout << "postprocessing took " << timingStats[3] << 
            " ms (" << timingStats[3]/1000 << " s)" << std::endl;
        }

        // save output image, use same type as input
        writeOutputImage(outImage, inImage.getOrigSampleType(), inImage.
        getOrigBpc(), ipParameters.m_outputImageUri);

    } catch (const common::CException& e) {
        std::cout << "Error: " << e.what() << std::endl;
        return -2;
    }   
    
    return 0;
}


/****************************************************************************/
