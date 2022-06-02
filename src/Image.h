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

#ifndef IMAGE_SIMPLEIMAGE_H
#define IMAGE_SIMPLEIMAGE_H

#include <type_traits>
#include "tiffio.h"
#include "Converter.h"
#include "Exception.h"
#include "Point2D.h"
#include "Point3D.h"

// for debugging, enables range checks etc.
#ifdef _DEBUG
    #define _DEBUG_IMAGE_CLASS
#endif
    

namespace image {


    /**
     * CImage stores image bitmaps, and reads or writes them in TIF format.
     * Internally, the bitmaps are stored interleaved with multiple channels
     * allowed. Virtual subimages can be created which share the data with
     * their parent base images.
     *
     * @author Dirk 'jtk' Frommholz
     * @date April 1, 2018
     */
    template<typename _SampleType>
    class CImage {
    public:


        /**
         * The pixel component type.
         */
        typedef _SampleType SampleType;


        /**
         * The iterator type.
         */
        typedef SampleType* Iterator;


    private:


        /**
         * The image width, height and number of channels.
         */
        data::CPoint3DInt m_dim;


        /**
         * The original sample type of the TIFF when it has been read
         * from mass storage (libTIFF constants).
         */
        data::pos_type m_origSampleType;


        /**
         * The original bits per channel of the TIFF when it has been read
         * from mass storage.
         */
        data::pos_type m_origBpc;


        /**
         * The virtual image width, height and number of channels.
         */
        data::CPoint3DInt m_vDim;


        /**
         * The horizontal, vertical and channel offset of the virtual image.
         */
        data::CPoint3DInt m_vOffset;


        /**
         * Flag indicating if an image is a virtual or a base image.
         */
        bool m_isVirtual;
       

        /**
         * The number of samples per row of the non-virtual image.
         */
        data::pos_type m_samplesPerRow;


        /**
         * The image bitmap.
         */
        SampleType* m_bitmap_p;


        /**
         * The IEC Rec. 709 grayscale weights.
         */
        float m_iecWeights_p[3];    


    private:


        /**
         * TranslateToBuffer converts buffered samples from a source 
         * to a target buffer. The source and target data types have to
         * be provided via the template arguments. Samples can
         * be down- and upscaled during conversion. Endianess is not 
         * altered however.
         *
         * @tparam S the sample type of the source buffer
         * @tparam T the sample type of the target buffer
         * @param f_targetBuffer_p the target buffer
         * @param f_origSampleBuffer_p the source buffer as from the image 
         * read
         * @param f_nrOfSamplesToTransfer the number of samples to be 
         * translated
         * @param f_scaleToTarget true if the source samples shall be scaled
         * to fully utilize the target bit depth
         */
        template<typename S, typename T>
        void translateToBuffer(uint8_t* f_targetBuffer_p, 
        uint8_t* f_origSampleBuffer_p, uint64_t f_nrOfSamplesToTransfer,
        bool f_scaleToTarget) {

            double scaleFactor=1.0;

            if (f_scaleToTarget) {
                scaleFactor=pow(2.0, ((long)(sizeof(T))-
                (long)(sizeof(S)))*8);
            }

            // for convenience, convert the strip buffer to the native 
            // sample format
            S* origSampleBuffer_p=reinterpret_cast<S*>(f_origSampleBuffer_p);
            T* targetBuffer_p=reinterpret_cast<T*>(f_targetBuffer_p);

            // assume the read position is aligned on sizeof(T) bytes
            for (uint64_t sample=0; sample<f_nrOfSamplesToTransfer;
            ++sample, ++targetBuffer_p) {
                *targetBuffer_p=static_cast<T>(scaleFactor*origSampleBuffer_p
                [sample]);
            } // sample
        }


        /**
         * TranslateBetweenSamples converts buffered samples from a source 
         * to a target buffer depending on the actual sample types.
         *
         * @param f_targetBuffer_p the target buffer
         * @param f_origSampleBuffer_p the source buffer as from the image 
         * read
         * @param f_nrOfSamplesToTransfer the number of samples to be 
         * translated
         * @param f_srcBpc the source bits per component
         * @param f_sampleFormat the sample format (libTIFF constants)
         * @param f_scaleToTarget true if the source samples shall be scaled
         * to fully utilize the target bit depth
         */
        void translateBetweenSamples(uint8_t* f_targetBuffer_p, 
        uint8_t* f_origSampleBuffer_p, uint64_t f_nrOfSamplesToTransfer,
        uint16_t f_srcBpc, uint16_t f_sampleFmt, bool f_scaleToTarget) {

           // convert buffer to target format        
            switch(f_srcBpc) {                
                case 8:
                    if (f_sampleFmt==SAMPLEFORMAT_UINT) {
                        translateToBuffer<uint8_t, SampleType>(
                        f_targetBuffer_p, f_origSampleBuffer_p, 
                        f_nrOfSamplesToTransfer, f_scaleToTarget);
                    } else {
                       translateToBuffer<int8_t, SampleType>(
                       f_targetBuffer_p, f_origSampleBuffer_p, 
                       f_nrOfSamplesToTransfer, f_scaleToTarget);
                    }
                break;
                case 16:
                    if (f_sampleFmt==SAMPLEFORMAT_UINT) {
                        translateToBuffer<uint16_t, SampleType>(
                        f_targetBuffer_p, f_origSampleBuffer_p, 
                        f_nrOfSamplesToTransfer, f_scaleToTarget);
                    } else {
                        translateToBuffer<int16_t, SampleType>(
                        f_targetBuffer_p, f_origSampleBuffer_p, 
                        f_nrOfSamplesToTransfer, f_scaleToTarget);
                    }
                break;
                case 32:
                    switch (f_sampleFmt) {
                        case SAMPLEFORMAT_UINT: 
                            translateToBuffer<uint32_t, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        case SAMPLEFORMAT_INT:
                            translateToBuffer<int32_t, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        case SAMPLEFORMAT_IEEEFP:
                            translateToBuffer<float, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        default:
                            throw common::CException(
                            "CImage::translateBetweenSamples: Unsupported "
                            "32 bit sample format");
                    }
                break;
                case 64:
                    switch (f_sampleFmt) {
                        case SAMPLEFORMAT_UINT: 
                            translateToBuffer<uint64_t, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        case SAMPLEFORMAT_INT:
                            translateToBuffer<int64_t, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        case SAMPLEFORMAT_IEEEFP:
                            translateToBuffer<double, SampleType>(
                            f_targetBuffer_p, f_origSampleBuffer_p, 
                            f_nrOfSamplesToTransfer, f_scaleToTarget);
                        break;
                        default:
                            throw common::CException(
                            "CImage::translateBetweenSamples: Unsupported "
                            "64 bit sample format");
                    }
                break;
                default:
                    throw common::CException(
                    "CImage::translateBetweenSamples: Unknown pixel type "
                    "conversion");
            }
        }


        /**
         * GrayscaleSummingRow reduces a multispectral row of samples to
         * an intensity by either taking the first channel on 1- or 2-
         * channel images or applying the Rec.709 formula to >=3-ch-images.
         *
         * @param f_src_p the source samples, must be a multiple of
         * f_spp
         * @param f_destIt the iterator to the destination inside the
         * virtual image where the samples ar to be placed at
         * @param f_nrOfSrcSamples the number of samples in f_src_p
         * @param f_spp the samples per pixel in f_src_p
         * @param f_doGrayscale apply grayscale summing either by taking a
         * single channel or by applying the Rec.709 weightings
         */    
        void grayscaleSummingRow(SampleType* f_src_p, Iterator f_destIt, 
        uint32_t f_nrOfSrcSamples, uint16_t f_spp, bool f_doGrayscale) {

            if ((f_src_p==0) || (f_destIt==0) || 
            (f_nrOfSrcSamples==0) || (f_spp==0)) {
                throw common::CException("Cannot process image row with "
                "zero samples per pixel or zero samples or invalid buffers");
            }

            if (f_doGrayscale) {

                // do grayscale summing
                switch (f_spp) {                    
                    case 1:
                    case 2:
                        // <3 channels ? simply take the first one
                        for (uint32_t p=0; p<f_nrOfSrcSamples; p+=f_spp, 
                        incX(f_destIt), f_src_p+=f_spp) {
                            *f_destIt=*f_src_p;
                        }
                    break;
                    default:
                        // put weighted sum of first 3 channels in one 
                        // channel (good for RGBI, RGBA etc.)
                        for (uint32_t p=0; p<f_nrOfSrcSamples; p+=f_spp, 
                        incX(f_destIt), f_src_p+=f_spp) {
                            *f_destIt=static_cast<SampleType>(
                                m_iecWeights_p[0]*(*f_src_p)+
                                m_iecWeights_p[1]*(*(f_src_p+1))+
                                m_iecWeights_p[2]*(*(f_src_p+2))
                            );
                        }
                } // switch

            } else {
                // no grayscale summing, plain copy to dest buffer                   
                for (uint32_t p=0; p<f_nrOfSrcSamples; p+=f_spp, incX(
                f_destIt)) {
                    for (uint16_t s=0; s<f_spp; ++s, ++f_src_p, incCh(
                    f_destIt)) {
                        *f_destIt=*f_src_p;
                    }

                    decCh(f_destIt, f_spp);                        
                }
            }
        }


        /**
         * Swap shallowly swaps the internals between this object and f_rhs.
         *
         * @param f_rhs the right hand side to shallowly swap the 
         * internals with
         */
        void swap(CImage& f_rhs) {
            std::swap(this->m_dim, f_rhs.m_dim);
            std::swap(this->m_vDim, f_rhs.m_vDim);
            std::swap(this->m_vOffset, f_rhs.m_vOffset);
            std::swap(this->m_isVirtual, f_rhs.m_isVirtual);
            std::swap(this->m_samplesPerRow, f_rhs.m_samplesPerRow);
            std::swap(this->m_bitmap_p, f_rhs.m_bitmap_p);
        }


        /**
         * SetDefaults initializes the class member variables. Called by the
         * constructors.
         */
        void setDefaults() {
            m_bitmap_p=0;
            m_isVirtual=false;
            m_samplesPerRow=0;
            m_origSampleType=-1;
            m_origBpc=0;

            // IEC grayscale summing weights
            m_iecWeights_p[0]=0.2126f;
            m_iecWeights_p[1]=0.7152f;
            m_iecWeights_p[2]=0.0722f;
        }


    public:


        /**
         * The constructor creates a new image and fills it with an initial
         * value.
         *
         * @param f_width the non-negative image width in pixels
         * @param f_height the non-negative image height in pixels
         * @param f_nrOfChannels the non-negative number of color channels
         * @param f_initialValue the initial value the color components of
         * the pixels are to be set to
         */
        CImage(data::pos_type f_width=0, data::pos_type f_height=0, data::
        pos_type f_nrOfChannels=0, const SampleType& f_initialValue=0) {
            setDefaults();           
            create(f_width, f_height, f_nrOfChannels, f_initialValue);
        }


        /**
         * The constructor accepting a CPoint3DInt object for the image
         * geometry.
         *
         * @param f_geom the image geometry
         * @param f_initialValue the initial value to fill the pixel color
         * components with
         */
        CImage(const data::CPoint3DInt& f_geom, const SampleType& 
        f_initialValue=0) {
            setDefaults();
            create(f_geom.m_x, f_geom.m_y, f_geom.m_z, f_initialValue);
        }


        /**
         * MakeVirtualImage creates a virtual image from the existing image.
         * The virtual image shares the data of the underlying bitmap and
         * hence will be valid as long as the underlying bitmap exists in
         * its original form.
         *
         * @param f_virtualImage the virtual image, normally be a CImage
         * instance built by the default constructor (if not it will get
         * destroyed first)
         * @param f_width the width of the virtual image
         * @param f_height the height of the virtual image
         * @param f_ch the channel count of the virtual image
         * @param f_offsetX the horizontal offset of the virtual image 
         * inside the underlying image
         * @param f_offsetY the vertical offset of the virtual image 
         * inside the underlying image
         * @param f_offsetCh the color channel offset of the virtual image 
         * inside the underlying image
         */
        void makeVirtualImage(CImage& f_virtualImage, data::pos_type 
        f_width, data::pos_type f_height, data::pos_type f_ch, 
        data::pos_type f_offsetX=0, data::pos_type f_offsetY=0, 
        data::pos_type f_offsetCh=0) {

            // destroy virtual image
            f_virtualImage.~CImage();

            // set dimensions
            f_virtualImage.m_dim=m_dim;
            f_virtualImage.m_vDim.setCoords(f_width, f_height, f_ch);
            f_virtualImage.m_vOffset.setCoords(f_offsetX, f_offsetY, 
            f_offsetCh);
                        
            // samples per row, storage is from the underlying image
            // slicing is done inside the iterator functions
            f_virtualImage.m_samplesPerRow=m_samplesPerRow;       
            f_virtualImage.m_bitmap_p=m_bitmap_p;

            // mark as virtual and return the shallow copy
            f_virtualImage.m_isVirtual=true;
        }


        /**
         * The assignment operator assigns a new image to this one.
         * Implemented using the copy constructor and swap idiom.
         *
         * @param f_other the new image
         */
        CImage& operator=(CImage f_other) {
            f_other.swap(*this); 
            return *this;
        }


        /**
         * The copy constructor creates a new non-virtual image from 
         * another virtual or non-virtual image.
         *
         * @param f_other the image to construct this one from
         */
        CImage(const CImage& f_other) {

            setDefaults();
            create(f_other.getWidth(), f_other.getHeight(), 
            f_other.getNrOfChannels());

            // copy content to this image
            f_other.copy(*this);
        }


        /**
         * Copy copies the data from one image to another one with
         * the same geometry and channel count. If the images have
         * different properties, the least common denominator will
         * be copied.
         *
         * @param f_dest the target image
         */
        void copy(CImage& f_dest) const {

            // copy this image to dest
            Iterator srcIt=getIterator();
            Iterator destIt=f_dest.getIterator();

            // least common denominator properties
            const data::CPoint3DInt lcd(
                common::min2(getWidth(), f_dest.getWidth()),
                common::min2(getHeight(), f_dest.getHeight()),
                common::min2(getNrOfChannels(), f_dest.getNrOfChannels())
            );

            for (data::pos_type y=0; y<lcd.m_y; ++y, incY(srcIt), f_dest.
            incY(destIt)) {
                for (data::pos_type x=0; x<lcd.m_x; ++x, incX(srcIt), 
                f_dest.incX(destIt)) {
                    for (data::pos_type ch=0; ch<lcd.m_z; ++ch, incCh(srcIt),
                    f_dest.incCh(destIt)) {
                        *destIt=*srcIt;
                    }

                    f_dest.decCh(destIt, lcd.m_z); 
                    decCh(srcIt, lcd.m_z); 
                }

                f_dest.decX(destIt, lcd.m_x); 
                decX(srcIt, lcd.m_x); 
            }
        }


        /**
         * Clear sets all components of all pixels to the same value.
         *
         * @param f_value the value to be set for all color components of
         * all image pixels
         */
        void clear(const SampleType& f_value=0) {

            Iterator it=getIterator();

            for (data::pos_type y=0; y<getHeight(); ++y, incY(it)) {
                for (data::pos_type x=0; x<getWidth(); ++x, incX(it)) {
                    for (data::pos_type ch=0; ch<getNrOfChannels(); ++ch, 
                    incCh(it)) {
                        *it=f_value;    
                    }

                    decCh(it, getNrOfChannels()); 
                }

                decX(it, getWidth()); 
            }
        }   


        /**
         * AdjustBrightnessContrast changes the brightness and contrast of
         * the image with saturation.
         *
         * @param f_scale the value each pixel component is to be multiplied
         * with, adjusts the contrast (executed first)
         * @param f_offset the value to be added to each pixel component, 
         * adjusts the brightness
         * @param f_min the minimum value to truncate underflows to
         * @param f_max the maximum value to truncate overflows to
         */
        void adjustBrightnessContrast(const double& f_scale, const
        double& f_offset, const SampleType& f_min, const SampleType& f_max) {

            Iterator it=getIterator();

            for (data::pos_type y=0; y<getHeight(); ++y, incY(it)) {
                for (data::pos_type x=0; x<getWidth(); ++x, incX(it)) {
                    for (data::pos_type ch=0; ch<getNrOfChannels(); ++ch, 
                    incCh(it)) {

                        const double val=static_cast<double>(*it+*it*f_scale+
                        f_offset);    

                        // saturate
                        *it=(val<f_min) ? f_min:static_cast<SampleType>(val);
                        *it=(val>f_max) ? f_max:static_cast<SampleType>(val);
                    }

                    decCh(it, getNrOfChannels()); 
                }

                decX(it, getWidth()); 
            }
        }   


        /**
         * Create creates a new image replacing the old bitmap. 
         *
         * @param f_width the new image width
         * @param f_height the new image height
         * @param f_nrOfChannels the number of channels of the new image
         * @param f_initialValue the initial value of the color components
         * of the new image
         *
         */
        void create(data::pos_type f_width, data::pos_type f_height, 
        data::pos_type f_nrOfChannels, const SampleType& f_initialValue=0) {

            // virtual images must be create with another function
            if (isVirtual()) {
                throw common::CException("Cannot reallocate virtual image");
            }

            // no virtual image, (re-)create            
            if (m_bitmap_p) {
                delete[] m_bitmap_p;
                m_bitmap_p=0;
            }

            const size_t imageBytes=static_cast<size_t>(
            f_width*f_height*f_nrOfChannels);

            if (imageBytes>0) {

                m_bitmap_p=new SampleType[imageBytes];

                // update metadata
                m_dim.setCoords(f_width, f_height, f_nrOfChannels);
                m_vDim=m_dim;
                m_samplesPerRow=m_dim.m_x*m_dim.m_z;

                // init
                clear(f_initialValue);
            }
        }


        /**
         * IsVirtual returns true if an image is a virtual image, that is,
         * a subset of another image bitmap.
         *
         * @return true if the image is virtual, false if not
         */
        bool isVirtual() const {
            return m_isVirtual;
        }


        /**
         * GetPixel returns the intensity of a particular color component of
         * a pixel.
         *
         * @param f_x the horizontal position in pixels inside the image
         * @param f_y the vertical position in pixels inside the image
         * @param f_ch the color component index
         * @return the intensity of the color component of the pixel
         */
        SampleType getPixel(data::pos_type f_x, data::pos_type f_y, 
        data::pos_type f_ch) const {
            return *getIterator(f_x, f_y, f_ch);
        }


        /**
         * GetPixel returns the intensity of a particular color component of
         * a pixel.
         *
         * @param f_pos the horizontal and vertical position in pixels 
         * inside the image
         * @param f_ch the color component index
         * @return the intensity of the color component of the pixel
         */
        SampleType getPixel(const data::CPoint2DInt& f_pos, data::pos_type 
        f_ch=0) const {
            return getPixel(f_pos.m_i, f_pos.m_j, f_ch);
        }


        /**
         * GetPixel returns the intensity of a particular color component of
         * a pixel.
         *
         * @param f_it the position and color component index
         * @return the intensity of the color component of the pixel
         */
        SampleType getPixel(const data::CPoint3DInt& f_it) 
        const {
            return *getIterator(f_it);
        }


        /**
         * SetPixel sets the intensity of a color component of a pixel
         * inside the image.
         *
         * @param f_it the position and color component index
         * @param f_val the intensity to be set
         */
        void setPixel(const data::CPoint3DInt& f_it, const SampleType& 
        f_val) {
            *getIterator(f_it)=f_val;
        }


        /**
         * SetPixel sets the intensity of a color component of a pixel
         * inside the image.
         *
         * @param f_x the horizontal position in pixels inside the image
         * @param f_y the vertical position in pixels inside the image
         * @param f_ch the color component index
         * @param f_val the intensity to be set
         */
        void setPixel(data::pos_type f_x, data::pos_type f_y, data::pos_type
        f_ch, SampleType f_val) {
            *getIterator(f_x, f_y, f_ch)=f_val;
        }


        /**
         * GetNrOfChannels returns the number of channels of the 
         * virtual image.
         *
         * @return the number of channels of the image
         */        
        data::pos_type getNrOfChannels() const {
            return m_vDim.m_z;
        }


        /**
         * GetWidth returns the width of the virtual image.
         *
         * @return the width of the virtual image
         */
        data::pos_type getWidth() const {
            return m_vDim.m_x;
        }


        /**
         * GetHeight returns the height of the virtual image.
         *
         * @return the height of the virtual image
         */
        data::pos_type getHeight() const {
            return m_vDim.m_y;
        }


        /**
         * GetDim returns the dimensions of the virtual image, i.e.
         * its width, height and color channel count.
         *
         * @return the image dimensions of the virtual image
         */
        const data::CPoint3DInt& getDim() const {
            return m_vDim;
        }        


        /**
         * GetDimXY returns the width and height of the image as
         * a 2D point.
         *
         * @return the 2D image dimension of the virtual image
         */
        const data::CPoint2DInt getDimXY() const {
            return data::CPoint2DInt(m_vDim.m_x, m_vDim.m_y);
        }       


        /**
         * GetOrigSampleType returns the sample type of an image 
         * loaded from mass storage as a libTIFF sample type.
         * This value will be negative if the image has not been
         * successfully loaded before.
         *
         * @return the original sample type as of libTIFF, or
         * a negative value if the image has not been read 
         * successfully
         */
        const data::pos_type getOrigSampleType() const {
            return m_origSampleType;
        }


        /**
         * GetOrigBpc returns the bits per color channel of an image 
         * loaded from mass storage. This value will be zero when
         * the image has not been successfully loaded before.
         *
         * @return the original bpc, or zero if the image has not been read 
         * successfully
         */
        const data::pos_type getOrigBpc() const {
            return m_origBpc;
        }


        /**
         * IsInImage returns true if the coordinates and color channel index
         * are within the range of the virtual image.
         *
         * @param f_x the horizontal position
         * @param f_y the vertical position
         * @param f_ch the color channel
         * @return true if the position and color channel are within the
         * limits of the image
         */
        bool isInImage(data::pos_type f_x, data::pos_type f_y, data::pos_type
        f_ch) const {
            return ((f_x>=0) && (f_x<getWidth()) && (f_y>=0) && (f_y<
            getHeight()) && (f_ch>=0) && (f_ch<getNrOfChannels()));
        }

    
        /**
         * IsInImage returns true if the coordinates passed are within 
         * the range of the virtual image.
         *
         * @param f_x the horizontal position
         * @param f_y the vertical position
         * @return true if the position is within the limits of the image
         */
        bool isInImage(int64_t f_x, int64_t f_y) const {
            return ((f_x>=0) && (f_x<getWidth()) && (f_y>=0) && (f_y<
            getHeight()));
        }


        /**
         * IsInBaseImage returns true if the coordinates and color channel 
         * passed are within the range of the underlying base image of the 
         * virtual image. Since the virtual image can be a smaller subset of
         * the base image this function may return true even though the
         * addressed intensity value is outside the virtual image.
         *
         * @param f_x the horizontal position
         * @param f_y the vertical position
         * @param f_ch the color channel
         * @return true if the position is within the limits of the image
         */
        bool isInBaseImage(data::pos_type f_x, data::pos_type f_y, 
        data::pos_type f_ch) const {
            return ((f_x>=0) && (f_x<m_dim.m_x) && (f_y>=0) && (f_y<
            m_dim.m_y) && (f_ch>=0) && (f_ch<m_dim.m_z));
        }


        /**
         * ProbeTiff will read the TIFF header and return basic image
         * parameters. All output parameters are optional hence just 
         * probing for the format is possible.
         * 
         * @param f_fileName the file name of the image to be probed
         * @param f_dim_p the image dimension found
         * @param f_bpc_p the bits per component found
         * @param f_sampleFormat_p the sample format (libTIFF constants)
         * @return true if the image has been probed successfully, false
         * on error (e.g. wrong file type, read error etc.)
         */
        static bool probeTiff(const std::string& f_fileName, data::
        CPoint3DInt* f_dim_p=0, data::pos_type* f_bpc_p=0, data::pos_type*
        f_sampleFormat_p=0) {

            TIFF* tif_p=TIFFOpen(f_fileName.c_str(), "r");
            if (tif_p) {
                
	            uint32_t imageHeight, imageWidth;
                uint16_t config, spp, bpc;
                uint16_t samplefmt=SAMPLEFORMAT_UINT;

                TIFFGetField(tif_p, TIFFTAG_IMAGEWIDTH, &imageWidth);
	            TIFFGetField(tif_p, TIFFTAG_IMAGELENGTH, &imageHeight);
	            TIFFGetField(tif_p, TIFFTAG_PLANARCONFIG, &config);
                TIFFGetField(tif_p, TIFFTAG_SAMPLESPERPIXEL, &spp);
                TIFFGetField(tif_p, TIFFTAG_BITSPERSAMPLE, &bpc);
                TIFFGetField(tif_p, TIFFTAG_SAMPLEFORMAT, &samplefmt);

                // fill output variables
                if (f_dim_p) {
                    f_dim_p->setCoords(imageWidth, imageHeight, spp);
                }

                if (f_bpc_p) {
                    *f_bpc_p=bpc;
                }

                if (f_sampleFormat_p) {
                    *f_sampleFormat_p=samplefmt;
                }

   	            TIFFClose(tif_p);
                return true;

            } 

            // probe error
            return false;
        }


        /**
         * ReadTiff reads a TIFF image from the given URI and sets up
         * an image for it. The TIFF samples can be scaled to the
         * full width of the target data type, however, no signed/unsigned
         * conversion will be done so expect errors if signed data goes
         * to unsigned pixel types.
         *
         * @param f_fileName the TIFF URI
         * @param f_scaleToTarget true if the samples are to be scaled to the
         * target sample type of the image, false if plain reading shall
         * be done
         * @param f_grayscale true if grayscale summing is to be applied to
         * the input image producing a one-channel intensity image of the
         * same size
         * @param f_create create a new image axactly as the one to be load
         * in case of a physical image; if false load into the existing image
         * (careful!)
         */
        void readTiff(const std::string& f_fileName, bool 
        f_scaleToTarget=false, bool f_grayscale=false, bool f_create=true) {

            // enable strip chopping
            TIFF* tif_p=TIFFOpen(f_fileName.c_str(), "rC");
            if (tif_p) {
                
	            uint32_t imageHeight, imageWidth;
                uint16_t config, spp, bpc;
                uint16_t samplefmt=SAMPLEFORMAT_UINT;
	            uint8_t* buf_p, *gsbuf_p;

                TIFFGetField(tif_p, TIFFTAG_IMAGEWIDTH, &imageWidth);
	            TIFFGetField(tif_p, TIFFTAG_IMAGELENGTH, &imageHeight);
	            TIFFGetField(tif_p, TIFFTAG_PLANARCONFIG, &config);
                TIFFGetField(tif_p, TIFFTAG_SAMPLESPERPIXEL, &spp);
                TIFFGetField(tif_p, TIFFTAG_BITSPERSAMPLE, &bpc);
                TIFFGetField(tif_p, TIFFTAG_SAMPLEFORMAT, &samplefmt);
	            buf_p=(uint8_t*)_TIFFmalloc(TIFFScanlineSize(tif_p));
                // buffer to hold pixels before grayscaling
                gsbuf_p=(uint8_t*)_TIFFmalloc(imageWidth*spp*
                sizeof(SampleType));

                // create internal image if not virtual and requested
                if (m_isVirtual || !f_create) {
                    if ((imageWidth>getDim().m_x) || (imageHeight>getDim().
                    m_y) || ((f_grayscale ? 1: spp)>getDim().m_z)) {
                        throw common::CException("Cannot load "+f_fileName+
                        " into virtual image, source image exceeds virtual "
                        "image dimensions");
                    }
                } else {
                    create(imageWidth, imageHeight, f_grayscale ? 1: spp);
                }

	            if (config==PLANARCONFIG_CONTIG) {
	                for (uint32_t row=0; row<imageHeight; ++row) {
		                TIFFReadScanline(tif_p, buf_p, row);

                        // scale/transform to SampleType
                        translateBetweenSamples((uint8_t*)gsbuf_p, 
                        buf_p, spp*imageWidth, bpc, samplefmt, 
                        f_scaleToTarget);

                        // grayscale summing or plain copying to target
                        // bitmap
                        grayscaleSummingRow((SampleType*)gsbuf_p, 
                        getIterator(0, row, 0), spp*imageWidth, spp, 
                        f_grayscale);
                    }
	            } 

                if (config==PLANARCONFIG_SEPARATE) {

                    // grayscale summing with planar images will consider
                    // the first three channels only, or the first one if 
                    // less
                    // weighting happens on the fly later
                    const uint16_t planesToRead=f_grayscale ?
                    (spp<3 ? 1 : 3) : spp;

                    // now interleave the planar color information
                    // according to libtiff
	                for (uint16_t s=0; s<planesToRead; ++s) {
		                for (uint32_t row=0; row<imageHeight; ++row) {
		                    TIFFReadScanline(tif_p, buf_p, row, s);

                            // slow but must be done this way in case of
                            // separate planes, see libTiff hint
                            translateBetweenSamples((uint8_t*)gsbuf_p,
                            buf_p, m_samplesPerRow, bpc, samplefmt,
                            f_scaleToTarget);

                            // put the currently read channel to the
                            // right image position - plain copy
                            Iterator it=getIterator(0, row, f_grayscale ?
                            0 : s);

                            for (uint32_t x=0; x<imageWidth; ++x) {

                                if (f_grayscale && (planesToRead==3)) {
                                    // weighting on-the-fly: init on 1st
                                    // component, weight-add on 2nd/3rd
                                    //
                                    // we may loose precision here
                                    *it=(s==0) ? 
                                    static_cast<SampleType>(*(gsbuf_p+x)*
                                    m_iecWeights_p[s]) :
                                    static_cast<SampleType>(*it+*(gsbuf_p+
                                    x)*m_iecWeights_p[s]);
                                } else {
                                    // plain copy    
                                    *it=*(gsbuf_p+x);
                                }

                                incX(it);
                            }                                
	                    } // for row
                    } // for s
                } // contiguousSamples

                // alter original sample type etc. on success only
                m_origBpc=bpc;
                m_origSampleType=samplefmt;

                _TIFFfree(gsbuf_p);
	            _TIFFfree(buf_p);
	            TIFFClose(tif_p);
            } else {
                throw common::CException("Cannot open image "+f_fileName+
                " for reading");
            } // tif_p
        }


        /**
         * WriteTiff writes the image as a TIFF bitmap.
         *
         * @tparam OutType the sample data type to be used
         * @param f_fileName the filename of the destination bitmap
         *
         */
        template<typename OutType=SampleType>
        void writeTiff(std::string f_fileName) {   

            // do BigTIFF ?
            const std::string writeFormatStr=getImageBytes()>4294967295ULL ?
            "w8" : "w4";

            TIFF* tif_p=TIFFOpen(f_fileName.c_str(), writeFormatStr.c_str());

            if (!tif_p) {
                throw common::CException("Cannot open image "+f_fileName+
                " for writing");
            }

            TIFFSetField(tif_p, TIFFTAG_IMAGEWIDTH, getWidth());
            TIFFSetField(tif_p, TIFFTAG_IMAGELENGTH, getHeight());
            TIFFSetField(tif_p, TIFFTAG_SAMPLESPERPIXEL, getNrOfChannels());
            TIFFSetField(tif_p, TIFFTAG_BITSPERSAMPLE, sizeof(OutType)*8);
            TIFFSetField(tif_p, TIFFTAG_ORIENTATION, ORIENTATION_TOPLEFT);
            TIFFSetField(tif_p, TIFFTAG_PLANARCONFIG, PLANARCONFIG_CONTIG);
            TIFFSetField(tif_p, TIFFTAG_PHOTOMETRIC, getNrOfChannels()!=1 ? 
            PHOTOMETRIC_RGB : PHOTOMETRIC_MINISBLACK);

            // sample format
            uint64_t tiffSampleFormat=SAMPLEFORMAT_UINT;
    
            if ( (typeid(OutType)==typeid(int8_t)) ||
            (typeid(OutType)==typeid(int16_t)) ||
            (typeid(OutType)==typeid(int32_t)) ||
            (typeid(OutType)==typeid(int64_t)) ) {
                tiffSampleFormat=SAMPLEFORMAT_INT;
            }

            if ( (typeid(OutType)==typeid(float)) ||
            (typeid(OutType)==typeid(double)) ) {
                tiffSampleFormat=SAMPLEFORMAT_IEEEFP;
            }

            TIFFSetField(tif_p, TIFFTAG_COMPRESSION, COMPRESSION_NONE);
            TIFFSetField(tif_p, TIFFTAG_SAMPLEFORMAT, tiffSampleFormat);    

            // add some merits
            const std::string softwareString=
            (std::string("Simple CImage class (LibTiff ")+common::
            CConverter::xToStr<uint32_t>(TIFFLIB_VERSION)+")").c_str();
            TIFFSetField(tif_p, TIFFTAG_SOFTWARE, softwareString.c_str());

            const std::string artistString("jTk!");
            TIFFSetField(tif_p, TIFFTAG_ARTIST, artistString.c_str());

            // write strips
            tsize_t lineBytes=getNrOfChannels()*getWidth()*sizeof(
            OutType);
            TIFFSetField(tif_p, TIFFTAG_ROWSPERSTRIP, 1);

            // now writing image to the file one strip at a time
            //
            // can we use plain writing ?
            // this requires the virtual image channel count to be equal to
            // its base image
            const bool havePlainWriting=(typeid(OutType)==typeid(SampleType))
            && (getNrOfChannels()==m_dim.m_z); 
            OutType* scanline_p=0;

            for (int32_t y=0; y<getHeight(); ++y) {

                if (havePlainWriting) {          
                    // getIterator() will take into account any offsets of
                    // the virtual image, so this is simple          
                    // must use the cast to make the compiler not complain
                    // about mismatching template args; in this case this
                    // statement never will get called
                    scanline_p=reinterpret_cast<OutType*>(getIterator(0, y, 0));
                } else {
                    // compose scanline
                    if (!scanline_p) {
                        scanline_p=new OutType[getWidth()*getNrOfChannels()];
                    }

                    // copy offset color components to local scanline (slow!)
                    for (data::pos_type x=0; x<getWidth(); ++x) {
                        for (data::pos_type ch=0; ch<getNrOfChannels(); 
                        ++ch) {
                            SampleType sample=getPixel(x, y, ch);

                            // saturate to output type
                            if ((typeid(OutType)==typeid(float)) || (typeid(
                            OutType)==typeid(double))) {
                                if (sample<-std::numeric_limits<OutType>::max()) {
                                    sample=-std::numeric_limits<OutType>::max();
                                }

                                if (sample>std::numeric_limits<OutType>::max()) {
                                    sample=std::numeric_limits<OutType>::max();
                                }
                            } else {
                                if (sample<std::numeric_limits<OutType>::min()) {
                                    sample=std::numeric_limits<OutType>::min();
                                }

                                if (sample>std::numeric_limits<OutType>::max()) {
                                    sample=std::numeric_limits<OutType>::max();
                                }
                            }

                            scanline_p[x*getNrOfChannels()+ch]=static_cast<
                            OutType>(sample);
                        }
                    }
                }

                if (TIFFWriteScanline(tif_p, scanline_p, y, 0)<0) {
                    throw common::CException("Failed to write virtual "
                    "TIF row #"+common::CConverter::xToStr(y)+
                    " mapping to physical row #"+common::CConverter::
                    xToStr(y+m_vOffset.m_y)+" to file "+f_fileName);
                }
            }

            // clean up if we had to compose the scanline content locally
            if ((!havePlainWriting) && (scanline_p)) {
                delete[] scanline_p;
            }

            TIFFClose(tif_p);
        }


        /**
         * GetImageBytes returns the number of bytes consumed by the 
         * virtual or underlying image.
         *
         * @param f_computeForPhysImage true if the memory consumption is
         * to be computed for the underlying physical image not the virtual
         * image parameters
         * @return the memory consumption
         */
        const uint64_t getImageBytes(const bool f_computeForPhysImage=
        false) const {
            if (f_computeForPhysImage) {
                return m_dim.m_x*m_dim.m_y*m_dim.m_z*sizeof(SampleType);
            } else {
                return getWidth()*getHeight()*getNrOfChannels()*sizeof(
                SampleType);
            }
        }

        //
        // resampling
        //


        /**
         * Resample resamples the image to fit the destination image which
         * must have been created already. The destination  image may 
         * have a different pixel data type, however, plain copying with 
         * rounding is applied. This function uses a bilinear resampling
         * filter.
         *
         * @param f_destImage the destination image, must have been created
         * already
         */
        template<typename DestSampleType>
        void resample(CImage<DestSampleType>& f_destImage) {            

            for (data::pos_type y=0; y<f_destImage.getHeight(); ++y) {

                // compute the two vertical positions inside the source 
                // image where to read from and interpolate
                const data::pos_type y0=(y*getHeight())/f_destImage.
                getHeight();
                const data::pos_type y1=(y0+1<getHeight()) ? y0+1 : y0;

                // the vertical color weights
                const double y1w=((static_cast<double>(y)*getHeight())/
                f_destImage.getHeight())-y0;
                const double y0w=1.0-y1w;

                for (data::pos_type x=0; x<f_destImage.getWidth(); ++x) {

                    // compute the two horizontal positions inside the source
                    // image where to read from and interpolate
                    const data::pos_type x0=(x*getWidth())/f_destImage.
                    getWidth();
                    const data::pos_type x1=(x0+1<getWidth()) ? x0+1 : x0;

                    // the horizontal color weights
                    const double x1w=((static_cast<double>(x)*getWidth())/
                    f_destImage.getWidth())-x0;
                    const double x0w=1.0-x1w;

                    for (data::pos_type ch=0; ch<f_destImage.
                    getNrOfChannels(); ++ch) {

                        // get the source color
                        const DestSampleType color=
                        static_cast<DestSampleType>(
                            x0w*y0w*getPixel(x0, y0, ch)+
                            x0w*y1w*getPixel(x0, y1, ch)+
                            x1w*y0w*getPixel(x1, y0, ch)+
                            x1w*y1w*getPixel(x1, y1, ch)
                        );

                        f_destImage.setPixel(x, y, ch, color);
                    } // ch
                } // x
            } // y
        }


        //
        // iterator interface - no range checking
        // operate on the non-virtual base image
        //


        /**
         * GetIterator returns an iterator to a pixel.
         *
         * @param f_x the initial x coordinate where the iterator is to
         * be placed on
         * @param f_y the initial y coordinate where the iterator is to
         * be placed on
         * @param f_ch the initial channel where the iterator is to
         * be placed on
         */
        Iterator getIterator(const data::pos_type f_x=0, const data::pos_type 
        f_y=0, const data::pos_type f_ch=0)  {

            #ifdef _DEBUG_IMAGE_CLASS
                if (!isInBaseImage(f_x, f_y, f_ch)) {
                    throw common::CException(
                    "CImage::getIterator (const): Accessing pixel x;y;ch="+
                    common::CConverter::xToStr(f_x)+";"+common::CConverter::
                    xToStr(f_y)+";"+common::CConverter::xToStr(f_ch)+
                    " is outside the physical image dimension of dimX;dimY;"
                    "dimZ="+common::CConverter::xToStr(m_dim.m_x)+";"+
                    common::CConverter::xToStr(m_dim.m_y)+";"+common::
                    CConverter::xToStr(m_dim.m_z));
                }
            #endif

            return m_bitmap_p+((m_dim.m_x*(f_y+m_vOffset.m_y)+f_x+
            m_vOffset.m_x)*m_dim.m_z+f_ch+m_vOffset.m_z);
        }


        /**
         * GetIterator returns a read-only iterator to the pixel of the 
         * virtual image given by its coordinates.
         *
         * @param f_x the initial horizontal position where the iterator is 
         * to be placed at
         * @param f_y the initial vertical position where the iterator is 
         * to be placed at
         * @param f_ch the color channel to place the iterator at
         */
        const Iterator getIterator(const data::pos_type f_x=0, const data::pos_type 
        f_y=0, const data::pos_type f_ch=0) const {

            #ifdef _DEBUG_IMAGE_CLASS
                if (!isInBaseImage(f_x, f_y, f_ch)) {
                    throw common::CException(
                    "CImage::getIterator (const): Accessing pixel x;y;ch="+
                    common::CConverter::xToStr(f_x)+";"+common::CConverter::
                    xToStr(f_y)+";"+common::CConverter::xToStr(f_ch)+
                    " is outside the physical image dimension of dimX;dimY;"
                    "dimZ="+common::CConverter::xToStr(m_dim.m_x)+";"+
                    common::CConverter::xToStr(m_dim.m_y)+";"+common::
                    CConverter::xToStr(m_dim.m_z));
                }
            #endif

            return m_bitmap_p+((m_dim.m_x*(f_y+m_vOffset.m_y)+f_x+
            m_vOffset.m_x)*m_dim.m_z+f_ch+m_vOffset.m_z);
        }


        /**
         * GetIterator returns a read-only iterator to the aggregate
         * pixel of the virtual image.
         *
         * @param f_p the initial position where the iterator is to be 
         * placed on
         */
        const Iterator getIterator(const data::CPoint3DInt& f_p) const {
            return getIterator(f_p.m_x, f_p.m_y, f_p.m_z);
        }


        /**
         * GetIterator returns an iterator to the aggregate pixel of the 
         * virtual image.
         *
         * @param f_p the initial position where the iterator is to be 
         * placed on
         */
        Iterator getIterator(const data::CPoint3DInt& f_p) {
            return getIterator(f_p.m_x, f_p.m_y, f_p.m_z);
        }


        /**
         * GetIterator returns a read-only iterator to the given pixel and 
         * channel of the virtual image.
         *
         * @param f_p the initial position where the iterator is to be 
         * placed on
         * @param f_ch the initial color channel
         */
        const Iterator getIterator(const data::CPoint2DInt& f_p, const data::
        pos_type f_ch=0) const {
            return getIterator(f_p.m_i, f_p.m_j, f_ch);
        }


        /**
         * GetIterator returns a writable iterator to the given pixel and 
         * channel of the virtual image.
         *
         * @param f_p the initial position where the iterator is to be 
         * placed on
         * @param f_ch the initial color channel
         */
        Iterator getIterator(const data::CPoint2DInt& f_p, const data::
        pos_type f_ch=0) {
            return getIterator(f_p.m_i, f_p.m_j, f_ch);
        }


        /**
         * IncX increments the horizontal position of the iterator
         * inside the virtual image
         *
         * @param f_it the iterator whose horizontal position is to 
         * be incremented
         */
        void incX(Iterator& f_it) const {
            f_it+=m_dim.m_z;
        }


        /**
         * IncX increments the horizontal position of the iterator by
         * the given number of pixels.
         *
         * @param f_it the iterator whose horizontal position is to be 
         * incremented
         * @param f_add the number of pixels to go back horizontally
         */
        void incX(Iterator& f_it, const data::pos_type f_add) const {
            f_it+=(m_dim.m_z*f_add);
        }


        /**
         * DecX decrements the horizontal position of the iterator
         * inside the virtual image
         *
         * @param f_it the iterator whose horizontal position is to be 
         * decremented
         */
        void decX(Iterator& f_it) const {
            f_it-=m_dim.m_z;
        }


        /**
         * DecX decrements the horizontal position of the iterator by
         * the given number of pixels.
         *
         * @param f_it the iterator whose horizontal position is to be 
         * decremented
         * @param f_sub the number of pixels to go back horizontally
         */
        void decX(Iterator& f_it, const data::pos_type f_sub) const {
            f_it-=(m_dim.m_z*f_sub);
        }


        /**
         * IncY increments the vertical position of the iterator.
         *
         * @param f_it the iterator whose vertical position is to
         * be incremented
         */
        void incY(Iterator& f_it) const {
            f_it+=m_samplesPerRow;
        }


        /**
         * IncY increments the vertical position of the iterator by the
         * given number of pixels.
         *
         * @param f_it the iterator whose vertical position is to
         * be incremented
         * @param f_add the number of rows to add
         */
        void incY(Iterator& f_it, const data::pos_type f_add) const {
            f_it+=(m_samplesPerRow*f_add);
        }


        /**
         * DecY decrements the vertical position of the iterator.
         *
         * @param f_it the iterator whose vertical position is to be 
         * decremented
         */
        void decY(Iterator& f_it) const {
            f_it-=m_samplesPerRow;
        }


        /**
         * DecY decrements the vertical position of the iterator by the
         * given number of pixels.
         *
         * @param f_it the iterator whose vertical position is to
         * be decremented
         * @param f_sub the number of rows to subtract
         */
        void decY(Iterator& f_it, const data::pos_type f_sub) const {
            f_it-=(m_samplesPerRow*f_sub);
        }


        /**
         * IncCh increments the channel of the iterator.
         *
         * @param f_it the iterator the channel is to be incremented for
         */
        void incCh(Iterator& f_it) const {
            ++f_it;
        }


        /**
         * DecCh decrements the channel of the iterator.
         *
         * @param f_it the iterator the channel is to be decremented for
         */
        void decCh(Iterator& f_it) const {
            --f_it;
        }


        /**
         * IncCh increments the channel of the iterator by the given number.
         *
         * @param f_it the iterator the channel is to be incremented for
         * @param f_add the increment value
         */
        void incCh(Iterator& f_it, const data::pos_type f_add) const {
            f_it+=f_add;
        }


        /**
         * DecCh decrements the channel of the iterator by the given number.
         *
         * @param f_it the iterator the channel is to be decremented for
         * @param f_sub the decrement value
         */
        void decCh(Iterator& f_it, const data::pos_type f_sub) const {
            f_it-=f_sub;
        }


        /**
         * The destructor
         */
        ~CImage() {

            // do not delete bitmap if virtual, this bitmap will
            // get destroyed when the underlying base image gets
            // destroyed
            if (m_bitmap_p && (!isVirtual())) {
                delete[] m_bitmap_p;
                m_bitmap_p=0;
            }
        }
    };


/**
 * The data type for images with 8 bits per pixel, unsigned.
 */
typedef CImage<uint8_t> CU8Image;


/**
 * The data type for images with 16 bits per pixel, unsigned.
 */
typedef CImage<uint16_t> CU16Image;


/**
 * The data type for single-precision floating point images.
 */
typedef CImage<float> CFloatImage;

} // namespace image

#endif // IMAGE_SIMPLEIMAGE_H
