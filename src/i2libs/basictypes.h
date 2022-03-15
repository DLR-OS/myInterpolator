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

#ifndef DATA_BASIC_TYPES_H
#define DATA_BASIC_TYPES_H

#include <cfloat>


/*
 * Basic data types used for the description of samplesddata (images,
 * audio, video). Build on the top of platform-specific definitions.
 *
 * If you change the typedefs below, make sure that the new types
 * are mutually compatible and keep the original semantics.
 */


namespace data {

//
// Original basic types: 
//
// - allow for precise samples
// - allow huge positions/size arguments, however, the user must take
//   care not to violate the 64/63 bit limit
//

#define COMMON_CPLUSPLUS_VERSION __cplusplus


/**
 * Sample data type, i.e. for pixel components.
 */
typedef double sample_type;
//typedef float sample_type;


/**
 * This one is to be used for integer positions and offsets. We 
 * explicitly allow for negative positions.
 */
typedef int64_t pos_type;


/**
 * Floating-point positions are needed to achive i.e. subpixel precision.
 * Negative values are allowed explicitly.
 */
typedef double fppos_type;


/**
 * Sizes and values that have to be positive.
 */
typedef uint64_t size_type;


} // namespace data

#endif

