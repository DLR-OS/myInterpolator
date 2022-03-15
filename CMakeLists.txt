#
# CMake file for myInterpolator (standalone application)
# (CPU version)
#

cmake_minimum_required(VERSION 3.17)
project(myInterpolator VERSION 1.0)
set(CMAKE_VERBOSE_MAKEFILE ON)
message(STATUS "Static library suffix is " ${CMAKE_STATIC_LIBRARY_SUFFIX})

#
# define directories etc.
#
set(PROJECT_SOURCE_DIR src)
set(SRCPATH ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(BINARY_DEBUG_POSTFIX "d")

#
# OpenMP support (optional, required for multithreading)
#
find_package(OpenMP)
if (OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif (OPENMP_FOUND)

#
# Require dot, treat the other components as optional
#
find_package(Doxygen)

#
# include directories
#

# include custom CMakeModules
set(CMAKE_MODULE_PATH $ENV{CMAKE_MODULE_PATH})

find_package(TIFF REQUIRED)
if (TIFF_FOUND)
  include_directories(${TIFF_INCLUDE_DIR})
endif (TIFF_FOUND)

# src dir and subdirs itself
include_directories(${PROJECT_SOURCE_DIR} ${PROJECT_SOURCE_DIR}/i2libs)
include(CheckIncludeFiles)

#
# categorized source files
#

# i2libs 
set(HEADERS_I2LIBS
    ${SRCPATH}/i2libs/basictypes.h
    ${SRCPATH}/i2libs/Converter.h
    ${SRCPATH}/i2libs/coolmath.h
    ${SRCPATH}/i2libs/Exception.h
    ${SRCPATH}/i2libs/fpcomp_base.h
    ${SRCPATH}/i2libs/fpcomp.h
    ${SRCPATH}/i2libs/Point3D.h
    ${SRCPATH}/i2libs/Progress.h
    ${SRCPATH}/i2libs/StdProgress.h
    ${SRCPATH}/i2libs/StopWatch.h
)

set(SOURCES_I2LIBS
    ${SRCPATH}/i2libs/Converter.cpp
    ${SRCPATH}/i2libs/coolmath.cpp
    ${SRCPATH}/i2libs/StopWatch.cpp
 )

SOURCE_GROUP("Header Files\\i2libs" FILES ${HEADERS_I2LIBS})
SOURCE_GROUP("Source Files\\i2libs" FILES ${SOURCES_I2LIBS})

# myInterpolator
set(HEADERS_MYINTERPOLATOR
    ${SRCPATH}/Image.h
    ${SRCPATH}/Interpolator2DOpenMP.h
    ${SRCPATH}/IpParameters.h
    ${SRCPATH}/Postprocess.h
)

set(SOURCES_MYINTERPOLATOR
    ${SRCPATH}/Interpolator2DOpenMP.cpp
    ${SRCPATH}/Postprocess.cpp
    ${SRCPATH}/myInterpolator.cpp
 )

SOURCE_GROUP("Header Files\\myInterpolator" FILES ${HEADERS_MYINTERPOLATOR})
SOURCE_GROUP("Source Files\\myInterpolator" FILES ${SOURCES_MYINTERPOLATOR})

#
# announce all sources
#

set(HEADERS
    ${HEADERS_I2LIBS}
    ${HEADERS_MYINTERPOLATOR}
)

set(SOURCES    
    ${SOURCES_I2LIBS}
    ${SOURCES_MYINTERPOLATOR}
)

#
# configure binary
#

# debug postfix to distinguish between the builds
set(CMAKE_DEBUG_POSTFIX ${BINARY_DEBUG_POSTFIX})

#
# build applications
#

# myInterpolator
add_executable(myInterpolator ${SOURCES} ${HEADERS})
set_target_properties(myInterpolator PROPERTIES DEBUG_POSTFIX ${BINARY_DEBUG_POSTFIX})
target_link_libraries(
    myInterpolator
    ${TIFF_DEBUG_LIBRARIES} ${JPEG_DEBUG_LIBRARIES} ${ZLIB_DEBUG_LIBRARIES}
    ${TIFF_LIBRARIES} ${JPEG_LIBRARIES} ${ZLIB_LIBRARIES}
)

# build binary distribution
message(STATUS "Install directory will be " ${CMAKE_INSTALL_PREFIX}/bin)
install(TARGETS myInterpolator RUNTIME DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

# documentation
if (DOXYGEN_FOUND)
    doxygen_add_docs(Doxygen ${SRCPATH} ALL COMMENT "Create documentation")
    install(DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/html DESTINATION ${CMAKE_INSTALL_PREFIX}/share)
endif(DOXYGEN_FOUND)

