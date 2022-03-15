#!/bin/bash

#
# build and install myInterpolator into a directory passed by the caller (Linux)
#

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo "usage: $0 <buildDir> <installDir>"
    exit -1
fi

BATDIR="$(dirname "`readlink -f "$0"`")"
BUILDDIR="$1"
INSTALLDIR="$2"
CMAKEMODULESDIR="${BATDIR}/CMakeModules"

# test if install directory exists
if [ ! -d "$INSTALLDIR" ]; then
    echo "Installation direction $INSTALLDIR must exist already"
    exit -2
fi

# temporarily enhance CMAKE_MODULE_PATH
OLD_CMAKE_MODULE_PATH="${CMAKE_MODULE_PATH}"
if [ "${CMAKE_MODULE_PATH}" == "" ]; then
    export CMAKE_MODULE_PATH="$CMAKEMODULESDIR"
else
    export CMAKE_MODULE_PATH="${CMAKE_MODULE_PATH};${CMAKEMODULESDIR}"
fi

# test for build system, use Makefiles if not available
if [ "${CMAKE_GENERATOR}" == "" ]; then
    CMAKE_GENERATOR="Unix Makefiles"
fi

echo "Using CMAKE_MODULE_PATH set to ${CMAKE_MODULE_PATH} for the build"

# build all libraries
# must do separate runs for each configuration unlike Win/VC++

for CURPROJECT in myInterpolator; do
    for CURCONFIG in debug release; do

        CURPROJECT_BUILDDIR="${BUILDDIR}/$CURPROJECT/${CURCONFIG}"
        mkdir -p "${CURPROJECT_BUILDDIR}"
        cmake -G"${CMAKE_GENERATOR}" -DCMAKE_BUILD_TYPE=${CURCONFIG} -DCMAKE_INSTALL_PREFIX="${INSTALLDIR}" -B "${CURPROJECT_BUILDDIR}" -S "${BATDIR}"
        cmake --build "${CURPROJECT_BUILDDIR}" --target install --parallel
    done
done

#
# hint to set environment variable
#
echo ""
echo "Don't forget to include ${CMAKEMODULESDIR} in your CMAKE_MODULE_PATH variable if you haven't done so yet!"

export CMAKE_MODULE_PATH=${OLD_CMAKE_MODULE_PATH}
