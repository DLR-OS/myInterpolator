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

#ifndef COMMON_STD_PROGRESS_H
#define COMMON_STD_PROGRESS_H

#include <iostream>
#include "Progress.h"


namespace common {


    /**
     * CStdProgress implements the progress indicator using iostream's 
     * std::cout object. It can be directly used for command line 
     * applications.
     *
     * @author Dirk "jtk" Frommholz
     * @date May 8, 2013
     */
    class CStdProgress: public CProgress {
        public:


            /**
             * The constructor.
             */
            CStdProgress():CProgress() { 
            }

    
            /**
             * ShowProgress shows and updates the progress indicator using
             * the iostream output interface.
             *
             * @param f_force show the progress no matter if it's time to
             */
            void showProgress(bool f_force=false) {
                if (needToVisualizeProgress() || f_force) {
                    std::cout << getOperation() << " ... " << 
                    static_cast<int>(100*getPercentage()) << 
                    "% (" << getElapsedTimeStr() << " passed, " << 
                    getRemainingTimeStr() << " remaining)" << std::endl;
                }
            }


            /**
             * The destructor.
             */
            ~CStdProgress() { 
            }
    };


} // namespace common


#endif
