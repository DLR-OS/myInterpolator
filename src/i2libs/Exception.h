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

#ifndef COMMON_EXCEPTION_H
#define COMMON_EXCEPTION_H


#include <string>
#include <exception>


namespace common {


    /**
     * CException is the basic exception class. It actually wraps the
     * std::exception class and transports it out of the std namespace
     * Additionally, it adds a message text.
     *
     * @author Dirk "jtk" Frommholz
     * @date April 22, 2008
     *
     */
    class CException {
    private:


        /**
         * The error message
         */
        std::string m_msg;


        /**
         * The location where the exception occured.
         */
        std::string m_location;


    public:


        /**
         * The constructor takes the string that further describes the error
         * and can be provided with the function's name that threw the 
         * CException (i.e. via the __func__ or _FILE_ macros).
         *
         * @param f_msg the exception message
         * @param f_location_p the location name (will be copied to the
         * CException instance)
         */
        CException(std::string  f_msg="", const char* f_location_p=0) {
            m_msg=f_msg;

            if (f_location_p) {
                m_location=*f_location_p;
            } else {
                m_location="(unknown function)";
            }
        }


        /**
         * What returns the error message.
         *
         * @return the error message
         */
        virtual const char* what() const throw() {
            return m_msg.c_str();
        }


        /**
         * Where returns the location where the exception originates from.
         * 
         * @return the location
         */
        virtual const char* where() const throw() {
            return m_location.c_str();
        }


         /**
          * The destructor.
          */
        virtual ~CException() { }
    }; 

} // namespace common


#endif
