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

#ifndef COMMON_STOPWATCH_H
#define COMMON_STOPWATCH_H

#include <string>
#include <time.h>
#include <cmath>
#include <limits>
#include "Converter.h"


namespace common {


    /**
     * CStopWatch returns the current time and measures time differences.
     *
     * @author Dirk "jtk" Frommholz
     * @date October 9th, 2009
     */
    class CStopWatch {
        private:


            /**
             * The base time.
             */
            time_t m_baseTime;


            /**
             * The last time stopped.
             */
            time_t m_stopTime;


        public:


            /**
             * The constructor.
             */
            CStopWatch();
    

            /**
             * Start starts the stopwatch and sets the time base.
             */
            void start();


            /**
             * Stop stops the stopwatch. This function can be called multiple
             * times.
             */
            void stop();


            /**
             * GetPassedTime returns the time in seconds that has 
             * passed between the last call of start() and the last call 
             * of stop().
             *
             * @return the time between start() and the last call to stop()
             * as a double
             */
            double getPassedTime() const;


            /**
             * GetPassedTimeStr returns the string with the time that has 
             * passed between the last call of start() and the last call of 
             * stop(). Precision is one second.
             *
             * @return the time between start() and the last call to stop()
             * as a string
             */
            std::string getPassedTimeStr() const;


            /**
             * SecsToHr converts f_secs to a human-readable text string
             * in dd:hh:mm:ss format.
             *
             * @param f_secs the seconds to be converted
             * @return f_secs in human-readable form
             */
            static std::string secsToHr(double f_secs);


            /**
             * GetPassedTimeHrStr returns a string with the time that has
             * passed in a human readable dd-hh-mm-ss format
             *
             * @return the time passed in a human readable format
             */
            std::string getPassedTimeStrHr() const {
                return secsToHr(getPassedTime());
            }


            /**
             * GetCurrentTimeStr returns the current system time.
             *
             * @return the current system time as a string
             */
            std::string getCurrentTimeStr() const;


            /**
             * GetCurrentTimeStrISO returns a string containing the
             * current system time (UTC) in ISO8601 format.
             *
             * @return the current UTC time as of ISO8601
             */
            std::string getCurrentTimeStrISO() const;


            /**
             * GetCurrentTime returns the current local time.
             *
             * @param f_time the current time on success
             * @return true if the time can be retrieved, false 
             * if not; in this case f_time is undefined
             */
            static bool getCurrentTime(std::tm& f_time);


            /**
             * GetDayOfYear returns the day of the year for the date
             * given or, if called with a partial parameter set, for
             * the current year, month and/or day in local time.
             *
             * @param f_year the year to get the day of year for 
             * @param f_month the month to get the day of year for [1..12]
             * @param f_day the day to get the day of year for [1..31]
             * @return the current day (.first) of the current year
             * (.second)
             */
            std::pair<int, int> getDayOfYear(int f_year=std::numeric_limits<
            int>::min(), int f_month=std::numeric_limits<int>::min(), 
            int f_day=std::numeric_limits<int>::min()) const;


            /**
             * The destructor.
             */
            ~CStopWatch();

    };


} // namespace common

#endif
