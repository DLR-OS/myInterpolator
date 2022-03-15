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

#include "StopWatch.h"


/****************************************************************************/


common::CStopWatch::CStopWatch() {
    time(&m_baseTime);
    time(&m_stopTime);
}


/****************************************************************************/
    

void common::CStopWatch::start() {
    time(&m_baseTime);
}


/****************************************************************************/


void common::CStopWatch::stop() {    
    time(&m_stopTime);
}


/****************************************************************************/


double common::CStopWatch::getPassedTime() const {
    return difftime(m_stopTime, m_baseTime);
}


/****************************************************************************/


std::string common::CStopWatch::getPassedTimeStr() const {
    return CConverter::xToStr(getPassedTime());    
}


/****************************************************************************/


std::string common::CStopWatch::secsToHr(double f_secs) {

    const uint64_t passedSecs=static_cast<uint64_t>(fabs(f_secs));
    const uint64_t days=passedSecs/86400;
    const uint64_t hours=(passedSecs%86400)/3600;
    const uint64_t minutes=(passedSecs%3600)/60;
    const uint64_t seconds=(passedSecs)%60;

    std::string result=(f_secs<0) ? "-" : "";
    result.reserve(30);

    if (days>0) {        
        result+=(days<10) ? "0": "";
        result+=CConverter::xToStr(days)+"d ";
    }

    if ((days>0) || (hours>0)) {
        result+=(hours<10) ? "0" : "";
        result+=CConverter::xToStr(hours)+"h ";
    }

    if ((days>0) || (hours>0) || (minutes>0)) {
        result+=(minutes<10) ? "0" : "";
        result+=CConverter::xToStr(minutes)+"m ";
    }

    result+=((seconds>0) && (seconds<10)) ? "0" : "";
    result+=CConverter::xToStr(seconds)+"s";

    return result;
}


/****************************************************************************/


std::string common::CStopWatch::getCurrentTimeStrISO() const {
    time_t rawTime;        
    time(&rawTime);

    struct tm utcTime=*gmtime(&rawTime);
    std::string utcStr=
    // yyyy
    CConverter::padStr(CConverter::xToStr(1900+utcTime.tm_year), '0', 4, 
    true)+"-"+
    // mm
    CConverter::padStr(CConverter::xToStr(utcTime.tm_mon), '0', 2, true)+"-"+
    // dd
    CConverter::padStr(CConverter::xToStr(utcTime.tm_mday), '0', 2, true)+
    "T"+
    CConverter::padStr(CConverter::xToStr(utcTime.tm_hour), '0', 2, true)+
    ":"+
    CConverter::padStr(CConverter::xToStr(utcTime.tm_min), '0', 2, true)+
    ":"+
    CConverter::padStr(CConverter::xToStr(utcTime.tm_sec), '0', 2, true)+
    "+00:00";

    return utcStr;
}


/****************************************************************************/


std::string common::CStopWatch::getCurrentTimeStr() const {
    time_t currTime;
    time(&currTime);
    return ctime(&currTime);
}


/****************************************************************************/


std::pair<int, int> common::CStopWatch::getDayOfYear(int f_year, int f_month,
int f_day) const {
   
    // get current time, override if requested by the caller
    std::time_t t=std::time(0);
    std::tm* tm_p=std::localtime(&t);

    if (tm_p) {

        std::tm tm=*tm_p;
       
        tm.tm_year=(f_year==std::numeric_limits<int>::min()) ? tm.tm_year : 
        f_year-1900;
        tm.tm_mon=(f_month==std::numeric_limits<int>::min()) ? tm.tm_mon : 
        f_month-1;
        tm.tm_mday=(f_day==std::numeric_limits<int>::min()) ? tm.tm_mday : 
        f_day;
    
        std::mktime(&tm); 
        std::pair<int, int> doy(tm.tm_yday, tm.tm_year+1900);
        return doy;
    }

    // error retrieving the time
    return std::pair<int, int>(1, 1900);
}


/****************************************************************************/


bool common::CStopWatch::getCurrentTime(std::tm& f_time) {
    std::time_t t=std::time(0);
    std::tm* tm_p=std::localtime(&t);

    if (tm_p) {
        f_time=*tm_p;
        return true;
    }
    
    return false;
}


/****************************************************************************/


common::CStopWatch::~CStopWatch() {
}


/****************************************************************************/
