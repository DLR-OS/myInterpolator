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

#include "Converter.h"

//****************************************************************************


unsigned common::CConverter::determineBase(const std::string& f_nr, 
std::string* f_strippedNr) {

    unsigned encBase=10;

    // strip sign if necessary
    std::string numStr=f_nr;
    bool isNegative=false;

    if (numStr.length()>0) {

        if ((f_nr[0]=='+') || (f_nr[0]=='-')) {
            numStr=f_nr.substr(1, f_nr.length()-1);
            isNegative=(f_nr[0]=='-');
        }

    	// set custom bases/remove encoded prefixes in numStr
    	switch (numStr[0]) {
    	    case '%':
    	    	encBase=2;
    	    	numStr=numStr.substr(1, numStr.length()-1);
    	    break;
    	    case '$':
    	    	encBase=16;
    	    	numStr=numStr.substr(1, numStr.length()-1);
    	    break;
    	    case '0':

                if (numStr.length()>1) {
                    if (numStr[1]=='x') {
                        encBase=16;
                    } else {
                        encBase=8;
                    }
                } else {
                   encBase=8;
                }
            
    	    	numStr=numStr.substr(1, numStr.length()-1);
            break;
    	}
    } else {
        // empty string -> 0
        numStr="0";
    }

    // output result
    if (f_strippedNr) {
        *f_strippedNr=isNegative ? '-'+numStr : numStr;
    }

    return encBase;
}


//****************************************************************************


unsigned common::CConverter::strToUInt(std::string f_nr, bool* f_error_p, 
unsigned f_base) {

    unsigned result=0;

    std::string strippedNr;
    const unsigned encBase=determineBase(f_nr, &strippedNr);

    // convert        
    char *errPtr_p=0;
    result=strtoul(strippedNr.c_str(), &errPtr_p, (f_base>0) ? f_base : 
    encBase);

    if (f_error_p) {
        *f_error_p=(errPtr_p!=0);
    }

    return result;
}



//****************************************************************************


double common::CConverter::strToDouble(const std::string& f_nr) {
    return strtod(f_nr.c_str(), 0);
}


//****************************************************************************


float common::CConverter::strToFloat(const std::string& f_nr) {
    return strtof(f_nr.c_str(), 0);
}


//****************************************************************************


std::string common::CConverter::uintToBin(unsigned f_nr, unsigned f_bits)
{
    std::string binStr="";

    if (f_nr==0) {
        binStr="0";
    }

    while (f_nr>0) {
    	if (f_nr%2) binStr='1'+binStr;
    	    else binStr='0'+binStr;
    	f_nr/=2;
    }

    while (binStr.length()<f_bits) {
        binStr='0'+binStr;
    }

    return binStr;
}


//****************************************************************************


std::string common::CConverter::doubleToStr(double f_d, std::streamsize 
f_prec, bool f_fixed) {

    // nicely looking version without exponentials and only as many
    // digits as necessary
    std::stringstream ss;
    ss << std::setprecision(f_prec) << std::setiosflags(std::ios_base::
    fixed) << f_d;

    if (f_fixed) {
        return ss.str();
    } else {
        // remove padding at end, special care for integers like 
        /// 1200.000 -> 1200

        std::string s=ss.str();
        std::size_t pos=s.find_last_not_of('0');
        if (pos != std::string::npos) {
            s.erase(pos+1);
        }

        pos=s.find_last_not_of('.');
        if(pos != std::string::npos) {
            s.erase(pos+1);
        }

        return s;
    }
}


//****************************************************************************

char* common::CConverter::upStrInPlace(char* f_str_p)
{
    unsigned c=0;
    while(f_str_p[c]!='\0') {
        f_str_p[c]=toupper(f_str_p[c]);
        ++c;
    }

    return f_str_p;
}


//****************************************************************************

std::string common::CConverter::upStr(const std::string& f_str) {

    std::string result(f_str);
    for (std::string::size_type c=0; c<result.length(); ++c) {
        result[c]=toupper(result[c]);
    }

    return result;
}


//****************************************************************************


std::string& common::CConverter::lowStrInPlace(std::string& f_str) {

    for (std::string::size_type c=0; c<f_str.length(); ++c) {
        f_str[c]=tolower(f_str[c]);
    }

    return f_str;
}


//****************************************************************************


std::string& common::CConverter::upStrInPlace(std::string& f_str) {

    for (std::string::size_type c=0; c<f_str.length(); ++c) {
        f_str[c]=toupper(f_str[c]);
    }

    return f_str;
}


//****************************************************************************


char* common::CConverter::lowStrInPlace(char* f_str_p) {
    unsigned c=0;
    while(f_str_p[c]!='\0') {
        f_str_p[c]=tolower(f_str_p[c]);
        ++c;
    }

    return f_str_p;
}


//****************************************************************************


std::string common::CConverter::lowStr(const std::string& f_str) {
    std::string result(f_str);
    for (std::string::size_type c=0; c<result.length(); ++c) {
        result[c]=tolower(result[c]);
    }

    return result;
}


//****************************************************************************


long int common::CConverter::strToInt(std::string f_nr, bool* f_error_p, 
unsigned f_base) {

    char* err_p;
    long int result=0;

    std::string strippedNr;
    const unsigned encBase=determineBase(f_nr, &strippedNr);

    // convert
    result=strtol(strippedNr.c_str(), &err_p, (f_base>0) ? f_base : encBase);
    if (f_error_p) {
        *f_error_p=(*err_p!=0);
    }

    return result;
}


//****************************************************************************


std::string common::CConverter::stripFileName(const std::string& f_fullPath)
{
    std::string::size_type lastSlashPosition=f_fullPath.find_last_of('/');

    if (lastSlashPosition==std::string::npos) {
        lastSlashPosition=f_fullPath.find_last_of('\\');
    }


    return (lastSlashPosition==std::string::npos) ?
    "" :
    f_fullPath.substr(0, lastSlashPosition+1);
}


//****************************************************************************


bool common::CConverter::isIntNumber(const std::string& f_str, bool 
f_allowSign) {

    // empty string is not an integer number
    if (f_str.length()==0) {
        return false;
    }

    // sign check
    std::string::size_type digitStart=0;

    if (f_allowSign) {
        if ( (f_str[0]!='+') && (f_str[0]!='-') &&
        !isdigit(f_str[0])) {
            return false;
        }

        // skip sign in check loop below if there is no sign
        // in a possibly signed number
        if (!isdigit(f_str[0])) {
            digitStart=1;
        }
    }

    for (std::string::size_type c=digitStart; c<f_str.length(); ++c) {
        if (!isdigit(f_str[c])) {
            return false;
        }        
    }

    // all tests passed, it's a number
    return true;
}


//****************************************************************************


std::string common::CConverter::padStr(const std::string& f_s, 
const char f_ch, unsigned f_targetLength, bool f_padFront) {

    std::string result(f_s);
    result.reserve(f_targetLength);

    while (result.length()<f_targetLength) {
        result=f_padFront ? f_ch+result : result+f_ch;
    }

    return result;
}


//****************************************************************************

