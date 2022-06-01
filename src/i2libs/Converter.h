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

#ifndef COMMON_CONVERTER_H
#define COMMON_CONVERTER_H

#include <string>
#include <cstdlib>
#include <cctype>
#include <sstream>
#include <iomanip>
#include <ctype.h>
#include <cstdio>

namespace common {


    /**
     * CConverter provides useful functions for data type conversion
     * and formatting.
     *
     * @author Dirk 'jtk' Frommholz
     * @date April 1, 2007
     */
    class CConverter {
        public:


            /**
             * DetermineBase determines the base of an (unsigned) integer
             * stored as text in f_nr by evaluating any encoded base
             * identifiers like 0x/$ (hex), 0 (octal), % (binary).
             * If no base encoding is found base 10 will be returned.
             *
             * @param f_nr the number string
             * @param f_strippedNr f_nr stripped by the base encoding
             * @return the base
             */
            static unsigned determineBase(const std::string& f_nr, std::
            string* f_strippedNr);


            /**
             * StrToUInt converts a string into an unsigned integer. The hex
             * (0x, $) octal (0) and binary (%) prefixes are obeyed if the
             * f_base remains set to 0 (aka auto mode). In auto-base mode
             * and if there is no prefix encoded base 10 will be assumed.
             *
             * @param f_nr the number string to be converted
             * @param f_error_p the optional error flag, will be set on 
             * conversion errors
             * @param f_base the base to be used; if zero the base encoded 
             * in f_nr will be used (or 10 if there is no base encoded in
             * f_nr); if positive the f_base will be enforced and any base 
             * prefixes inside f_nr will be removed and ignored
             * @return the integer representation of f_nr if there where no
             * errors
             */
            static unsigned strToUInt(std::string f_nr, bool* f_error_p=0, 
            unsigned f_base=0);


            /**
             * PadStr pads the given string f_s by putting the character 
             * f_ch in the front or back of f_s until the target length 
             * of the string is reached. If the length of f_s equals 
             * or exceeds the target length on entry this function 
             * does nothing.
             * 
             * @param f_s the string to be padded
             * @param f_ch the character to be placed in front or back of 
             * f_s until the target length is reached
             * @param f_targetLength the target length of f_s after padding
             * @param f_padFront set true to add copies of f_ch at the front
             * of f_s and false to add copies of f_ch at the end of f_s
             */
            static std::string padStr(const std::string& f_s, 
            const char f_ch, unsigned f_targetLength, bool f_padFront);
            

            /**
             * UintToBin converts an unsigned integer to a binary string.
             * The number of bits desired can be specified too; in case
             * the conversion result has less than f_bits leading zero
             * bits will be added. However, if the conversion result does
             * not fit into f_bits the f_bits value will be violated.           
             *
             * @param f_nr the number to be converted
             * @param f_bits the number of bits to return
             */
            static std::string uintToBin(unsigned f_nr, unsigned f_bits=0);


            /**
             * XToStr converts all common data types to a string. Optionally,
             * a leading character may be specified. If so, the length of the
             * result must be also given.
             *
             * @param f_a the value to be converted to a string
             * @param f_leadingCharacter the optional leading character
             * @param f_width the with of the resulting string if a leading
             * character is given (should be larger than the value a in its
             * string representation)
             */
            template<typename T>
            static std::string xToStr(T f_a, char f_leadingCharacter=0, 
            unsigned f_width=0) {
                std::stringstream ss;

                if(f_width==0) {
                    ss << f_a;
                } else {
                    ss << std::setfill(f_leadingCharacter);
                    ss << std::setw(f_width) << f_a;
                }

                return ss.str();
            }


            /**
             * StrToX converts string s to another data type which has to be
             * specified by the template parameter.
             *
             * @param f_str the string to be converted
             * @return the converted string
             */
            template<typename T>
            static T strToX(std::string f_str) {
                std::stringstream ss;
                ss << f_str;
                T result;
                ss >> result;
                return result;
            }

    
            /**
             * DoubleToStr converts floating-point numbers to strings. The
             * precision can be specified (usually this is the maximum
             * meaningful number of digits).
             *
             * @param f_d the number to be converted to a string
             * @param f_prec the precision (meaningful digits)
             * @param f_fixed force fixed notation, i.e. trailing zeros;
             * if false up to f_fixed digits will be printed as necessary
             * before the scientific notation is used
             * @return the converted number d as a string
             */
            static std::string doubleToStr(double f_d, std::streamsize 
            f_prec=12, bool f_fixed=false);


            /**
             * StrToDouble convertes a string to a double-precision floating
             * point number.
             *
             * @param f_nr the number string
             * @return the floating point conversion result
             */
            static double strToDouble(const std::string& f_nr);


            /**
             * StrTofloat convertes a string to a single-precision floating
             * point number.
             *
             * @param f_nr the number string
             * @return the floating point conversion result
             */
            static float strToFloat(const std::string& f_nr);


            /**
             * StrToInt converts a string to an signed integer with respect
             * to the given base. If no base is given, 10 (decimal) is 
             * assumed. See strToUInt.
             *
             * @param f_nr the number string that is to be converted
             * @param f_error_p holds true if there was a conversion error, 
             * and in this case the result is undefined
             * @param f_base the base that shall be used if the number does
             * not contain a base identifier itself
             * @return the integer number that is contained in f_nr
             */
            static long int strToInt(std::string f_nr, bool* f_error_p=0, 
            unsigned f_base=10);

            /**
             * UpStr converts the letters in the given C string to
             * uppercase in-place using the toupper() function on all 
             * characters of the C string. For convenience, a pointer to the
             * converted C string is returned.
             *
             * @param f_str_p the C string to be converted to uppercase in-place
             * @return the converted string
             */
            static char* upStrInPlace(char* f_str_p);


            /**
             * LowStr converts the letters in the given C string to
             * lowercase in-place using the tolower() function on all 
             * characters of the C string. For convenience, a pointer to the
             * converted C string is returned.
             *
             * @param f_str_p the C string to be converted to lowercase 
             * in-place
             * @return the converted string
             */
            static char* lowStrInPlace(char* f_str_p);


            /**
             * UpStr converts the letters in the given C++ string to
             * uppercase and returns the result.
             *
             * @param f_str the C++ string to be converted to uppercase 
             * @return the converted string
             */
            static std::string upStr(const std::string& f_str);


            /**
             * UpStrInPlace converts the letters in the given C++ string to
             * uppercase in place and returns a reference to the result
             * for convenience.
             *
             * @param f_str the C++ string to be converted to uppercase 
             * in-place
             * @return a reference to the converted string
             */
            static std::string& upStrInPlace(std::string& f_str);


            /**
             * LowStr converts the letters in the given C++ string to
             * lowercase and returns the result.
             *
             * @param f_str the C++ string to be converted to lowercase
             * @return the converted string
             */
            static std::string lowStr(const std::string& f_str);


            /**
             * LowStrInPlace converts the letters in the given C++ string to
             * lowercase in place and returns a reference to the result
             * for convenience.
             *
             * @param f_str the C++ string to be converted to lowercase 
             * in-place
             * @return a reference to the converted string
             */
            static std::string& lowStrInPlace(std::string& f_str);


            /**
             * StripFileName returns the path portion of the Unix-style or 
             * Windows style full path stored in f_fullPath including the 
             * trailing (back)slash. For example, /mnt/media/1.txt is returned 
             * as /mnt/media/ by this function.
             *
             * @return the path without the file name of f_fullPath
             */
            static std::string stripFileName(const std::string& f_fullPath);


            /**
             * StringBufferLength returns the first position of the zero 
             * character in the array f_buffer_p. If no C string terminator is 
             * found the function returns f_maxLen which is to be provided by 
             * the caller.
             *
             * @param f_buffer_p the character buffer
             * @param f_maxLen the length to be returned if f_buffer_p does     
             * not contain the C string terminator (zero character)
             */
            template<typename CharT>
            static std::string::size_type stringBufferLength(CharT* 
            f_buffer_p, std::string::size_type f_maxLen) {

                for (std::string::size_type p=0; p<f_maxLen; ++p) {
                    if (f_buffer_p[p]=='\0') {
                        return p;
                    }
                }

                return f_maxLen;
            }


            /**
             * RemoveLeadTrailSpaces removes leading and trailing spaces from
             * the string f_s.
             *
             * @param f_s the string to be cleaned from leading and trailing 
             * spaces
             */
            static void removeLeadTrailSpaces(std::string& f_s) {
        
                std::string::size_type firstNonSpace=std::string::npos;
                std::string::size_type lastNonSpace=std::string::npos;

                // find first non-space    
                for (std::string::size_type c=0; c<f_s.length(); ++c) {
                    if (!isspace(static_cast<unsigned char>(f_s[c]))) {
                        firstNonSpace=c;
                        break;
                    }
                }

                // cut spaces from front
                if (firstNonSpace!=std::string::npos) {
                    if (firstNonSpace>0) {
                        f_s=f_s.erase(0, firstNonSpace);
                    }
                } else {
                    // have whitespaces only, return with empty string
                    f_s="";
                    return;
                }

                // now, find last non-space on modified string
                for (int64_t c=static_cast<int64_t>(f_s.length()-1); c>=0; --c) {
                    if (!isspace(static_cast<unsigned char>(f_s[c]))) {
                        lastNonSpace=c;
                        break;
                    }
                }

                // cut space from the end
                if (lastNonSpace!=std::string::npos) {
                    if (lastNonSpace+1<f_s.length()) {
                        f_s=f_s.erase(lastNonSpace+1);
                    }
                }
            }


            /**
             * IsNumber returns true if the string passed consists of a
             * sequence of digits only, optionally preceeded by a sign.
             * Non-digit characters including whitespaces will make the 
             * function return false. The empty string is not considered 
             * an integer number, hence the function will return false in 
             * this case.
             *
             * @param f_str the string to examine
             * @param f_allowSign true if f_str may have a sign as the
             * first character
             * @return true if f_str contains digits only (optionally preceeded
             * by a sign), false if not
             */
            static bool isIntNumber(const std::string& f_str, bool 
            f_allowSign=false);

    };

} // namespace

#endif
