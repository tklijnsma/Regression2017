// Minimalistic GRBMaker for semiparametric regression
// Based on work by Jean-Baptiste Sauvan and Josh Bendavid

#ifndef UTILITIES_H
#define UTILITIES_H

#include <string>
#include <iostream>
#include <vector>
#include <sstream>

/**
 * @brief  Separates a string into tokens.
 * @param  str : the string to tokenize
 * @param  delimiters : delimites tokens
 * @return tokens
 */
void tokenize(const std::string& str,
        std::vector<std::string>& tokens,
        const std::string& delimiter = " ");

std::string intToString(int n);

void findAndReplace(std::string& sInput, std::string sFind, std::string sReplace );

void strip(std::string& sInput);


/**
 * @brief  Converts a string into base types.
 * @param  s : string to convert	 
 * @return t : converted string
 */
    template <class T>
bool fromString(T& t, 
        const std::string& s, 
        std::ios_base& (*f)(std::ios_base&) = std::dec)
{
    std::istringstream iss(s);
    return !(iss >> f >> t).fail();
}




#endif


