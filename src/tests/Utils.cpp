#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <string>

// Utility: turn a number into a bit string.
std::string felt::stringifyBitmask(long mask, short length)
{
	std::string str;
	for (unsigned bitIdx = 0; bitIdx < length; bitIdx++)
		str += std::to_string(1 & (mask >> (length-1-bitIdx)));
	return str;
}
