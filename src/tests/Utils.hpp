#ifndef SRC_TESTS_UTILS_HPP_
#define SRC_TESTS_UTILS_HPP_

#include <string>
#include <iomanip>

#include "Felt/Grid.hpp"


namespace felt
{
	/// Utility: turn a vector into a string.
	template <class T>
	std::string stringifyVector(const T& p, const INT& prec = 3)
	{
		std::stringstream str;
		str << "(";
		for (UINT i = 0; i < p.size(); i++)
		{
			if (i != 0)
				str << ", ";
			str << std::setprecision(prec) << p(i);
		}
		str << ")";
		return str.str();
	}

	/// Utility: turn a number into a bit string.
	std::string stringifyBitmask(long mask, short length = 8);

	/// Utility: take a slice of a 3D grid and return a tabulated string.
	// Utility: take a slice of a 3D grid and return a tabulated string.
	template <typename T>
	std::string stringifyGridSlice(
		const Grid<T,3>& grid, UINT axis_plane = 2, INT axis_plane_offset = 0
	) {
		const Vec3u& dims = grid.dims();
		const Vec3i& offset = grid.offset();
		std::stringstream strGrid;
		UINT axis_1 = (axis_plane+1)%3;
		UINT axis_2 = (axis_plane+2)%3;
		INT z = axis_plane_offset;
		for (INT x = offset(axis_1); x < (INT)dims(axis_1) + offset(axis_1);
			x++)
		{
			strGrid << std::endl << "|";
			for (INT y = offset(axis_2); y < (INT)dims(axis_2) + offset(axis_2);
				y++)
			{
				Vec3i pos;
				pos(axis_plane) = axis_plane_offset;
				pos(axis_1) = x;
				pos(axis_2) = y;
				strGrid << std::setw(5) << (FLOAT)grid(pos) << " |";
			}
		}
		strGrid << std::endl;
		return strGrid.str();
	}
}

#endif /* SRC_TESTS_UTILS_HPP_ */