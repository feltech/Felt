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
	template <class Derived>
	std::string stringifyGridSlice(
		const GridBase<Derived>& grid, UINT axis_plane = 2, INT axis_plane_offset = 0
	) {
		using GridType = GridBase<Derived>;
		using VecDu = typename GridBase<Derived>::VecDu;
		using VecDi = typename GridBase<Derived>::VecDi;

		const VecDu& size = grid.size();
		const VecDi& offset = grid.offset();
		std::stringstream strGrid;
		UINT axis_1 = (axis_plane+1) % GridType::Dims;
		UINT axis_2 = (axis_plane+2) % GridType::Dims;
		INT z = axis_plane_offset;
		for (INT x = offset(axis_1); x < (INT)size(axis_1) + offset(axis_1);
			x++)
		{
			strGrid << std::endl << "|";
			for (INT y = offset(axis_2); y < (INT)size(axis_2) + offset(axis_2);
				y++)
			{
				VecDi pos;
				if (axis_plane < pos.size())
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
