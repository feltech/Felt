#ifndef SRC_TESTS_UTILS_HPP_
#define SRC_TESTS_UTILS_HPP_

#include <string>
#include <iomanip>

#include "catch.hpp"

#include "Felt/Grid.hpp"


namespace felt
{
	/// Utility: turn a number into a bit string.
	std::string stringify_bitmask(long mask, short length = 8);

	/// Utility: take a slice of a 3D grid and return a tabulated string.
	template <class Derived>
	std::string stringify_grid_slice(
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
			strGrid << std::endl;
			for (INT y = offset(axis_2); y < (INT)size(axis_2) + offset(axis_2);
				y++)
			{
				VecDi pos;
				if (axis_plane < pos.size())
					pos(axis_plane) = axis_plane_offset;
				pos(axis_1) = x;
				pos(axis_2) = y;
				strGrid << std::setw(5) << (FLOAT)grid(pos) << ",";
			}
		}
		strGrid << std::endl;
		return strGrid.str();
	}


	/**
	 * Copy of Catch::Detail::Approx to handle vector types.
	 */
	template <class VecType>
	class ApproxVecImpl
	{
	public:
		using ThisType = ApproxVecImpl<VecType>;

		explicit ApproxVecImpl (VecType value)
		:   m_epsilon( std::numeric_limits<float>::epsilon()*100 ),
			m_scale( 1.0 ),
			m_value( value )
		{}

		ApproxVecImpl( ThisType const& other )
		:   m_epsilon( other.m_epsilon ),
			m_scale( other.m_scale ),
			m_value( other.m_value )
		{}

		ThisType operator()( VecType value ) {
			Approx approx( value );
			approx.epsilon( m_epsilon );
			approx.scale( m_scale );
			return approx;
		}

		friend bool operator == ( VecType lhs, ThisType const& rhs ) {
			bool is_equal = true;
			for (UINT i = 0; i < lhs.size(); i++)
				is_equal &= fabs(lhs(i) - rhs.m_value(i)) < rhs.m_epsilon * (
					rhs.m_scale + std::max( fabs(lhs(i)), fabs(rhs.m_value(i)) )
				);
			return is_equal;
		}

		friend bool operator == ( ThisType const& lhs, VecType rhs ) {
			return operator==( rhs, lhs );
		}

		friend bool operator != ( VecType lhs, ThisType const& rhs ) {
			return !operator==( lhs, rhs );
		}

		friend bool operator != ( ThisType const& lhs, VecType rhs ) {
			return !operator==( rhs, lhs );
		}

		ThisType& epsilon( double newEpsilon ) {
			m_epsilon = newEpsilon;
			return *this;
		}

		ThisType& scale( double newScale ) {
			m_scale = newScale;
			return *this;
		}

		std::string toString() const {
			std::ostringstream oss;
			oss << "Approx( " << std::endl << m_value << std::endl << " )";
			return oss.str();
		}

	private:
		double m_epsilon;
		double m_scale;
		VecType m_value;
	};

	/**
	 * Wrap ApproxVec implementation in a function to enable template parameter deduction.
	 *
	 * @param value vector.
	 * @return ApproxVecImpl instance.
	 */
	template <class VecType>
	ApproxVecImpl<VecType> ApproxVec (VecType value)
	{
		return ApproxVecImpl<VecType>(value);
	}
}


namespace Catch
{
	template<class VecType>
	struct StringMaker< felt::ApproxVecImpl<VecType> >
	{
		static std::string convert( felt::ApproxVecImpl<VecType> const& value_ )
		{
			return value_.toString();
		}
	};
}

#endif /* SRC_TESTS_UTILS_HPP_ */
