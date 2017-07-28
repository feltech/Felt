#ifndef SRC_TESTS_UTILS_HPP_
#define SRC_TESTS_UTILS_HPP_

#include <string>
#include <iomanip>

#include "../catch.hpp"

#include "Felt/Grid.hpp"


namespace Felt
{
	/// Utility: take a slice of a 3D grid and return a tabulated string.
	template <class GridType>
	std::string stringify_grid_slice(
		const GridType& grid, UINT axis_plane = 2, INT axis_plane_offset = 0
	) {
		static constexpr UINT t_dims = Impl::Traits<GridType>::t_dims;
		using VecDi = typename Felt::VecDi<t_dims>;

		const VecDi& size = grid.size();
		const VecDi& offset = grid.offset();
		std::stringstream strGrid;
		UINT axis_1 = (axis_plane+1) % t_dims;
		UINT axis_2 = (axis_plane+2) % t_dims;
		INT z = axis_plane_offset;
		for (INT x = offset(axis_1); x < (INT)size(axis_1) + offset(axis_1); x++)
		{
			strGrid << std::endl;
			for (INT y = offset(axis_2); y < (INT)size(axis_2) + offset(axis_2); y++)
			{
				VecDi pos;
				if (axis_plane < pos.size())
					pos(axis_plane) = axis_plane_offset;
				pos(axis_1) = x;
				pos(axis_2) = y;
				strGrid << std::setw(5) << (FLOAT)grid.get(pos) << ",";
			}
		}
		strGrid << std::endl;
		return strGrid.str();
	}

	// https://wjngkoh.wordpress.com/2015/03/04/c-hash-function-for-eigen-matrix-and-vector/
	// Hash function for Eigen matrix and vector.
	// The code is from `hash_combine` function of the Boost library. See
	// http://www.boost.org/doc/libs/1_55_0/doc/html/hash/reference.html#boost.hash_combine .
	template<typename T>
	struct matrix_hash : std::unary_function<T, size_t> {
		std::size_t operator()(T const& matrix) const {
			// Note that it is oblivious to the storage order of Eigen matrix (column- or
			// row-major). It will give you the same hash value for two different matrices if they
			// are the transpose of each other in different storage order.
			size_t seed = 0;
			for (Eigen::Index i = 0; i < matrix.size(); ++i) {
				auto elem = *(matrix.data() + i);
				seed ^= std::hash<typename T::Scalar>()(elem) + 0x9e3779b9 + (seed << 6) +
					(seed >> 2);
			}
			return seed;
		}
	};
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
	struct StringMaker< Felt::ApproxVecImpl<VecType> >
	{
		static std::string convert( Felt::ApproxVecImpl<VecType> const& value_ )
		{
			return value_.toString();
		}
	};
}

#endif /* SRC_TESTS_UTILS_HPP_ */
