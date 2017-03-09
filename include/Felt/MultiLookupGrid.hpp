#ifndef LOOKUPGRID_HPP_
#define LOOKUPGRID_HPP_

#include "LookupGridBase.hpp"

namespace felt
{

/**
 * Standard lookup grid.
 *
 * Holds a set of tracking lists storing grid positions, and a corresponding grid storing tuples of
 * list indices.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class MultiLookupGrid : public LookupGridBase< MultiLookupGrid<D, N> >
{
public:
	using ThisType = MultiLookupGrid<D, N>;
	using Base = LookupGridBase<ThisType>;
	using typename Base::PosArray;
public:
	using Base::LookupGridBase;
	using Base::list;
private:
	/// Make private, forcing use of overrides that specify list id.
	const PosArray& list() const {}
};


/**
 * Standard traits for all classes CRTP derived from LookupGridBase
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct DefaultMultiLookupGridTraits : DefaultGridTraits<VecDu<N>, D >
{
	/// Null index data type.
	using NULL_IDX_TYPE = VecDu<N>;
	/// Null index grid value in data array.
	static const NULL_IDX_TYPE NULL_IDX_DATA;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;
};


template <UINT D, UINT N>
const VecDu<N> DefaultMultiLookupGridTraits<D, N>::NULL_IDX_DATA = (
	VecDu<N>::Constant(std::numeric_limits<UINT>::max())
);


/**
 * Traits for MultiLookupGrid.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<MultiLookupGrid<D, N> > : DefaultMultiLookupGridTraits<D, N>
{
	using ThisType = MultiLookupGrid<D, N>;
	/// Set as eagerly constructed.
	static const Laziness IsLazy = Laziness::EAGER;
};


}
#endif /* MAPPEDGRID_HPP_ */
