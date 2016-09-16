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
class EagerMultiLookupGrid : public EagerLookupGridBase<EagerMultiLookupGrid<D, N> >
{
public:
	using ThisType = EagerMultiLookupGrid<D, N>;
	using Base = EagerLookupGridBase<ThisType>;
	using typename Base::PosArray;
public:
	using Base::EagerLookupGridBase;
	using Base::list;
private:
	const PosArray& list() const {}
};


/**
 * Traits for StaticLookupGridBase.
 *
 * Just forward the traits defined in subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
struct GridTraits< EagerLookupGridBase<Derived> > : GridTraits<Derived>
{};



/**
 * Traits for LazyLookupGridBase.
 *
 * Just forward the traits defined in subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
struct GridTraits< LazyLookupGridBase<Derived> > : GridTraits<Derived>
{};


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
struct GridTraits<EagerMultiLookupGrid<D, N> > : DefaultMultiLookupGridTraits<D, N>
{
	using ThisType = EagerMultiLookupGrid<D, N>;
};


}
#endif /* MAPPEDGRID_HPP_ */
