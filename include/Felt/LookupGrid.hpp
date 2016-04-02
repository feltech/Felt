#ifndef LOOKUPGRID_HPP_
#define LOOKUPGRID_HPP_

#include "LookupGridBase.hpp"

namespace felt
{

/**
 * Base class for static lookup grid.
 *
 * @tparam Derived CRTP derived class.
 */
template <class Derived>
class StaticLookupGridBase : public LookupGridBase<StaticLookupGridBase<Derived>, false>
{
public:
	using Base = LookupGridBase<StaticLookupGridBase<Derived>, false>;
	using Base::LookupGridBase;
};

/**
 * Base class for lazy lookup grid.
 *
 * @tparam Derived CRTP derived class.
 */
template <class Derived>
class LazyLookupGridBase : public LookupGridBase<LazyLookupGridBase<Derived>, true>
{
public:
	using Base = LookupGridBase<LazyLookupGridBase<Derived>, true>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::PosArray;
	using Traits = GridTraits<LazyLookupGridBase<Derived> >;

	using Base::LookupGridBase;
	using Base::Base::is_active;
	using Base::is_active;

	/**
	 * @copydoc LazyGridBase::deactivate
	 *
	 * Additionally frees the tracking list(s).
	 */
	void deactivate()
	{
		Base::deactivate();
		for (PosArray& list : this->m_a_pos)
		{
			list.clear();
			list.shrink_to_fit();
		}
	}
};


template <class Derived, bool IsLazy> const UINT
LookupGridBase<Derived, IsLazy>::NULL_IDX = std::numeric_limits<UINT>::max();


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
class LookupGrid : public StaticLookupGridBase<LookupGrid<D, N> >
{
public:
	using ThisType = LookupGrid<D, N>;
	using Base = StaticLookupGridBase<ThisType>;
	using Base::StaticLookupGridBase;
};


/**
 * Lazy lookup grid - only initialised on activation, otherwise returns NULL_IDX when queried.
 *
 * @snippet test_MappedGrid.cpp LazyLookupGrid initialisation
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class LazyLookupGrid : public LazyLookupGridBase<LazyLookupGrid<D, N> >
{
public:
	using ThisType = LazyLookupGrid<D, N>;
	using Base = LazyLookupGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using Base::LazyLookupGridBase;
	using Base::reset;
};


/**
 * Traits for StaticLookupGridBase.
 *
 * Just forward the traits defined in subclasses.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
struct GridTraits< StaticLookupGridBase<Derived> > : GridTraits<Derived>
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
struct DefaultLookupGridTraits : DefaultGridTraits<VecDu<N>, D >
{
	/// Null index data type.
	using NULL_IDX_TYPE = VecDu<N>;
	/// Null index grid value in data array.
	static const NULL_IDX_TYPE NULL_IDX_DATA;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;
};

template <UINT D, UINT N>
const VecDu<N> DefaultLookupGridTraits<D, N>::NULL_IDX_DATA = (
	VecDu<N>::Constant(std::numeric_limits<UINT>::max())
);


/**
 * Traits for LookupGrid.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LookupGrid<D, N> > : DefaultLookupGridTraits<D, N>
{
	using ThisType = LookupGrid<D, N>;
};


/**
 * Traits for LazyLookupGrid.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LazyLookupGrid<D, N> > : DefaultLookupGridTraits<D, N>
{
	using ThisType = LazyLookupGrid<D, N>;
};


}
#endif /* MAPPEDGRID_HPP_ */
