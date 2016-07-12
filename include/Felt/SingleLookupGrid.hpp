#ifndef INCLUDE_FELT_SingleLOOKUPGrid_HPP_

#define INCLUDE_FELT_SingleLOOKUPGrid_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include <Felt/MultiLookupGrid.hpp>


namespace felt
{

/**
 * Similar to MultiLookupGrid but grid nodes store only a single list index.
 *
 * Useful in cases where grid nodes cannot be in more than one list.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <class Derived, Laziness IsLazy>
class SingleLookupGridBase : public LookupGridBase<SingleLookupGridBase<Derived, IsLazy>, IsLazy>
{
public:
	friend class LookupGridBase<SingleLookupGridBase<Derived, IsLazy>, IsLazy>;
	using ThisType = SingleLookupGridBase<Derived, IsLazy>;
	using Base = LookupGridBase<ThisType, IsLazy>;
	using Base::NULL_IDX;
	using typename Base::LeafType;
	using typename Base::VecDu;
	using typename Base::VecDi;
protected:
	using Base::m_a_pos;
public:
	using Base::LookupGridBase;

	/**
	 * Return true if position currently tracked for given list id, false
	 * otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if position is tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos_) const
	{
		return this->get(pos_) != Base::NULL_IDX;
	}

	/**
	 * Add position to tracking list and store index in tracking list in grid.
	 *
	 * Overloads LookupGridBase to place lookup index in first (and only) element of tuple of
	 * indices at the grid position.
	 *
	 * @param pos_ position in grid to track.
	 * @param arr_idx_ tracking list id to add position to.
	 * @return true if grid node set and position added to list, false if grid node was already set
	 * so position already in a list.
	 */
	bool add (const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		return Base::add(pos_, arr_idx_, 0u);
	}

	/**
	 * For given tracking list, set all lookup grid nodes to NULL index and
	 * clear the list.
	 *
	 * Overloads LookupGridBase to reset first (and only) element of tuple
	 * of indices for each grid position referenced in given tracking list.
	 *
	 * @param arr_idx_ tracking list id.
	 */
	void reset (const UINT arr_idx_ = 0)
	{
		Base::reset(arr_idx_, 0u);
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear all lists.
	 */
	void reset_all ()
	{
		for (UINT idx = 0; idx < m_a_pos.size(); idx++)
			Base::reset(idx, 0u);
	}

	/**
	 * Remove an element from a tracking list by index and set it's
	 * corresponding grid node to NULL index.
	 *
	 * Overloads LookupGridBase to set the first (and only) element of
	 * tuple of indices in the grid to NULL index.
	 *
	 * @param idx_ index in tracking list.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const UINT idx_, const UINT arr_idx_ = 0)
	{
		const VecDi& pos = this->list(arr_idx_)[idx_];
		Base::remove(idx_, pos, arr_idx_, 0);
	}

	/**
	 * Look up tracking list index in grid, remove from list and set grid
	 * node to NULL index.
	 *
	 * Overloads LookupGridBase to set the first (and only) element of
	 * tuple of indices in the grid to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		const UINT idx = idx_from_pos(pos_);
		Base::remove(idx, pos_, arr_idx_, 0);
	}

protected:

	/**
	 * Reset the entire grid to null indices.
	 */
	void clear ()
	{
		this->fill(NULL_IDX);
	}

	/**
	 * Get index in a tracking list from position.
	 *
	 * @param pos_ position in grid to find index data.
	 * @param arr_idx_ tracking list id.
	 * @return
	 */
	UINT& idx_from_pos(const VecDi& pos_, const UINT arr_idx_ = 0)
	{
		return this->get(pos_);
	}

	/**
	 * @copydoc idx_from_pos(const VecDi&,const UINT)
	 */
	const UINT idx_from_pos(const VecDi& pos_, const UINT arr_idx_) const
	{
		return this->get(pos_);
	}
};

/**
 * Base class for static lazy lookup grid with overlapping tracking lists.
 *
 * @tparam Derived CRTP derived class.
 */
template <class Derived>
class StaticSingleLookupGridBase : public SingleLookupGridBase <
	StaticSingleLookupGridBase<Derived>, Laziness::EAGER
> {
public:
	using Base = SingleLookupGridBase<StaticSingleLookupGridBase<Derived>, Laziness::EAGER>;
	using Base::SingleLookupGridBase;
};


/**
 * Base class for lazy lookup grid with non-overlapping tracking lists.
 *
 * @tparam Derived CRTP derived class.
 */
template <class Derived>
class LazySingleLookupGridBase : public SingleLookupGridBase <
	LazySingleLookupGridBase<Derived>, Laziness::LAZY
> {
public:
	using ThisType = LazySingleLookupGridBase<Derived>;
	using Base = SingleLookupGridBase<LazySingleLookupGridBase<Derived>, Laziness::LAZY>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using Base::Base::Base::is_active;
	using Base::is_active;
	using Traits = GridTraits< LazySingleLookupGridBase<Derived> >;

	LazySingleLookupGridBase() : Base()
	{}

	/**
	 * Construct lazy lookup grid, initialising the background value to NULL index.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 */
	LazySingleLookupGridBase(const VecDu& size_, const VecDi& offset_)
	{
		this->init(size_, offset_, Traits::NULL_IDX_DATA);
	}

	/**
	 * Copy constructor to augment base class.
	 *
	 * @param other_ grid to copy.
	 */
	LazySingleLookupGridBase(const ThisType& other_) : Base(other_)
	{
		this->m_background = other_.m_background;
	}

	/**
	 * Move constructor to augment base class.
	 *
	 * @param other_ grid to move.
	 */
	LazySingleLookupGridBase(ThisType&& other_) : Base(other_)
	{
		this->m_background = std::move(other_.m_background);
	}

	/**
	 * Copy assigment to augment base class.
	 *
	 * @param other_ grid to copy.
	 */
	ThisType& operator=(const ThisType& other_)
	{
		this->m_background = std::move(other_.m_background);
		return static_cast<ThisType&>(Base::operator=(other_));
	}

	/**
	 * Move assigment operator to augment base class.
	 *
	 * @param other_ grid to move.
	 */
	ThisType& operator=(ThisType&& other_)
	{
		this->m_background = std::move(other_.m_background);
		return static_cast<ThisType&>(Base::operator=(other_));
	}
};


/**
 * Standard traits for all classes CRTP derived from LookupGridBase
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N>
struct DefaultSingleLookupGridTraits : DefaultGridTraits<UINT, D >
{
	/// Null index grid value in data array.
	static const UINT NULL_IDX_DATA;
	/// Number of tracking lists taken from template parameter.
	static const UINT NumLists = N;
};

template <UINT D, UINT N>
const UINT DefaultSingleLookupGridTraits<D, N>::NULL_IDX_DATA = std::numeric_limits<UINT>::max();


/**
 * Concrete definition of standard SingleLookupGrid from SingleLookupGridBase.
 *
 * A simple stub to expose SingleLookupGridBase, which is also CRTP derived from elsewhere
 * - @see SingleTrackedPartitionedGrid.
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class SingleLookupGrid : public StaticSingleLookupGridBase<SingleLookupGrid<D, N> >
{
public:
	using ThisType = SingleLookupGrid<D, N>;
	using Base = StaticSingleLookupGridBase<ThisType>;
	using Base::StaticSingleLookupGridBase;
};


/**
 * Concrete definition of LazySingleLookupGrid from LazySingleLookupGridBase.
 *
 * A simple stub to expose LazySingleLookupGridBase, which is also CRTP derived from elsewhere
 * - @see LazySingleTrackedPartitionedGrid.
 *
 * @snippet test_MappedGrid.cpp LazySingleLookupGrid initialisation
 *
 * @tparam D the dimension of the grid.
 * @tparam N the number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class LazySingleLookupGrid : public LazySingleLookupGridBase<LazySingleLookupGrid<D, N> >
{
public:
	using ThisType = LazySingleLookupGrid<D, N>;
	using Base = LazySingleLookupGridBase<ThisType>;
	using Base::LazySingleLookupGridBase;
};


/**
 * Traits for SingleLookupGridBase.
 *
 * Just forward the traits defined for SingleLookupGridBase subclasses.
 */
template <class Derived, Laziness IsLazy>
struct GridTraits< SingleLookupGridBase<Derived, IsLazy> > : GridTraits<Derived>
{};


/**
 * Traits for LazySingleLookupGridBase.
 *
 * Just forward the traits defined for SingleLookupGridBase subclasses.
 */
template <class Derived>
struct GridTraits<LazySingleLookupGridBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits for LazySingleLookupGridBase.
 *
 * Just forward the traits defined for SingleLookupGridBase subclasses.
 */
template <class Derived>
struct GridTraits<StaticSingleLookupGridBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits for SingleLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<SingleLookupGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The derived type.
	using ThisType = SingleLookupGrid<D, N>;
};


/**
 * Traits for LazySingleLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LazySingleLookupGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The derived type.
	using ThisType = LazySingleLookupGrid<D, N>;
};

} // End namespace felt


#endif /* INCLUDE_FELT_SingleLOOKUPGrid_HPP_ */
