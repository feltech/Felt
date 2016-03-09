#ifndef LOOKUPGRID_HPP_
#define LOOKUPGRID_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include "Grid.hpp"

namespace felt
{

/** @defgroup LookupGrids
 *
 *  Bi-directional lookup between tracking list and grid locations - grids index lists and lists
 *  index grids.
 *  @{
 */


/**
 * Base class for a lookup grid.
 *
 * Array elements store grid positions and grid nodes store array indices.  Each grid node holds an
 * n-tuple of indices, and there are n arrays associated with this lookup grid.  An alternative
 * would be to use multiple grids each tracking a single list, but a multi-list structure allows
 * related lists to be spatially tracked in contiguous memory.
 *
 * @tparam Derived the CRTP derived class
 */
template <class Derived, bool IsLazy=false>
class LookupGridBase : public GridBase<LookupGridBase<Derived, IsLazy>, IsLazy>
{
public:
	using ThisType = LookupGridBase<Derived, IsLazy>;
	using Traits = GridTraits<Derived>;
	using DerivedType = typename Traits::ThisType;
	/// Type of data stored in grid nodes. For lookup grids this is an N-tuple of array indices.
	using LeafType = typename Traits::LeafType;
	/// GridBase base class.
	using Base = GridBase<ThisType, IsLazy>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::PosArray;
	/// An array index indicating a NULL index (nothing pointed to).
	static const UINT	NULL_IDX;
	/// Number of lists tracked in this grid.
	static const UINT 	NUM_LISTS = Traits::NumLists;

protected:
	/// N-tuple of lists of grid positions - the tracking lists.
	std::array<PosArray, NUM_LISTS>	m_a_pos;
	/// Mutex for use by other classes where multiple threads hold a reference to this grid.
	std::mutex						m_mutex;
public:

	/**
	 * Default (empty) destructor.
	 */
	~LookupGridBase() {}

	/**
	 * Explicitly defined default constructor.
	 */
	LookupGridBase () = default;

	/**
	 * Construct a lookup grid of given size and spatial offset.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 */
	LookupGridBase (const VecDu& size_, const VecDi& offset_ = VecDi::Zero()) : Base()
	{
		this->init(size_, offset_);
	}

	/**
	 * Copy constructor.
	 *
	 * Required since destructor defined and mutex is non-copyable.
	 *
	 * @param other_ grid to copy.
	 */
	LookupGridBase(const ThisType& other_)
	{
		m_a_pos = other_.m_a_pos;
		this->m_data = other_.m_data;
		this->m_offset = other_.m_offset;
		this->m_size = other_.m_size;
	}

	/**
	 * Move constructor,
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to move.
	 */
	LookupGridBase(ThisType&& other_)
	{
		m_a_pos = std::move(other_.m_a_pos);
		this->m_data = std::move(other_.m_data);
		this->m_offset = std::move(other_.m_offset);
		this->m_size = std::move(other_.m_size);
	}

	/**
	 * Copy assigment operator.
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to copy.
	 */
	ThisType& operator=(const ThisType& other_)
	{
		m_a_pos = other_.m_a_pos;
		this->m_data = other_.m_data;
		this->m_offset = other_.m_offset;
		this->m_size = other_.m_size;
		return *this;
	}

	/**
	 * Move assigment operator.
	 *
	 * Required since destructor defined and mutex is non-movable.
	 *
	 * @param other_ grid to move.
	 */
	ThisType& operator=(ThisType&& other_)
	{
		m_a_pos = std::move(other_.m_a_pos);
		this->m_data = std::move(other_.m_data);
		this->m_offset = std::move(other_.m_offset);
		this->m_size = std::move(other_.m_size);
		return *this;
	}

	/**
	 * Get mutex associated with this lookup grid - only used externally.
	 *
	 * Adding and removing elements from the tracking list/grid is not thread safe, so if multiple
	 * threads use this lookup grid, this mutex must be used.
	 *
	 * @return mutex for use by external objects.
	 */
	std::mutex& mutex ()
	{
		return m_mutex;
	}

	/**
	 * Reshape grid and fill with NULL indices.
	 *
	 * @param size_ new size of the grid.
	 */
	void size (const VecDu& size_)
	{
		Base::size(size_);
		self->clear();
	}

	using Base::size;

	/**
	 * Get tracking list by id.
	 *
	 * @param arr_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	inline PosArray& list (const UINT arr_idx_ = 0)
	{
		return m_a_pos[arr_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param arr_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	inline const PosArray& list (const UINT arr_idx_ = 0) const
	{
		return m_a_pos[arr_idx_];
	}

	/**
	 * Return true if position currently tracked for given list id, false
	 * otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @param arr_idx_ id of tracking list to test.
	 * @return true if grid position tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos_, const UINT arr_idx_ = 0) const
	{
		return cself->idx_from_pos(pos_, arr_idx_) != NULL_IDX;
	}

	/**
	 * Add position to tracking list and store index in tracking list in
	 * grid.
	 *
	 * Does nothing if grid node already set with an index in tracking list.
	 *
	 * @param pos_ position in grid to track.
	 * @param arr_idx_ id of tracking list to add to.
	 * @return true if grid node set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool add (const VecDi& pos, const UINT arr_idx_ = 0)
	{
		return add(pos, arr_idx_, arr_idx_);
	}

	/**
	 * Remove an element from a tracking list by index and set it's
	 * corresponding grid node to NULL index.
	 *
	 * @param idx_ index in tracking list.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const UINT idx, const UINT arr_idx_ = 0)
	{
		const VecDi& pos = this->list(arr_idx_)[idx];
		remove(idx, pos, arr_idx_, arr_idx_);
	}

	/**
	 * Look up tracking list index in grid, remove from list and set grid
	 * node to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 */
	void remove (const VecDi& pos, const UINT arr_idx_ = 0)
	{
		const UINT idx = this->get(pos)[arr_idx_];
		if (idx == NULL_IDX)
			return;
		remove(idx, pos, arr_idx_, arr_idx_);
	}

	/**
	 * Clear tracking list and reset every grid point to NULL index.
	 *
	 * @param arr_idx_ tracking list id to clear.
	 */
	void reset (const UINT arr_idx_ = 0)
	{
		reset(arr_idx_, arr_idx_);
	}

protected:

	/**
	 * Reset the entire grid to null indices.
	 */
	void clear ()
	{
		this->fill(LeafType::Constant(NULL_IDX));
	}

	/**
	 * Add pos to tracking list and set pos in grid to index in tracking list.
	 *
	 * If a grid node has a non-NULL index then does nothing.
	 *
	 * The arr_idx and lookup_idx can be different values - used in
	 * subclass overrides.  See LookupGrid.
	 *
	 * @param pos_ position in lookup grid.
	 * @param arr_idx_ tracking list id.
	 * @param lookup_idx_ lookup grid id.
	 * @return true if grid node was set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool add(const VecDi& pos_, const UINT arr_idx_, const UINT lookup_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		this->assert_pos_bounds(pos_, "add: ");
		#endif

		UINT& idx = self->idx_from_pos(pos_, lookup_idx_);
		// Do not allow duplicates.
		if (idx != NULL_IDX)
			return false;
		// idx is by reference, so this sets the value in grid.
		idx = this->list(arr_idx_).size();
		list(arr_idx_).push_back(pos_);
		return true;
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear all lists.
	 *
	 * @param arr_idx_ tracking list id.
	 */
	void reset_all ()
	{
		for (UINT idx = 0; idx < m_a_pos.size(); idx++)
			reset(idx, idx);
	}

	/**
	 * For given tracking list, set all lookup grid nodes to NULL index and
	 * clear the list.
	 *
	 * @param arr_idx_ tracking list id.
	 * @param lookup_idx_ lookup grid node tuple index.
	 */
	void reset(const UINT arr_idx_, const UINT lookup_idx_)
	{
		for (VecDi pos : m_a_pos[arr_idx_])
			self->idx_from_pos(pos, lookup_idx_) = NULL_IDX;
		list(arr_idx_).clear();
	}

	/**
	 * Get index in a tracking list from position.
	 *
	 * @param pos_ position in grid to find index data.
	 * @param arr_idx_ tracking list id.
	 * @return
	 */
	UINT& idx_from_pos(const VecDi& pos_, const UINT arr_idx_)
	{
		return this->get(pos_)[arr_idx_];
	}

	/**
	 * @copydoc idx_from_pos(const VecDi&,const UINT)
	 */
	const UINT idx_from_pos(const VecDi& pos_, const UINT arr_idx_) const
	{
		return this->get(pos_)[arr_idx_];
	}

	/**
	 * Remove pos at index idx in the array and set lookup at pos to NULL index.
	 *
	 * NOTE: idx passed by value since it changes indirectly via grid lookup.
	 *
	 * @param idx_ index of element in tracking list to remove.
	 * @param pos_ position in grid matching index in tracking list to remove.
	 * @param arr_idx_ tracking list id to remove element from.
	 * @param lookup_idx_ index in tuple stored at grid node that references element to remove.
	 */
	void remove(
		const UINT idx_, const VecDi& pos_, const UINT arr_idx_, const UINT lookup_idx_
	) {
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		this->assert_pos_bounds(pos_, "remove: ");
		#endif

		// Set index lookup to null value.
		self->idx_from_pos(pos_, lookup_idx_) = NULL_IDX;

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const UINT size = this->list(arr_idx_).size();
		if (idx_ < size - 1)
		{
			// Duplicate last element into this index.
			const VecDi& pos_last = this->list(arr_idx_)[size - 1];
			this->list(arr_idx_)[idx_] = pos_last;
			// Set the lookup grid to reference the new index in the array.
			self->idx_from_pos(pos_last, lookup_idx_) = idx_;
		}
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		list(arr_idx_).pop_back();
	}
};

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
	using Base::Base::is_active;
	using Base::is_active;
	using Traits = GridTraits<LazyLookupGridBase<Derived> >;

	/**
	 * Construct lazy lookup grid, initialising the background value to NULL index.
	 *
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 */
	LazyLookupGridBase(const VecDu& size_, const VecDi& offset_)
	{
		this->init(size_, offset_, Traits::NULL_IDX_DATA);
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
 * Traits for GridBase to understand LookupGridBase.
 *
 * Just forward the traits defined for LookupGridBase subclasses.
 *
 * @tparam Derived the CRTP derived class.
 * @tparam IsLazy true if the direct ancestor is a LazyGridBase.
 */
template <class Derived, bool IsLazy>
struct GridTraits< LookupGridBase<Derived, IsLazy> > : GridTraits<Derived>
{};


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
	/// Null index grid value in data array.
	static const VecDu<N> NULL_IDX_DATA;
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


/** @} */ // End group LookupGrid.

}
#endif /* MAPPEDGRID_HPP_ */
