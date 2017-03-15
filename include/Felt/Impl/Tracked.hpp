#ifndef FELT_IMPL_TRACKED_HPP_
#define FELT_IMPL_TRACKED_HPP_

#include <Felt/Util.hpp>
#include <Felt/Impl/Base.hpp>
#include <Felt/public/Grid.hpp>

namespace Felt
{
namespace Impl
{
namespace Tracked
{


template <class Derived>
class Activator : protected Grid::Activator<Derived>
{
protected:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Base class.
	using Base =  Felt::Impl::Grid::Activator<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;

protected:
	/**
	 * Allocate the internal data array and lookup grid.
	 */
	void activate()
	{
		Base::activate();
		nself->m_grid_lookup.activate();
	}
	/**
	 * Destroy the internal data array and lookup grid.
	 */
	void deactivate()
	{
		nself->m_data.clear();
		nself->m_data.shrink_to_fit();
		nself->m_grid_lookup.deactivate();
	}
};


template <class Derived>
class LazySingleByValue : private Base<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;

	/// Dimensions of the grid.
	static constexpr UINT Dims = TraitsType::Dims;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;

	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;

	using LeafType = typename TraitsType::LeafType;

protected:
	using Lookup = LazySingleLookupGrid<Dims, NumLists>;
	using PosArray = typename Lookup::PosArray;


protected:
	Lookup m_grid_lookup;

protected:
	LazySingleByValue(const VecDi& size, const VecDi& offset) :
		m_grid_lookup{size, offset}
	{}

	/**
	 * Destroy the internal and lookup data.
	 *
	 * @snippet test_Grid.cpp LazyGrid deactivation
	 */
	void deactivate()
	{
		nself->m_data.clear();
		nself->m_data.shrink_to_fit();
		m_grid_lookup.deactivate();
	}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	Lookup& lookup()
	{
		return m_grid_lookup;
	}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	const Lookup& lookup() const
	{
		return m_grid_lookup;
	}

	/**
	 * Get list of active grid points from lookup grid.
	 *
	 * @param list_idx_ tracking list id.
	 * @return the tracking list of active grid positions from the internal lookup grid.
	 */
	PosArray& list(const UINT list_idx_)
	{
		return m_grid_lookup.list(list_idx_);
	}

	/**
	 * Get list of active grid points from lookup grid.
	 *
	 * @param list_idx_ tracking list id.
	 * @return the tracking list of active grid positions from the internal lookup grid.
	 */
	const PosArray& list(const UINT list_idx_) const
	{
		return m_grid_lookup.list(list_idx_);
	}

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in lookup grid and position added to
	 * tracking list, false if grid node was already set so position already
	 * in a list.
	 */
	bool add(const VecDi& pos_, const LeafType val_, const UINT list_idx_)
	{
		nself->set(pos_, val_);
		return add(pos_, list_idx_);
	}

	/**
	 * Add a position to the lookup grid with given tracking list ID.
	 *
	 * @param pos_ position in the grid to add.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in lookup grid and position added to
	 * tracking list, false if grid node was already set so position already
	 * in a list.
	 */
	bool add(const VecDi& pos_, const UINT list_idx_)
	{
		return m_grid_lookup.add(pos_, list_idx_);
	}

	/**
	 * Set every active grid node (those referenced by lookup grid) to background value and reset
	 * the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices and it's tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param list_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const UINT list_idx_)
	{
		for (VecDi pos : m_grid_lookup.list(list_idx_))
			nself->set(pos, cself->m_background);
		m_grid_lookup.reset(list_idx_);
	}

	/**
	 * Get tracking list index from lookup grid, remove from list and set lookup
	 * grid node to NULL index.
	 *
	 * @param pos_ position in lookup grid.
	 * @param list_idx_ tracking list id.
	 */
	void remove (const VecDi& pos_, const UINT list_idx_)
	{
		m_grid_lookup.remove(pos_, list_idx_);
	}

	/**
	 * Check if position currently tracked for given list ID.
	 *
	 * @param pos_ position in grid to query.
	 * @param list_idx_ tracking list id to query.
	 * @return true if position currently tracked, false otherwise.
	 */
	const bool is_active (const VecDi& pos_, const UINT list_idx_) const
	{
		return m_grid_lookup.is_active(pos_, list_idx_);
	}
};

} // Lookup
} // Impl
} // Felt

#endif /* FELT_IMPL_TRACKED_HPP_ */
