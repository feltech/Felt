#ifndef FELT_IMPL_TRACKED_HPP_
#define FELT_IMPL_TRACKED_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Util.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Tracked
{


template <class Derived>
class Activator : protected Grid::Activator<Derived>
{
private:
	/// Base class.
	using Base = Felt::Impl::Mixin::Grid::Activator<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Type of data to store in grid nodes.
	using LeafType = typename TraitsType::LeafType;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;

protected:
	using Base::Activator;
	using Base::is_active;

	/**
	 * Allocate the internal data array and lookup grid.
	 */
	void activate()
	{
		Base::activate();
		pself->m_grid_lookup.activate();
	}

	void deactivate(LeafType background_)
	{
		this->m_background = background_;
		deactivate();
	}

	/**
	 * Destroy the internal data array and lookup grid.
	 */
	void deactivate()
	{
		pself->m_data.clear();
		pself->m_data.shrink_to_fit();
		pself->m_grid_lookup.deactivate();
	}
};


template <class Derived>
class Resize : protected Grid::Resize<Derived>
{
private:
	using Base = Grid::Resize<Derived>;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	using Base::Resize;

	void resize(const VecDi& size_, const VecDi& offset_)
	{
		Base::resize(size_, offset_);
		pself->m_grid_lookup.resize(size_, offset_);
	}
};


template <class Derived>
class ByRef
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	static constexpr UINT Dims = TraitsType::Dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;
	using LeafType = typename TraitsType::LeafType;
protected:
	/**
	 * Set value in grid at given position and track position to lookup grid.
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
	bool track(const LeafType& val_, const PosIdx pos_idx_, const UINT list_idx_)
	{
		pself->get(pos_idx_) = val_;
		return pself->lookup().track(pos_idx_, list_idx_);
	}

	/**
	 * @copydoc track(const VecDi&, const LeafType&, const UINT)
	 */
	bool track(LeafType&& val_, const PosIdx pos_idx_, const UINT list_idx_)
	{
		pself->get(pos_idx_) = val_;
		return pself->lookup().track(pos_idx_, list_idx_);
	}
};


template <class Derived>
class ByValue
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	static constexpr UINT Dims = TraitsType::Dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;
	using LeafType = typename TraitsType::LeafType;
protected:
	/**
	 * Set value in grid at given position and track position to lookup grid.
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
	bool track(const LeafType val_, const PosIdx pos_idx_, const UINT list_idx_)
	{
		pself->set(pos_idx_, val_);
		return pself->lookup().track(pos_idx_, list_idx_);
	}
};


template <class Derived>
class LookupInterface
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;

	/// Dimensions of the grid.
	static constexpr UINT Dims = TraitsType::Dims;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;

	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;

	using LeafType = typename TraitsType::LeafType;
	using LookupType = typename TraitsType::LookupType;

protected:
	using PosArray = typename LookupType::PosArray;

protected:
	LookupType m_grid_lookup;


protected:
	LookupInterface(LookupType&& grid_lookup_)
		: m_grid_lookup(grid_lookup_)
	{}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	LookupType& lookup()
	{
		return m_grid_lookup;
	}

	/**
	 * Get lookup grid.
	 *
	 * @return the internal lookup grid tracking active grid positions.
	 */
	const LookupType& lookup() const
	{
		return m_grid_lookup;
	}
};


template <class Derived>
class Resetter
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	static constexpr UINT Dims = TraitsType::Dims;
	static constexpr ListIdx NumLists = TraitsType::NumLists;
	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;
protected:

	/**
	 * Set every active grid node (those referenced by lookup grid) to background value and reset
	 * the lookup grid.
	 *
	 * Lookup grid will then be full of NULL indices and it's tracking list(s) will be empty.
	 *
	 * @param list_idx_ tracking list id to cycle over and clear.
	 */
	void reset()
	{
		for (ListIdx list_idx = 0; list_idx < NumLists; list_idx++)
		{
			for (const PosIdx pos_idx : pself->m_grid_lookup.list(list_idx))
				pself->set(pos_idx, pself->m_background);
		}
		pself->m_grid_lookup.reset();
	}
};


} // Tracked
} // Mixin
} // Impl
} // Felt

#endif /* FELT_IMPL_TRACKED_HPP_ */
