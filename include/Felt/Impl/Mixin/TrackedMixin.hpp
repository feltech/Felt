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


template <class TDerived>
class Activate : protected Grid::Activate<TDerived>
{
private:
	/// Base class.
	using Base = Felt::Impl::Mixin::Grid::Activate<TDerived>;
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits::Leaf;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

protected:
	using Base::Activate;
	using Base::is_active;

	/**
	 * Allocate the internal data array and lookup grid.
	 */
	void activate()
	{
		Base::activate();
		pself->m_grid_lookup.activate();
	}

	/**
	 * Destroy the internal data array and lookup grid and change the background value
	 *
	 * @param background_ new background value to report when (now-inactive) grid is queried.
	 */
	void deactivate(Leaf background_)
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


template <class TDerived>
class Resize : protected Grid::Resize<TDerived>
{
private:
	using Base = Grid::Resize<TDerived>;
	/// Dimension of the grid.
	static const Dim t_dims = Traits<TDerived>::t_dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:

	void resize(const VecDi& size_, const VecDi& offset_)
	{
		Base::resize(size_, offset_);
		pself->m_grid_lookup.resize(size_, offset_);
	}
};


template <class TDerived>
class LookupInterface
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;

	/// Dimensions of the grid.
	static constexpr UINT t_dims = Traits::t_dims;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;

	using Leaf = typename Traits::Leaf;
	using Lookup = typename Traits::Lookup;

protected:
	Lookup m_grid_lookup;


protected:
	LookupInterface() {};

	LookupInterface(Lookup&& grid_lookup_)
		: m_grid_lookup(grid_lookup_)
	{}


	/**
	 * Serialisation hook for cereal library.
	 *
	 * @param ar
	 */
	template<class Archive>
	void serialize(Archive & ar)
	{
		ar(m_grid_lookup);
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
};


namespace SingleList
{

template <class TDerived>
class ByRef
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr UINT t_dims = Traits::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	using Leaf = typename Traits::Leaf;
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
	bool track(Leaf&& val_, const PosIdx pos_idx_)
	{
		pself->get(pos_idx_) = std::forward<Leaf>(val_);
		return pself->lookup().track(pos_idx_);
	}

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
	bool track(const Leaf& val_, const PosIdx pos_idx_)
	{
		pself->get(pos_idx_) = val_;
		return pself->lookup().track(pos_idx_);
	}
};


template <class TDerived>
class ByValue
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr UINT t_dims = Traits::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	using Leaf = typename Traits::Leaf;
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
	bool track(const Leaf val_, const PosIdx pos_idx_)
	{
		pself->set(pos_idx_, val_);
		return pself->lookup().track(pos_idx_);
	}
};


template <class TDerived>
class Reset
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr Dim t_dims = Traits::t_dims;
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
		for (const PosIdx pos_idx : pself->lookup().list())
			pself->get(pos_idx) = pself->background();
		pself->m_grid_lookup.reset();
	}
};

} // SingleList.


/**
 * Multiple tracking lists.
 */
namespace MultiList
{

template <class TDerived>
class ByRef
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr UINT t_dims = Traits::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	using Leaf = typename Traits::Leaf;
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
	bool track(Leaf&& val_, const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		pself->get(pos_idx_) = std::forward<Leaf>(val_);
		return pself->lookup().track(pos_idx_, list_idx_);
	}

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
	bool track(const Leaf& val_, const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		pself->get(pos_idx_) = val_;
		return pself->lookup().track(pos_idx_, list_idx_);
	}
};


template <class TDerived>
class ByValue
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr UINT t_dims = Traits::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	using Leaf = typename Traits::Leaf;
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
	bool track(const Leaf val_, const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		pself->set(pos_idx_, val_);
		return pself->lookup().track(pos_idx_, list_idx_);
	}
};


template <class TDerived>
class LookupInterface : protected Tracked::LookupInterface<TDerived>
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	/// Base class
	using Base = Tracked::LookupInterface<TDerived>;

	/// Dimensions of the grid.
	static constexpr UINT t_dims = Traits::t_dims;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;

	using Leaf = typename Traits::Leaf;
	using Lookup = typename Traits::Lookup;


protected:
	LookupInterface() {};


	LookupInterface(Lookup&& grid_lookup_) : Base(std::forward<Lookup>(grid_lookup_))
	{}

	/**
	 * Alias to access lookup grid's tracking lists.
	 *
	 * @param list_idx index of tracking list.
	 * @return reference to list.
	 */
	const PosIdxList& list(const TupleIdx list_idx) const
	{
		return this->m_grid_lookup.list(list_idx);
	}

	/**
	 * Alias to access lookup grid's tracking lists.
	 *
	 * @param list_idx index of tracking list.
	 * @return reference to list.
	 */
	PosIdxList& list(const TupleIdx list_idx)
	{
		return this->m_grid_lookup.list(list_idx);
	}
};


template <class TDerived>
class Reset
{
private:
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	static constexpr Dim t_dims = Traits::t_dims;
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;
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
		for (TupleIdx list_idx = 0; list_idx < t_num_lists; list_idx++)
		{
			for (const PosIdx pos_idx : pself->lookup().list(list_idx))
				pself->set(pos_idx, pself->background());
		}
		pself->m_grid_lookup.reset();
	}
};

} // MultiList.

} // Tracked
} // Mixin
} // Impl
} // Felt

#endif /* FELT_IMPL_TRACKED_HPP_ */
