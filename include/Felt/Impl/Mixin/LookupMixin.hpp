#ifndef FELT_IMPL_LOOKUP_HPP
#define FELT_IMPL_LOOKUP_HPP

#include <Felt/Impl/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Lookup
{

template <class TDerived>
class Base
{
protected:
	/// Dimension of the grid.
	static const Dim t_dims = Traits<TDerived>::t_dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<t_dims>;
};


namespace SingleList
{

template <class TDerived>
class Activate : protected Grid::Activate<TDerived>
{
protected:
	/// CRTP derived class.
	using Base =  Grid::Activate<TDerived>;
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	// Base class methods.
	using Base::Activate;
	using Base::activate;
	using Base::is_active;

	/**
	 * Destroy the internal data array.
	 */
	void deactivate()
	{
		pself->m_data.clear();
		pself->m_data.shrink_to_fit();
		pself->m_list_pos_idxs.clear();
		pself->m_list_pos_idxs.shrink_to_fit();
	}
};


template <class TDerived>
class SingleIdx : protected Base<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Base class
	using Base = Lookup::Base<TDerived>;
	/// Integer vector
	using typename Base::VecDi;
protected:
	/// List of position vectors, each of which have a corresponding grid node storing it's index.
	PosIdxList	m_list_pos_idxs;

	/**
	 * Get tracking list.
	 *
	 * @return tracking list.
	 */
	const PosIdxList& list() const
	{
		return m_list_pos_idxs;;
	}
	/**
	 * Get tracking list.
	 *
	 * @return tracking list.
	 */
	PosIdxList& list()
	{
		return m_list_pos_idxs;;
	}

	/**
	 * Return true if position currently tracked for given list id, false otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_tracked (const PosIdx pos_idx_) const
	{
		return pself->get(pos_idx_) != Felt::null_idx;
	}

	/**
	* Add pos to tracking list and set pos in grid to index in tracking list.
	*
	* If a grid node has a non-NULL index then does nothing.
	*
	* @param pos_ position in lookup grid.
	* @return true if grid node was set and position added to list,
	* false if grid node was already set so position already in a list.
	*/
	bool track(const PosIdx pos_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_idx_, "track: ");
		#endif

		ListIdx idx = pself->get(pos_idx_);
		// Do not allow duplicates.
		if (idx != Felt::null_idx)
		{
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			if (m_list_pos_idxs.size()-1 < idx)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pself->index(pos_idx_)) <<
					" detected as a duplicate, since " << idx << " is not " << Felt::null_idx <<
					", but the list is not that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
			#endif
			return false;
		}

		pself->set(pos_idx_, ListIdx(m_list_pos_idxs.size()));
		m_list_pos_idxs.push_back(pos_idx_);

		return true;
	}


	/**
	 * Remove pos from the array and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to untrack.
	 */
	void untrack(const PosIdx pos_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_idx_bounds(pos_idx_, "untrack: ");
		#endif
		// By reference, so we don't have to query again.
		ListIdx& idx_at_pos = pself->ref(pos_idx_);

		if (idx_at_pos == Felt::null_idx)
			return;

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const ListIdx last_idx = ListIdx(m_list_pos_idxs.size()) - 1;
		if (idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const PosIdx pos_idx_last = m_list_pos_idxs[last_idx];
			m_list_pos_idxs[idx_at_pos] = m_list_pos_idxs[last_idx];
			// Set the lookup grid to reference the new index in the array.
			pself->set(pos_idx_last, idx_at_pos);
		}
		// NULL out the value in the grid now that we're done with it.
		idx_at_pos = Felt::null_idx;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		m_list_pos_idxs.pop_back();
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear the list.
	 */
	void reset()
	{
		for (const PosIdx pos_idx : m_list_pos_idxs)
			pself->set(pos_idx, Felt::null_idx);
		m_list_pos_idxs.clear();
	}
};

} // Single.


namespace MultiList
{

template <class TDerived>
class Activate : protected Grid::Activate<TDerived>
{
protected:
	/// CRTP derived class.
	using Base =  Grid::Activate<TDerived>;
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	// Base class methods.
	using Base::Activate;
	using Base::activate;
	using Base::is_active;

	/**
	 * Destroy the internal data array.
	 */
	void deactivate()
	{
		pself->m_data.clear();
		pself->m_data.shrink_to_fit();
		for (TupleIdx list_idx = 0; list_idx < t_num_lists; list_idx++)
		{
			pself->m_a_list_pos_idxs[list_idx].clear();
			pself->m_a_list_pos_idxs[list_idx].shrink_to_fit();
		}
	}
};

template <class TDerived>
class SingleIdx : protected Base<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Base class
	using Base = Lookup::Base<TDerived>;
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;
	/// Integer vector
	using typename Base::VecDi;
	/// Number of tracking lists.
	static const TupleIdx t_num_lists = Traits::t_num_lists;
protected:
	/// N-tuple of lists of grid positions - the tracking lists.
	Tuple<PosIdxList, t_num_lists>	m_a_list_pos_idxs;
protected:

	/**
	 * Serialisation hook for cereal library.
	 *
	 * @param ar
	 */
	template<class Archive>
	void serialize(Archive & ar)
	{
		ar(m_a_list_pos_idxs);
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	PosIdxList& list (const TupleIdx list_idx_)
	{
		return m_a_list_pos_idxs[list_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	const PosIdxList& list (const TupleIdx list_idx_) const
	{
		return m_a_list_pos_idxs[list_idx_];
	}

	/**
	 * Return true if position currently tracked for given list id, false otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_tracked (const VecDi& pos_) const
	{
		return pself->get(pos_) != Felt::null_idx;
	}

	/**
	 * Add pos to tracking list and set pos in grid to index in tracking list.
	 *
	 * If a grid node has a non-NULL index then does nothing.
	 *
	 * @param pos_ position in lookup grid.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node was set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool track(const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		#ifdef FELT_DEBUG_ENABLED
		pself->assert_pos_bounds(pos_idx_, "track: ");
		#endif

		const ListIdx idx = pself->get(pos_idx_);
		// Do not allow duplicates.
		if (idx != Felt::null_idx)
		{
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			bool found = false;

			for (UINT list_idx = 0; list_idx < m_a_list_pos_idxs.size() && !found; list_idx++)
				if (m_a_list_pos_idxs[list_idx].size() >= idx)
					found = true;

			if (!found)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pself->index(pos_idx_)) <<
					" detected as a duplicate, since " << idx << " is not " << Felt::null_idx <<
					", but no list is that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
			#endif
			return false;
		}
		PosIdxList& list_to_update = m_a_list_pos_idxs[list_idx_];
		pself->set(pos_idx_, ListIdx(list_to_update.size()));
		list_to_update.push_back(pos_idx_);
		return true;
	}

	/**
	 * Remove pos from tracking list and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to untrack.
	 * @param list_idx_ tracking list id to untrack element from.
	 */
	void untrack(const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		#ifdef FELT_DEBUG_ENABLED
		pself->assert_pos_bounds(pos_idx_, "untrack: ");
		#endif

		// Get reference to list index stored at given position in grid.
		ListIdx& list_idx_at_pos = pself->ref(pos_idx_);

		if (list_idx_at_pos == Felt::null_idx)
			return;

		// Get a reference to the tracking list that we now need to update.
		PosIdxList& list_to_update = m_a_list_pos_idxs[list_idx_];

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const ListIdx last_idx = ListIdx(list_to_update.size()) - 1;
		if (list_idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const PosIdx pos_idx_last = list_to_update[last_idx];
			list_to_update[list_idx_at_pos] = pos_idx_last;
			// Set the lookup grid to reference the new index in the array.
			pself->set(pos_idx_last, list_idx_at_pos);
		}
		// NULL out the old value in the grid now that we're done with it.
		list_idx_at_pos = Felt::null_idx;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		list_to_update.pop_back();
	}


	/**
	 * Set all lookup grid nodes to NULL index and clear the lists.
	 */
	void reset()
	{
		for (TupleIdx list_idx = 0; list_idx < t_num_lists; list_idx++)
		{
			PosIdxList& list_pos_idxs = m_a_list_pos_idxs[list_idx];
			for (const PosIdx pos_idx : list_pos_idxs)
				pself->set(pos_idx, Felt::null_idx);
			list_pos_idxs.clear();
		}
	}
};


template <class TDerived>
class MultiIdx : private Base<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Base class
	using Base = Lookup::Base<TDerived>;
	/// Traits of derived class.
	using Traits = Impl::Traits<TDerived>;

	/// Integer vector
	using typename Base::VecDi;

	/// Number of tracking lists.
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;
	/// Tuple of indices stored at grid nodes.
	using IndexTuple = typename Traits::Leaf;

protected:
	static const IndexTuple s_null_idxs;
private:
	/// List of position vectors, each of which have a corresponding grid node storing it's index.
	/// N-tuple of lists of grid positions - the tracking lists.
	Tuple<PosIdxList, t_num_lists>	m_a_list_pos_idxs;
protected:

	/**
	 * Serialisation hook for cereal library.
	 *
	 * @param ar
	 */
	template<class Archive>
	void serialize(Archive & ar)
	{
		ar(m_a_list_pos_idxs);
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	PosIdxList& list (const TupleIdx list_idx_)
	{
		return m_a_list_pos_idxs[list_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	const PosIdxList& list (const TupleIdx list_idx_) const
	{
		return m_a_list_pos_idxs[list_idx_];
	}

	/**
	 * Check whether a given position is tracked by given tracking list.
	 *
	 * @param pos_idx_ position in grid to query.
	 * @param list_idx_ tracking list id.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_tracked (const PosIdx pos_idx_, const TupleIdx list_idx_) const
	{
		return pself->get(pos_idx_)[list_idx_] != Felt::null_idx;
	}

	/**
	 * Check whether a given position is tracked by any tracking list.
	 *
	 * @param pos_idx_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_tracked (const PosIdx pos_idx_) const
	{
		return pself->get(pos_idx_) != s_null_idxs;
	}

	/**
	 * Add pos to tracking list and set pos in grid to index in tracking list.
	 *
	 * If a grid node has a non-NULL index then does nothing.
	 *
	 * @param pos_ position in lookup grid.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node was set and position added to list,
	 * false if grid node was already set so position already in a list.
	 */
	bool track(const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_idx_, "track: ");
		#endif
		ListIdx& idx = pself->get(pos_idx_)[list_idx_];
		// Do not allow duplicates.
		if (idx != Felt::null_idx)
		{
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			bool found = false;

			for (UINT list_idx = 0; list_idx < m_a_list_pos_idxs.size() && !found; list_idx++)
				if (m_a_list_pos_idxs[list_idx].size() >= idx)
					found = true;

			if (!found)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pself->index(pos_idx_)) <<
					" detected as a duplicate, since " << idx << " is not " << Felt::null_idx <<
					", but no list is that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
			#endif
			return false;
		}
		PosIdxList& list_to_update = m_a_list_pos_idxs[list_idx_];
		// Update value in grid at appropriate tuple index.
		idx = ListIdx(list_to_update.size());
		list_to_update.push_back(pos_idx_);
		return true;
	}

	/**
	 * Remove pos from tracking list and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to untrack.
	 * @param list_idx_ tracking list id to untrack element from.
	 */
	void untrack(const PosIdx pos_idx_, const TupleIdx list_idx_)
	{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_idx_, "untrack: ");
#endif

		// Set index lookup to null value.
		ListIdx& idx_at_pos = pself->get(pos_idx_)[list_idx_];

		if (idx_at_pos == Felt::null_idx)
			return;

		// Get a reference to the tracking list that we now need to update.
		PosIdxList& list_to_update = m_a_list_pos_idxs[list_idx_];

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const ListIdx last_idx = ListIdx(list_to_update.size()) - 1;
		if (idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const PosIdx pos_idx_last = list_to_update[last_idx];
			list_to_update[idx_at_pos] = pos_idx_last;
			// Set the lookup grid to reference the new index in the array.
			pself->get(pos_idx_last)[list_idx_] = idx_at_pos;
		}
		// NULL out the old value in the grid now that we're done with it.
		idx_at_pos = Felt::null_idx;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		list_to_update.pop_back();
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear the lists.
	 */
	void reset()
	{
		for (TupleIdx list_idx = 0; list_idx < t_num_lists; list_idx++)
		{
			PosIdxList& list_pos_idxs = m_a_list_pos_idxs[list_idx];
			for (const PosIdx pos_idx : list_pos_idxs)
				pself->get(pos_idx)[list_idx] = Felt::null_idx;
			list_pos_idxs.clear();
		}
	}

};

template <class TDerived>
const typename Traits<TDerived>::Leaf
MultiIdx<TDerived>::s_null_idxs = IndexTuple::Constant(Felt::null_idx);
} // Multi

} // Lookup
} // Mixin
} // Impl
} // Felt

# endif // FELT_IMPL_LOOKUP_HPP
