#ifndef FELT_IMPL_LOOKUP_HPP
#define FELT_IMPL_LOOKUP_HPP

#include <Felt/public/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/Grid.hpp>

namespace Felt
{

/// Value to store in lookup grid when a node doesn't reference any list.
static constexpr UINT	NULL_IDX = std::numeric_limits<UINT>::max();

namespace Impl
{
namespace Mixin
{
namespace Lookup
{

template <class Derived>
class Base
{
protected:
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Integer vector.
	using VecDi = Felt::VecDi<Dims>;
	/// List of position vectors.
	using PosArray = std::vector<VecDi>;
};


template <class Derived>
class Activator : protected Impl::Mixin::Grid::Activator<Derived>
{
protected:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// CRTP derived class.
	using Base =  Felt::Impl::Mixin::Grid::Activator<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;

protected:
	using Base::activate;

	/**
	 * Destroy the internal data array.
	 */
	void deactivate()
	{
		nself->m_data.clear();
		nself->m_data.shrink_to_fit();
		for (UINT list_idx = 0; list_idx < NumLists; list_idx++)
		{
			nself->m_a_list_pos[list_idx].clear();
			nself->m_a_list_pos[list_idx].shrink_to_fit();
		}
	}
};


template <class Derived>
class Simple : protected Base<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Base class
	using BaseType = Base<Derived>;

	/// Integer vector
	using typename BaseType::VecDi;
	/// List of position vectors.
	using typename BaseType::PosArray;

private:
	/// List of position vectors, each of which have a corresponding grid node storing it's index.
	PosArray	m_list_pos;

protected:


protected:

	/**
	 * Get tracking list.
	 *
	 * @return tracking list.
	 */
	const PosArray& list() const
	{
		return m_list_pos;;
	}
	/**
	 * Get tracking list.
	 *
	 * @return tracking list.
	 */
	PosArray& list()
	{
		return m_list_pos;;
	}

	/**
	 * Return true if position currently tracked for given list id, false otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_active (const VecDi& pos_) const
	{
		return cself->get(pos_) != NULL_IDX;
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
	bool add(const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "add: ");
		#endif

		UINT idx = cself->get(pos_);
		// Do not allow duplicates.
		if (idx != NULL_IDX)
		{
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			if (!m_list_pos.size() > idx)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pos_) << " detected as a duplicate, since " <<
					idx << " is not " << NULL_IDX << ", but no list is that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
			#endif
			return false;
		}

		nself->set(pos_, m_list_pos.size());
		m_list_pos.push_back(pos_);

		return true;
	}


	/**
	 * Remove pos from the array and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to remove.
	 */
	void remove(const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "remove: ");
		#endif
		// By reference, so we don't have to query again.
		UINT& idx_at_pos = nself->ref(pos_);

		if (idx_at_pos == NULL_IDX)
			return;

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const UINT last_idx = m_list_pos.size() - 1;
		if (idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const VecDi& pos_last = m_list_pos[last_idx];
			m_list_pos[idx_at_pos] = m_list_pos[last_idx];
			// Set the lookup grid to reference the new index in the array.
			nself->set(pos_last, idx_at_pos);
		}
		// NULL out the value in the grid now that we're done with it.
		idx_at_pos = NULL_IDX;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		m_list_pos.pop_back();
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear the list.
	 */
	void reset()
	{
		for (const VecDi& pos : m_list_pos)
			nself->set(pos, NULL_IDX);
		m_list_pos.clear();
	}
};


template <class Derived>
class Single : protected Base<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Base class
	using BaseType = Base<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Integer vector
	using typename BaseType::VecDi;
	/// Number of tracking lists.
	static const UINT NumLists = TraitsType::NumLists;
protected:
	/// List of position vectors.
	using typename BaseType::PosArray;
protected:
	/// N-tuple of lists of grid positions - the tracking lists.
	std::array<PosArray, NumLists>	m_a_list_pos;
protected:
	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	PosArray& list (const UINT list_idx_)
	{
		return m_a_list_pos[list_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	const PosArray& list (const UINT list_idx_) const
	{
		return m_a_list_pos[list_idx_];
	}

	/**
	 * Return true if position currently tracked for given list id, false otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_active (const VecDi& pos_) const
	{
		return cself->get(pos_) != NULL_IDX;
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
	bool add(const VecDi& pos_, const UINT list_idx_)
	{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "add: ");
#endif

		const UINT idx = cself->get(pos_);
		// Do not allow duplicates.
		if (idx != NULL_IDX)
		{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			bool found = false;

			for (UINT list_idx = 0; list_idx < m_a_list_pos.size() && !found; list_idx++)
				if (m_a_list_pos[list_idx].size() > idx)
					found = true;

			if (!found)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pos_) << " detected as a duplicate, since " <<
				idx << " is not " << NULL_IDX << ", but no list is that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
#endif
			return false;
		}
		PosArray& list_to_update = m_a_list_pos[list_idx_];
		nself->set(pos_, list_to_update.size());
		list_to_update.push_back(pos_);
		return true;
	}

	/**
	 * Remove pos from tracking list and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to remove.
	 * @param list_idx_ tracking list id to remove element from.
	 */
	void remove(const VecDi& pos_, const UINT list_idx_)
	{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "remove: ");
#endif

		// Set index lookup to null value.
		UINT& idx_at_pos = nself->ref(pos_);

		if (idx_at_pos == NULL_IDX)
			return;

		// Get a reference to the tracking list that we now need to update.
		PosArray& list_to_update = m_a_list_pos[list_idx_];

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const UINT last_idx = list_to_update.size() - 1;
		if (idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const VecDi& pos_last = list_to_update[last_idx];
			list_to_update[idx_at_pos] = pos_last;
			// Set the lookup grid to reference the new index in the array.
			nself->set(pos_last, idx_at_pos);
		}
		// NULL out the old value in the grid now that we're done with it.
		idx_at_pos = NULL_IDX;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		list_to_update.pop_back();
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear the list.
	 *
	 * @param list_idx_ tracking list id.
	 */
	void reset(const UINT list_idx_)
	{
		for (const VecDi& pos : m_a_list_pos[list_idx_])
			nself->set(pos, NULL_IDX);
		m_a_list_pos[list_idx_].clear();
	}
};


template <class Derived>
class Multi : private Base<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Base class
	using BaseType = Base<Derived>;
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;

	/// Integer vector
	using typename BaseType::VecDi;
	/// List of position vectors.
	using typename BaseType::PosArray;

	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;
	/// Tuple of indices stored at grid nodes.
	using IndexTuple = typename TraitsType::LeafType;

protected:
	static const IndexTuple NULL_IDX_TUPLE;

private:

	/// List of position vectors, each of which have a corresponding grid node storing it's index.
	/// N-tuple of lists of grid positions - the tracking lists.
	std::array<PosArray, NumLists>	m_a_list_pos;

protected:

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	PosArray& list (const UINT list_idx_)
	{
		return m_a_list_pos[list_idx_];
	}

	/**
	 * Get tracking list by id.
	 *
	 * @param list_idx_ id of tracking list to get.
	 * @return tracking list at given index.
	 */
	const PosArray& list (const UINT list_idx_) const
	{
		return m_a_list_pos[list_idx_];
	}

	/**
	 * Return true if position currently tracked for given list id, false otherwise.
	 *
	 * @param pos_ position in grid to query.
	 * @return true if grid position tracked, false otherwise.
	 */
	bool is_active (const VecDi& pos_, const UINT list_idx_) const
	{
		return cself->get(pos_)[list_idx_] != NULL_IDX;
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
	bool add(const VecDi& pos_, const UINT list_idx_)
	{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "add: ");
#endif
		UINT& idx = nself->get(pos_)[list_idx_];
		// Do not allow duplicates.
		if (idx != NULL_IDX)
		{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			bool found = false;

			for (UINT list_idx = 0; list_idx < m_a_list_pos.size() && !found; list_idx++)
				if (m_a_list_pos[list_idx].size() > idx)
					found = true;

			if (!found)
			{
				std::stringstream sstr;
				sstr << "Position " << Felt::format(pos_) << " detected as a duplicate, since " <<
				idx << " is not " << NULL_IDX << ", but no list is that big";
				std::string str = sstr.str();
				throw std::domain_error(str);
			}
#endif
			return false;
		}
		PosArray& list_to_update = m_a_list_pos[list_idx_];
		// Update value in grid at appropriate tuple index.
		idx = list_to_update.size();
		list_to_update.push_back(pos_);
		return true;
	}

	/**
	 * Remove pos from tracking list and set lookup at pos to NULL index.
	 *
	 * @param pos_ position in grid matching index in tracking list to remove.
	 * @param list_idx_ tracking list id to remove element from.
	 */
	void remove(const VecDi& pos_, const UINT list_idx_)
	{
#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "remove: ");
#endif

		// Set index lookup to null value.
		UINT& idx_at_pos = nself->get(pos_)[list_idx_];

		if (idx_at_pos == NULL_IDX)
			return;

		// Get a reference to the tracking list that we now need to update.
		PosArray& list_to_update = m_a_list_pos[list_idx_];

		// If this is not the last remaining position in the array, then
		// we must move the last position to this position and update the
		// lookup grid.
		const UINT last_idx = list_to_update.size() - 1;
		if (idx_at_pos < last_idx)
		{
			// Duplicate last element into this index.
			const VecDi& pos_last = list_to_update[last_idx];
			list_to_update[idx_at_pos] = pos_last;
			// Set the lookup grid to reference the new index in the array.
			nself->get(pos_last)[list_idx_] = idx_at_pos;
		}
		// NULL out the old value in the grid now that we're done with it.
		idx_at_pos = NULL_IDX;
		// Remove the last element in the array (which is at this point
		// either the last remaining element or a duplicate).
		list_to_update.pop_back();
	}

	/**
	 * Set all lookup grid nodes to NULL index and clear the list.
	 *
	 * @param list_idx_ tracking list id.
	 */
	void reset(const UINT list_idx_)
	{
		for (const VecDi& pos : m_a_list_pos[list_idx_])
			nself->get(pos) = NULL_IDX_TUPLE;
		m_a_list_pos[list_idx_].clear();
	}
};

template <class Derived>
const typename Traits<Derived>::LeafType
Multi<Derived>::NULL_IDX_TUPLE = IndexTuple::Constant(NULL_IDX);

} // Lookup
} // Mixin
} // Impl
} // Felt

# endif // FELT_IMPL_LOOKUP_HPP
