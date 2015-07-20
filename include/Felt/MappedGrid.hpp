#ifndef MAPPEDGRID_HPP_
#define MAPPEDGRID_HPP_

#include <array>
#include <sstream>
#include <mutex>
#include "Grid.hpp"

namespace felt
{
	/**
	 * Base class for a lookup grid, where array elements store grid positions
	 * and grid nodes store array indices.
	 *
	 * @tparam D the dimension of the grid.
	 * @tparam N the number of tracking lists to use.
	 * @tparam I the type representing an index in the grid (e.g. 3-tuple if
	 * 3x tracking lists are used).
	 * @tparam V the type returned when querying the grid.
	 */
	template <UINT D, UINT N, typename I, typename V>
	class LookupGridBase : public GridBase<I, D, V>
	{
	public:
		typedef I								Idx_t;
		typedef V								Val_t;
		typedef LookupGridBase<D, N, I, V>		ThisType;
		typedef GridBase<Idx_t, D, Val_t>		Base;
		typedef typename Base::VecDu			VecDu;
		typedef typename Base::VecDi			VecDi;
		typedef typename Base::PosArray			PosArray;
	protected:
		std::array<PosArray, N>	m_a_pos;
		std::mutex				m_mutex;
	public:

		static const Idx_t				 			NULL_IDX_TUPLE;
		static const UINT				 			NULL_IDX;

		static inline const UINT num_lists()
		{
			return N;
		}

		virtual ~LookupGridBase() {}

		LookupGridBase() : Base()
		{}

		/**
		 * Copy constructor.
		 *
		 * Required since destructor defined and mutex is non-copyable.
		 *
		 * @param other
		 */
		LookupGridBase(const ThisType& other)
		{
			m_a_pos = other.m_a_pos;
			this->m_data = other.m_data;
			this->m_offset = other.m_offset;
			this->m_dims = other.m_dims;
		}

		/**
		 * Move constructor,
		 *
		 * Required since destructor defined and mutex is non-movable.
		 *
		 * @param other
		 */
		LookupGridBase(ThisType&& other)
		{
			m_a_pos = std::move(other.m_a_pos);
			this->m_data = std::move(other.m_data);
			this->m_offset = std::move(other.m_offset);
			this->m_dims = std::move(other.m_dims);
		}

		/**
		 * Copy assigment operator.
		 *
		 * Required since destructor defined and mutex is non-movable.
		 *
		 * @param other
		 */
		ThisType& operator=(const ThisType& other)
		{
			m_a_pos = other.m_a_pos;
			this->m_data = other.m_data;
			this->m_offset = other.m_offset;
			this->m_dims = other.m_dims;
			return *this;
		}

		/**
		 * Move assigment operator.
		 *
		 * Required since destructor defined and mutex is non-movable.
		 *
		 * @param other
		 */
		ThisType& operator=(ThisType&& other)
		{
			m_a_pos = std::move(other.m_a_pos);
			this->m_data = std::move(other.m_data);
			this->m_offset = std::move(other.m_offset);
			this->m_dims = std::move(other.m_dims);
			return *this;
		}

		/**
		 * Construct a lookup grid with given dims and offset.
		 *
		 * @param dims
		 * @param offset
		 */
		LookupGridBase(const VecDu& dims, const VecDi& offset)
		: Base()
		{
			this->init(dims, offset);
		}

		/**
		 * Get mutex.
		 *
		 * @return
		 */
		std::mutex& mutex ()
		{
			return m_mutex;
		}

		/**
		 * Reshape grid and fill with NULL indices.
		 *
		 * @param dims_new
		 */
		virtual void dims (const VecDu& dims_new)
		{
			Base::dims(dims_new);
			this->fill(NULL_IDX_TUPLE);
		}

		VecDu dims ()
		{
			return Base::dims();
		}

		/**
		 * Get tracking list by id.
		 *
		 * @param arr_idx
		 * @return
		 */
		virtual inline PosArray& list (const UINT& arr_idx = 0)
		{
			return m_a_pos[arr_idx];
		}
		/**
		 * Get tracking list by id.
		 *
		 * @param arr_idx
		 * @return
		 */
		virtual inline const PosArray& list (const UINT& arr_idx = 0) const
		{
			return m_a_pos[arr_idx];
		}

		/**
		 * Return true if position currently tracked for given list id, false
		 * otherwise.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return
		 */
		const bool is_active (const VecDi& pos, const UINT& arr_idx = 0) const
		{
			return this->get(pos)[arr_idx] != NULL_IDX;
		}

		/**
		 * Add position to tracking list and store index in tracking list in
		 * grid.
		 *
		 * Does nothing if grid node already set with an index in tracking list.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return true if grid node set and position added to list,
		 * false if grid node was already set so position already in a list.
		 */
		virtual bool add (const VecDi& pos, const UINT& arr_idx = 0)
		{
			return add(pos, arr_idx, arr_idx);
		}

		/**
		 * Clear tracking list and reset every grid point to NULL index.
		 *
		 * @param arr_idx
		 */
		virtual void reset (const UINT& arr_idx = 0)
		{
			reset(arr_idx, arr_idx);
		}


		/**
		 * Remove an element from a tracking list by index and set it's
		 * corresponding grid node to NULL index.
		 *
		 * @param idx index in tracking list.
		 * @param arr_idx tracking list id.
		 */
		virtual void remove (const UINT& idx, const UINT& arr_idx = 0)
		{
			const VecDi& pos = this->list(arr_idx)[idx];
			remove(idx, pos, arr_idx, arr_idx);
		}

		/**
		 * Look up tracking list index in grid, remove from list and set grid
		 * node to NULL index.
		 *
		 * @param pos position in lookup grid.
		 * @param arr_idx tracking list id.
		 */
		virtual void remove (const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = this->get_internal(pos)[arr_idx];
			if (idx == NULL_IDX)
				return;
			remove(idx, pos, arr_idx, arr_idx);
		}

	protected:

		/**
		 * Add pos to tracking list and set pos in grid to index in tracking
		 * list.
		 *
		 * If a grid node has a non-NULL index then does nothing.
		 *
		 * The arr_idx and lookup_idx can be different values - used in
		 * subclass overrides.  See LookupGrid.
		 *
		 * @param pos position in lookup grid.
		 * @param arr_idx tracking list id.
		 * @param lookup_idx lookup grid id.
		 * @return true if grid node was set and position added to list,
		 * false if grid node was already set so position already in a list.
		 */
		bool add(
			const VecDi& pos, const UINT& arr_idx, const UINT& lookup_idx
		) {
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			this->assert_pos_bounds(pos, "add: ");
			#endif

			UINT& idx = this->get_internal(pos)[lookup_idx];
			// Do not allow duplicates.
			if (idx != NULL_IDX)
				return false;
			// idx is by reference, so this sets the value in grid.
			idx = this->list(arr_idx).size();
			this->list(arr_idx).push_back(pos);
			return true;
		}

		/**
		 * For given tracking list, set all lookup grid nodes to NULL index and
		 * clear the list.
		 *
		 * @param arr_idx tracking list id.
		 * @param lookup_idx lookup grid node tuple index.
		 */
		void reset(
			const UINT& arr_idx, const UINT& lookup_idx
		) {
			for (VecDi pos : this->list(arr_idx))
				this->get_internal(pos)[lookup_idx] = NULL_IDX;
			this->list(arr_idx).clear();
		}


		/**
		 * Remove pos at index idx in the array and set lookup at pos to NULL
		 * index.
		 *
		 * NOTE: idx passed by value since it changes indirectly via grid
		 * lookup.
		 * @param idx
		 * @param pos
		 */
		void remove(
			const UINT idx, const VecDi& pos, const UINT& arr_idx,
			const UINT& lookup_idx
		) {
			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
			this->assert_pos_bounds(pos, "remove: ");
			#endif

			// Set index lookup to null value.
			this->get_internal(pos)[lookup_idx] = NULL_IDX;

			// If this is not the last remaining position in the array, then
			// we must move the last position to this position and update the
			// lookup grid.
			const UINT& size = this->list(arr_idx).size();
			if (idx < size - 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = this->list(arr_idx)[size - 1];
				this->list(arr_idx)[idx] = pos_last;
				// Set the lookup grid to reference the new index in the array.
				this->get_internal(pos_last)[lookup_idx] = idx;
			}
			// Remove the last element in the array (which is at this point
			// either the last remaining element or a duplicate).
			this->list(arr_idx).pop_back();
		}
	};


	template <UINT D, UINT N, typename V, typename I>
	const typename LookupGridBase<D, N, V, I>::Idx_t
	LookupGridBase<D, N, V, I>::NULL_IDX_TUPLE =
		LookupGridBase<D, N, V, I>::Idx_t::Constant(
			std::numeric_limits<UINT>::max()
		);

	template <UINT D, UINT N, typename V, typename I>
	const UINT
	LookupGridBase<D, N, V, I>::NULL_IDX = std::numeric_limits<UINT>::max();


	/**
	 * Standard lookup grid.
	 *
	 * Holds a set of tracking lists storing grid positions, and a corresponding
	 * grid storing tuples of list indices.
	 *
	 * @tparam D the dimension of the grid.
	 * @tparam N the number of tracking lists to use.
	 */
	template <UINT D, UINT N=1>
	class LookupGrid
	: public LookupGridBase<
	  D, N, Eigen::Matrix<UINT, N, 1>, Eigen::Matrix<UINT, N, 1>
	>
	{
	protected:
		typedef LookupGridBase<
			D, N, Eigen::Matrix<UINT, N, 1>, Eigen::Matrix<UINT, N, 1>
		> Base_t;
		typedef typename Base_t::VecDu			VecDu;
		typedef typename Base_t::VecDi			VecDi;

	public:
		typedef typename Eigen::Matrix<UINT, N, 1>	Val_t;

		LookupGrid() : Base_t()
		{}


		LookupGrid(const VecDu& dims, const VecDi& offset)
		: Base_t(dims, offset)
		{}


		/**
		 * Return tuple of indices stored at given position in grid.
		 *
		 * Overloads get method of GridBase to simply return value stored.
		 *
		 * @param pos
		 * @return
		 */
		Val_t& get (const VecDi& pos)
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}
		const Val_t& get (const VecDi& pos) const
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}


	};


	/**
	 * Similar to LookupGrid but grid nodes store only a single list index.
	 *
	 * Useful in cases where grid nodes cannot be in more than one list.
	 *
	 * @tparam D the dimension of the grid.
	 * @tparam N the number of tracking lists to use.
	 */
	template <UINT D, UINT N=1>
	class SharedLookupGrid
	: public LookupGridBase<D, N, Eigen::Matrix<UINT, 1, 1>, UINT>
	{
	public:
		typedef SharedLookupGrid<D, N>						ThisType;
		typedef Eigen::Matrix<UINT, 1, 1>					Idx_t;
		typedef UINT										Val_t;
		typedef LookupGridBase<D, N, Idx_t, Val_t>			Base;
		typedef typename Base::VecDu						VecDu;
		typedef typename Base::VecDi						VecDi;
	public:

		SharedLookupGrid() : Base()
		{}


		SharedLookupGrid(const VecDu& dims, const VecDi& offset)
		: Base(dims, offset)
		{}

		/**
		 * Return index in associated list of grid node.
		 *
		 * Overloads get method of GridBase to return the first (and only)
		 * element in the tuple of indices.
		 *
		 * @param pos
		 * @return
		 */
		Val_t& get (const VecDi& pos)
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx)[0];
		}
		const Val_t& get (const VecDi& pos) const
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx)[0];
		}

		/**
		 * Return true if position currently tracked for given list id, false
		 * otherwise.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return
		 */
		const bool is_active (const VecDi& pos) const
		{
			return this->get(pos) != Base::NULL_IDX;
		}
		/**
		 * Add position to tracking list and store index in tracking list in
		 * grid.
		 *
		 * Overloads LookupGridBase to place lookup index in first (and only)
		 * element of tuple of indices at the grid position.
		 *
		 * @param pos position in grid.
		 * @param arr_idx tracking list id.
		 * @return
		 */
		bool add (const VecDi& pos, const UINT& arr_idx = 0)
		{
			return Base::add(pos, arr_idx, 0u);
		}

		/**
		 * For given tracking list, set all lookup grid nodes to NULL index and
		 * clear the list.
		 *
		 * Overloads LookupGridBase to reset first (and only) element of tuple
		 * of indices for each grid position referenced in given tracking list.
		 *
		 * @param arr_idx tracking list id.
		 */
		void reset (const UINT& arr_idx = 0)
		{
			Base::reset(arr_idx, 0u);
		}

		/**
		 * Remove an element from a tracking list by index and set it's
		 * corresponding grid node to NULL index.
		 *
		 * Overloads LookupGridBase to set the first (and only) element of
		 * tuple of indices in the grid to NULL index.
		 *
		 * @param idx index in tracking list.
		 * @param arr_idx tracking list id.
		 */
		void remove (const UINT& idx, const UINT& arr_idx = 0)
		{
			const VecDi& pos = this->list(arr_idx)[idx];
			Base::remove(idx, pos, arr_idx, 0);
		}

		/**
		 * Look up tracking list index in grid, remove from list and set grid
		 * node to NULL index.
		 *
		 * Overloads LookupGridBase to set the first (and only) element of
		 * tuple of indices in the grid to NULL index.
		 *
		 * @param pos position in lookup grid.
		 * @param arr_idx tracking list id.
		 */
		void remove (const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = this->get_internal(pos)[0];
			Base::remove(idx, pos, arr_idx, 0);
		}
	};


	/**
	 * Base class for a tracking grid, where grid nodes store arbitrary values
	 * and active nodes are tracked by a lookup grid.
	 *
	 * @tparam T the type of value stored in grid nodes.
	 * @tparam D the dimension of the grid.
	 * @tparam N the number of tracking lists to use.
	 * @tparam G the lookup grid type - either LookupGrid or SharedLookupGrid.
	 * @See TrackedGrid and SharedTrackedGrid.
	 */
	template <typename T, UINT D, UINT N, class G>
	class TrackedGridBase : public Grid<T, D>
	{
	public:
		typedef G								Lookup;
	protected:
		typedef	Grid<T, D>						Base;
		typedef typename Lookup::PosArray		PosArray;

		Lookup									m_grid_lookup;
	public:
		typedef typename Base::VecDu			VecDu;
		typedef typename Base::VecDi			VecDi;


		TrackedGridBase() : Base(), m_grid_lookup()
		{}


		TrackedGridBase(const VecDu& dims, const VecDi& offset)
		: Base(), m_grid_lookup()
		{
			this->init(dims, offset);
		}

		/**
		 * Get size of grid.
		 *
		 * @return
		 */
		const VecDu& dims () const
		{
			return Base::dims();
		}

		/**
		 * Reshape both grid of values and lookup grid.
		 *
		 * @param dims_new
		 */
		void dims (const VecDu& dims_new)
		{
			Base::dims(dims_new);
			m_grid_lookup.dims(dims_new);
		}

		/**
		 * Get grid's spatial offset.
		 *
		 * @return
		 */
		const VecDi offset () const
		{
			return Base::offset();
		}

		/**
		 * Set offset of grid of values and lookup grid.
		 *
		 * @param offset_new
		 */
		void offset (const VecDi& offset_new)
		{
			Base::offset(offset_new);
			m_grid_lookup.offset(offset_new);
		}

		/**
		 * Get lookup grid.
		 *
		 * @return
		 */
		Lookup& lookup()
		{
			return m_grid_lookup;
		}
		const Lookup& lookup() const
		{
			return m_grid_lookup;
		}

		/**
		 * Get list of active grid points from lookup grid.
		 *
		 * @param arr_idx tracking list id.
		 * @return
		 */
		PosArray& list(const UINT& arr_idx = 0)
		{
			return m_grid_lookup.list(arr_idx);
		}
		const PosArray& list(const UINT& arr_idx = 0) const
		{
			return m_grid_lookup.list(arr_idx);
		}

		/**
		 * Set value in grid at given position and add position to lookup grid.
		 *
		 * Will set value regardless whether lookup grid already set for given
		 * position + tracking list.
		 *
		 * @param pos position in grid.
		 * @param val value to set.
		 * @param arr_idx tracking list id.
		 * @return true if grid node set in lookup grid and position added to
		 * tracking list, false if grid node was already set so position already
		 * in a list.
		 */
		bool add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			this->get(pos) = val;
			return add(pos, arr_idx);
		}

		/**
		 * Add a position to the lookup grid.
		 *
		 * @param pos position in the grid to add.
		 * @param arr_idx tracking list id.
		 * @return true if grid node set in lookup grid and position added to
		 * tracking list, false if grid node was already set so position already
		 * in a list.
		 */
		bool add(const VecDi& pos, const UINT& arr_idx = 0)
		{
			return m_grid_lookup.add(pos, arr_idx);
		}

		/**
		 * Set every active grid node (i.e. those referenced by lookup grid)
		 * to given value and reset the lookup grid.
		 *
		 * Lookup grid will then be full of NULL indices and it's tracking
		 * list(s) will be empty.
		 *
		 * @param val value to set in main grid.
		 * @param arr_idx tracking list id to cycle over and clear.
		 */
		void reset(const T& val, const UINT& arr_idx = 0)
		{
			for (VecDi pos : m_grid_lookup.list(arr_idx))
				this->get(pos) = val;
			reset(arr_idx);
		}

		/**
		 * Reset a tracking list on the lookup grid.
		 *
		 * @param arr_idx tracking list to clear.
		 */
		void reset(const UINT& arr_idx = 0)
		{
			m_grid_lookup.reset(arr_idx);
		}

		/**
		 * Remove an element from a tracking list by index and set it's
		 * corresponding grid node in the tracking grid to NULL index.
		 *
		 * @param idx index in tracking list.
		 * @param arr_idx tracking list id.
		 */
		void remove(const UINT& idx, const UINT& arr_idx = 0)
		{
			m_grid_lookup.remove(idx, arr_idx);
		}

		/**
		 * Look up tracking list index in grid, remove from list and set lookup
		 * grid node to NULL index.
		 *
		 * @param pos position in lookup grid.
		 * @param arr_idx tracking list id.
		 */
		void remove (const VecDi& pos, const UINT& arr_idx = 0)
		{
			m_grid_lookup.remove(pos, arr_idx);
		}

		/**
		 * Return true if position currently tracked for given list id, false
		 * otherwise.
		 *
		 * @param pos
		 * @param arr_idx
		 * @return
		 */
		const bool is_active (const VecDi& pos, const UINT& arr_idx = 0) const
		{
			return m_grid_lookup.is_active(pos, arr_idx);
		}
	};


	/**
	 * Standard TrackedGridBase with multiple lookup indices per grid node, one
	 * for each tracking list.
	 */
	template <typename T, UINT D, UINT N=1>
	using TrackedGrid = TrackedGridBase<T, D, N, LookupGrid<D, N> >;


	/**
	 * TrackedGridBase where each node of the lookup grid stores only a single
	 * list index.
	 *
	 * Useful where a grid node can only be in one of the tracking lists.
	 */
	template <typename T, UINT D, UINT N=1>
	using SharedTrackedGrid = TrackedGridBase<T, D, N, SharedLookupGrid<D, N> >;
}
#endif /* MAPPEDGRID_HPP_ */
