#ifndef MAPPEDGRID_HPP_
#define MAPPEDGRID_HPP_

#include "Grid.hpp"

namespace felt
{
	template <
		typename T, UINT D, UINT N=1, typename R=T
	>
	class ArrayMappedGrid : public felt::GridBase<T, D, R>
	{
	protected:
		typedef GridBase<T, D, R>	Grid_t;
		typedef typename Grid_t::VecDu		VecDu;
		typedef typename Grid_t::VecDi		VecDi;
	public:
		typedef std::vector<VecDi, Eigen::aligned_allocator<T> >	PosArray;
	protected:
		std::array<PosArray, N>	m_aPos;

	public:

		virtual ~ArrayMappedGrid () {}
		/**
		 * Initialise a grid with given dimension and offset.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		ArrayMappedGrid (
			const VecDu& dims, const VecDi& offset
		) : Grid_t(dims, offset)
		{
		}


		ArrayMappedGrid () : Grid_t()
		{
		}


		virtual R& get (const VecDi& pos)
		{
			return this->get_internal(pos);
		}


		virtual const R& get (const VecDi& pos) const
		{
			return this->get_internal(pos);
		}


		virtual inline const PosArray& list(const UINT& arr_idx = 0) const
		{
			return m_aPos[arr_idx];
		}


		virtual inline PosArray& list(const UINT& arr_idx = 0)
		{
			return m_aPos[arr_idx];
		}


		UINT add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			(*this)(pos) = val;
			list(arr_idx).push_back(pos);
			return list(arr_idx).size() - 1;
		}

		void reset(const T& val, const UINT& arr_idx = 0)
		{
			for (VecDi pos : list(arr_idx))
				(*this)(pos) = val;
			list(arr_idx).clear();
		}


		void remove(const UINT& idx, const UINT& arr_idx = 0)
		{
			const UINT& size = list(arr_idx).size();
			if (size > 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = list(arr_idx)[size - 1];
				list(arr_idx)[idx] = pos_last;
			}
			list(arr_idx).pop_back();
		}
	};


	template <UINT D, UINT N, typename I, typename V>
	class PosArrayMappedGridBase
	: public GridBase<I, D, V>
	{
	public:
		typedef I								Idx_t;
		typedef V								Val_t;
	protected:
		typedef GridBase<Idx_t, D, Val_t>		Base_t;
		typedef typename Base_t::VecDu			VecDu;
		typedef typename Base_t::VecDi			VecDi;
	public:
		typedef std::vector<
			VecDi, Eigen::aligned_allocator<VecDi>
		> PosArray;
	protected:
		std::array<PosArray, N>	m_aPos;
	public:

		static const Idx_t				 			NULL_IDX_TUPLE;
		static const UINT				 			NULL_IDX;

		static inline const UINT num_lists()
		{
			return N;
		}

		virtual ~PosArrayMappedGridBase() {}

		PosArrayMappedGridBase() : Base_t()
		{}


		PosArrayMappedGridBase(const VecDu& dims, const VecDi& offset)
		: Base_t()
		{
			this->init(dims, offset);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			Base_t::init(dims, offset, delta);
			this->fill(NULL_IDX_TUPLE);
		}


		virtual inline const PosArray& list(const UINT& arr_idx = 0) const
		{
			return m_aPos[arr_idx];
		}


		virtual inline PosArray& list(const UINT& arr_idx = 0)
		{
			return m_aPos[arr_idx];
		}


		virtual UINT add(const VecDi& pos, const UINT& arr_idx = 0) = 0;


		virtual void reset(const UINT& arr_idx = 0) = 0;


		virtual void remove(const UINT& idx, const UINT& arr_idx = 0) = 0;


		virtual void remove(const VecDi& pos, const UINT& arr_idx = 0) = 0;

	protected:

		UINT add(
			const VecDi& pos, const UINT& arr_idx, const UINT& lookup_idx
		) {
			UINT& idx = this->get_internal(pos)[lookup_idx];
			// Do not allow duplicates.
			if (idx != NULL_IDX)
				return idx;
			idx = this->list(arr_idx).size();
			this->list(arr_idx).push_back(pos);
			return idx;
		}


		void reset(
			const UINT& arr_idx, const UINT& lookup_idx
		) {
			for (VecDi pos : this->list(arr_idx))
				this->get_internal(pos)[lookup_idx] = NULL_IDX;
			this->list(arr_idx).clear();
		}


		/**
		 * Remove pos at index idx in the array and set lookup at pos to null
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
			// Set index lookup to null value.
			this->get_internal(pos)[lookup_idx] = NULL_IDX;

			// If this is not the last remaining position in the array, then
			// we must move the last position to this position and update the
			// lookup grid.
			const UINT& size = this->list(arr_idx).size();
			if (size > 1)
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
	const typename PosArrayMappedGridBase<D, N, V, I>::Idx_t
	PosArrayMappedGridBase<D, N, V, I>::NULL_IDX_TUPLE =
		PosArrayMappedGridBase<D, N, V, I>::Idx_t::Constant(
			std::numeric_limits<UINT>::max()
		);

	template <UINT D, UINT N, typename V, typename I>
	const UINT
	PosArrayMappedGridBase<D, N, V, I>::NULL_IDX = std::numeric_limits<UINT>::max();


	template <UINT D, UINT N=1>
	class PosArrayMappedGrid
	: public PosArrayMappedGridBase<
	  D, N, Eigen::Matrix<UINT, N, 1>, Eigen::Matrix<UINT, N, 1>
	>
	{
	protected:
		typedef PosArrayMappedGridBase<
			D, N, Eigen::Matrix<UINT, N, 1>, Eigen::Matrix<UINT, N, 1>
		> Base_t;
		typedef typename Base_t::VecDu			VecDu;
		typedef typename Base_t::VecDi			VecDi;

	public:
		typedef typename Eigen::Matrix<UINT, N, 1>	Val_t;

		virtual ~PosArrayMappedGrid() {}

		PosArrayMappedGrid() : Base_t()
		{}


		PosArrayMappedGrid(const VecDu& dims, const VecDi& offset)
		: Base_t(dims, offset)
		{}

		const Val_t& get (const VecDi& pos) const
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}

		Val_t& get (const VecDi& pos)
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx);
		}


		UINT add (const VecDi& pos, const UINT& arr_idx = 0)
		{
			return Base_t::add(pos, arr_idx, arr_idx);
		}


		void reset (const UINT& arr_idx = 0)
		{
			Base_t::reset(arr_idx, arr_idx);
		}


		void remove (const UINT& idx, const UINT& arr_idx = 0)
		{
			const VecDi& pos = this->list(arr_idx)[idx];
			Base_t::remove(idx, pos, arr_idx, arr_idx);
		}


		void remove (const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = this->get_internal(pos)[arr_idx];
			Base_t::remove(idx, pos, arr_idx, arr_idx);
		}
	};


	template <UINT D, UINT N=1>
	class PosArrayMappedSharedGrid
	: public PosArrayMappedGridBase<D, N, Eigen::Matrix<UINT, 1, 1>, UINT>
	{
	public:
		typedef Eigen::Matrix<UINT, 1, 1>					Idx_t;
		typedef UINT										Val_t;
	protected:
		typedef PosArrayMappedGridBase<D, N, Idx_t, Val_t>	Base_t;
		typedef typename Base_t::VecDu						VecDu;
		typedef typename Base_t::VecDi						VecDi;

	public:
		virtual ~PosArrayMappedSharedGrid() {}

		PosArrayMappedSharedGrid() : Base_t()
		{}


		PosArrayMappedSharedGrid(const VecDu& dims, const VecDi& offset)
		: Base_t(dims, offset)
		{}

		const Val_t& get (const VecDi& pos) const
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx)[0];
		}

		Val_t& get (const VecDi& pos)
		{
			const UINT& idx = this->index(pos);
			return this->data()(idx)[0];
		}

		UINT add (const VecDi& pos, const UINT& arr_idx = 0)
		{
			return Base_t::add(pos, arr_idx, 0u);
		}


		void reset (const UINT& arr_idx = 0)
		{
			Base_t::reset(arr_idx, 0u);
		}


		void remove (const UINT& idx, const UINT& arr_idx = 0)
		{
			const VecDi& pos = this->list(arr_idx)[idx];
			Base_t::remove(idx, pos, arr_idx, 0);
		}


		void remove (const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = this->get_internal(pos)[0];
			Base_t::remove(idx, pos, arr_idx, 0);
		}
	};
}
#endif /* MAPPEDGRID_HPP_ */
