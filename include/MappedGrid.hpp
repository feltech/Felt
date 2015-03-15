#ifndef MAPPEDGRID_HPP_
#define MAPPEDGRID_HPP_

#include "Grid.hpp"

namespace felt
{
	template <
		typename T, UINT D, UINT N=1
	>
	class ArrayMappedGrid : public felt::Grid<T, D>
	{
	protected:
		typedef Grid<T, D>				Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		typedef std::vector<VecDi, Eigen::aligned_allocator<T>>	PosArray;
	protected:
		PosArray	m_aPos[N];

	public:

		/**
		 * Initialise a grid with given dimension and offset.
		 *
		 * @param dims
		 * @param offset
		 * @param delta
		 */
		ArrayMappedGrid (
			const VecDu& dims, const VecDi& offset
		) : Grid_t::Grid(dims, offset)
		{
		}


		ArrayMappedGrid () : Grid_t::Grid()
		{
		}


		const PosArray& list(const UINT& arr_idx = 0) const
		{
			return m_aPos[arr_idx];
		}


		PosArray& list(const UINT& arr_idx = 0)
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
			const UINT& size = m_aPos[arr_idx].size();
			if (size > 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = m_aPos[arr_idx][size - 1];
				m_aPos[arr_idx][idx] = pos_last;
			}
			list(arr_idx).pop_back();
		}
	};


	template <UINT D, UINT N=1>
	class PosArrayMappedGrid : public ArrayMappedGrid<UINT, D, N>
	{
	protected:
		typedef ArrayMappedGrid<UINT, D, N> Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		static const UINT NULL_IDX;

		static inline const UINT num_lists()
		{
			return N;
		}

		PosArrayMappedGrid(const VecDu& dims, const VecDi& offset)
		: Grid_t(dims, offset)
		{
			this->fill(NULL_IDX);
		}

		UINT add(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = (*this)(pos);
			// Do not allow duplicates.
			if (idx != NULL_IDX)
				return idx;
			return this->add(pos, this->list(arr_idx).size(), arr_idx);
		}

		void reset(const UINT& arr_idx = 0)
		{
			reset(NULL_IDX, arr_idx);
		}

		void remove(const UINT& idx, const UINT& arr_idx = 0)
		{
			const VecDi& pos = this->list(arr_idx)[idx];
			remove(idx, pos);
		}

		void remove(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const UINT& idx = (*this)(pos);
			remove(idx, pos, arr_idx);
		}

	protected:
		UINT add(const VecDi& pos, const UINT& val, const UINT& arr_idx)
		{
			return Grid_t::add(pos, val, arr_idx);
		}

		void reset(const UINT& val, const UINT& arr_idx)
		{
			Grid_t::reset(val, arr_idx);
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
		void remove(const UINT idx, const VecDi& pos, const UINT& arr_idx = 0)
		{
			// Set index lookup to null value.
			(*this)(pos) = NULL_IDX;

			// If this is not the last remaining position in the array, then
			// we must move the last position to this position and update the
			// lookup grid.
			const UINT& size = this->m_aPos[arr_idx].size();
			if (size > 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = this->list(arr_idx)[size - 1];
				this->list(arr_idx)[idx] = pos_last;
				// Set the lookup grid to reference the new index in the array.
				(*this)(pos_last) = idx;
			}
			// Remove the last element in the array (which is at this point
			// either the last remaining element or a duplicate).
			this->list(arr_idx).pop_back();
		}
	};

	template <UINT D, UINT N> const UINT
	PosArrayMappedGrid<D, N>::NULL_IDX = std::numeric_limits<UINT>::max();
}
#endif /* MAPPEDGRID_HPP_ */
