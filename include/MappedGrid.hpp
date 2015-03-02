#ifndef MAPPEDGRID_HPP_
#define MAPPEDGRID_HPP_

#include "Grid.hpp"

namespace felt
{
	template <typename T, UINT D, typename A = std::allocator<T>>
	class ArrayMappedGrid : public felt::Grid<T, D>
	{
	protected:
		typedef Grid<T, D>				Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		typedef std::vector<VecDi, A>		PosArray;
	protected:
		PosArray	m_aPos;

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


		const PosArray& list() const
		{
			return m_aPos;
		}


		UINT add(const VecDi& pos, const T& val)
		{
			(*this)(pos) = val;
			return add(pos);
		}


		UINT add(const VecDi& pos)
		{
			m_aPos.push_back(pos);
			return m_aPos.size() - 1;
		}


		void reset(const T& val)
		{
			for (VecDi pos : m_aPos)
				(*this)(pos) = val;
			reset();
		}


		void reset()
		{
			m_aPos.clear();
		}


		void remove(const UINT& idx)
		{
			const UINT& size = m_aPos.size();
			if (size > 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = m_aPos[size - 1];
				m_aPos[idx] = pos_last;
			}
			m_aPos.pop_back();
		}
	};


	template <UINT D>
	class PosArrayMappedGrid
	: public ArrayMappedGrid<UINT, D>
	{
	protected:
		typedef ArrayMappedGrid<UINT, D> Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		static const UINT NULL_IDX;

		PosArrayMappedGrid(const VecDu& dims, const VecDi& offset)
		: Grid_t(dims, offset)
		{
			this->fill(NULL_IDX);
		}

		UINT add(const VecDi& pos)
		{
			const UINT idx = (*this)(pos);
			// Do not allow duplicates.
			if (idx != NULL_IDX)
				return idx;
			return this->add(pos, this->m_aPos.size());
		}

		void reset()
		{
			reset(NULL_IDX);
		}

		void remove(const UINT& idx)
		{
			const VecDi& pos = this->m_aPos[idx];
			remove(idx, pos);
		}

		void remove(const VecDi& pos)
		{
			const UINT& idx = (*this)(pos);
			remove(idx, pos);
		}

	private:
		UINT add(const VecDi& pos, const UINT& val)
		{
			return Grid_t::add(pos, val);
		}

		void reset(const UINT& val)
		{
			Grid_t::reset(NULL_IDX);
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
		void remove(const UINT idx, const VecDi& pos)
		{
			// Set index lookup to null value.
			(*this)(pos) = NULL_IDX;

			// If this is not the last remaining position in the array, then
			// we must move the last position to this position and update the
			// lookup grid.
			const UINT& size = this->m_aPos.size();
			if (size > 1)
			{
				// Duplicate last element into this index.
				const VecDi& pos_last = this->m_aPos[size - 1];
				this->m_aPos[idx] = pos_last;
				// Set the lookup grid to reference the new index in the array.
				(*this)(pos_last) = idx;
			}
			// Remove the last element in the array (which is at this point
			// either the last remaining element or a duplicate).
			this->m_aPos.pop_back();
		}
	};

	template <UINT D> const UINT
	PosArrayMappedGrid<D>::NULL_IDX = std::numeric_limits<UINT>::max();



	template <typename T, UINT D>
	class SpatiallyPartitionedArray : public felt::Grid<std::vector<T>, D>
	{
	protected:
		typedef PosArrayMappedGrid<D> ActivePartitionGrid;
		typedef Grid<std::vector<T>, D> Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;

		ActivePartitionGrid m_grid_active;
	public:
		static const UINT NULL_IDX;


		SpatiallyPartitionedArray(const VecDu& dims, const VecDi& offset)
		: Grid_t(dims, offset), m_grid_active(dims, offset)
		{

		}


		UINT add(const VecDi& pos, const T& val)
		{
			(*this)(pos).push_back(val);
			return m_grid_active.add(pos);
		}

		ActivePartitionGrid& active()
		{
			return m_grid_active;
		}


		const ActivePartitionGrid& active() const
		{
			return m_grid_active;
		}


		void reset()
		{
			for (VecDi pos : m_grid_active.list())
				(*this)(pos).clear();
			m_grid_active.reset();
		}
	};

	template <typename T, UINT D> const UINT
	SpatiallyPartitionedArray<T, D>::NULL_IDX = PosArrayMappedGrid<D>::NULL_IDX;

}
#endif /* MAPPEDGRID_HPP_ */
