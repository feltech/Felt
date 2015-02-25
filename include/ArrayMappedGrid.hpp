/*
 * ArrayMappedGrid.hpp
 *
 *  Created on: 23 Feb 2015
 *      Author: dave
 */

#ifndef ARRAYMAPPEDGRID_HPP_
#define ARRAYMAPPEDGRID_HPP_

#include "Grid.hpp"

namespace felt
{
	template <typename T, UINT D>
	class ArrayMappedGrid : public felt::Grid<T, D>
	{
	protected:
		typedef Grid<T, D>				Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		typedef std::vector<VecDi>		PosArray;
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


		const PosArray& list()
		{
			return m_aPos;
		}


		void add(const VecDi& pos, T val)
		{
			(*this)(pos) = val;
			this->add(pos);
		}


		void add(const VecDi& pos)
		{
			m_aPos.push_back(pos);
		}


		void reset(const T& val)
		{
			for (VecDi pos : m_aPos)
				(*this)(pos) = val;
			m_aPos.clear();
		}
	};


	template <typename T, UINT D>
	class GridMappedArray : public std::vector<T>
	{
	protected:
		typedef Grid<UINT, D>			IdxGrid;
		typedef std::vector<T>			Array_t;
		typedef typename IdxGrid::VecDu	VecDu;
		typedef typename IdxGrid::VecDi	VecDi;

		IdxGrid 	m_grid_idx;
	public:
		static const UINT NULL_IDX;

		GridMappedArray(const VecDu& dims, const VecDi& offset)
		: m_grid_idx(dims, offset)
		{
			this->reserve(dims.norm());
			m_grid_idx.fill(NULL_IDX);
		}

		void add(const VecDi& pos, const T& val)
		{
			this->push_back(val);
			m_grid_idx(pos) = this->size() - 1;
		}

		const UINT idx(const VecDi pos)
		{
			return m_grid_idx(pos);
		}
	};

	template <typename T, UINT D>
	const UINT
	GridMappedArray<T, D>::NULL_IDX = std::numeric_limits<UINT>::max();

}




#endif /* ARRAYMAPPEDGRID_HPP_ */
