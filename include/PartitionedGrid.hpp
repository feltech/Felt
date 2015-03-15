#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <stdexcept>
#include "MappedGrid.hpp"

namespace felt
{
	template <
		typename T, UINT D, UINT P, class G=Grid<T, D>, class PG=Grid<G, D>
	>
	class PartitionedGrid : public G
	{
	protected:
		typedef G			Base_t;
	public:
		typedef PG			PartGrid_t;
	protected:
		typedef typename Base_t::VecDu			VecDu;
		typedef typename Base_t::VecDi			VecDi;

		PartGrid_t	m_grid_parts;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

	public:

		PartitionedGrid () : Base_t(), m_grid_parts()
		{}

		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base_t(), m_grid_parts()
		{

			if ((dims / P) * P != dims)
			{
				throw std::invalid_argument(
					"dims must be divisible by partition size."
				);
			}

			this->init(dims, offset, delta);
		}


		const VecDi pos_to_partn (const VecDi& pos) const
		{
			return (
				(pos - this->offset()).array() / idims_child.array()
			).matrix() + m_grid_parts.offset();
		}


		const UINT index_parent (const VecDi& pos) const
		{
			return Base_t::index(
				pos, m_grid_parts.dims(), m_grid_parts.offset()
			);
		}


		const VecDi index_parent (const UINT& idx) const
		{
			return Base_t::index(
				idx, m_grid_parts.dims(), m_grid_parts.offset()
			);
		}
		

		PartGrid_t& parts ()
		{
			return m_grid_parts;
		}


		const PartGrid_t& parts () const
		{
			return m_grid_parts;
		}

		/**
		 * Get grid dimensions.
		 *
		 * @return
		 */
		const VecDu& dims () const
		{
			return Base_t::dims();
		}

		/**
		 * Reshape grid.
		 *
		 * @param vec_NewDims
		 * @return
		 */
		void dims (const VecDu& dims_grid)
		{
			this->m_vec_dims = dims_grid;

			VecDu dims_parent = (
				dims_grid.array() / udims_child.array()
			).matrix();

			m_grid_parts.dims(dims_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Base_t& part = m_grid_parts.data()(idx);
				part.dims(udims_child);
			}
		}


		void offset (const VecDi& offset_grid)
		{
			Base_t::offset(offset_grid);

			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_parts.offset(offset_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Base_t& part = m_grid_parts.data()(idx);
				const VecDi& pos_part = this->index_parent(idx);
				const VecDi& offset_part = (
					(
						(pos_part - m_grid_parts.offset()).array()
						* idims_child.array()
					).matrix() + offset_grid
				);

				part.offset(offset_part);
			}
		}


		const VecDi offset () const
		{
			return Base_t::offset();
		}


		/**
		 * Fill with a single value.
		 *
		 * @param val
		 */
		void fill (const T& val)
		{
			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Base_t& part = m_grid_parts.data()(idx);
				part.fill(val);
			}
		}


		T& get (const VecDi& pos)
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			Base_t& part = m_grid_parts.data()(idx_parent);
			return part.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			const Base_t& part = m_grid_parts.data()(idx_parent);
			return part.get(pos);
		}

	};


	template <typename T, UINT D, UINT P, class G, class PG>
		const typename PartitionedGrid<T, D, P, G, PG>::VecDu
			PartitionedGrid<T, D, P, G, PG>::udims_child(P, P, P);
	template <typename T, UINT D, UINT P, class G, class PG>
		const typename PartitionedGrid<T, D, P, G, PG>::VecDi
			PartitionedGrid<T, D, P, G, PG>::idims_child(P, P, P);



	/**
	 * Specialisation for ArrayMappedGrid.
	 */
	template <typename T, UINT D, UINT P, UINT N>
	class MappedPartitionedGrid
	: public PartitionedGrid <
		T, D, P, ArrayMappedGrid<T, D, N>,
		ArrayMappedGrid<ArrayMappedGrid<T, D, N>, D, N>
	>
	{
	protected:
		typedef PartitionedGrid <
			T, D, P, ArrayMappedGrid<T, D, N>,
			ArrayMappedGrid<ArrayMappedGrid<T, D, N>, D, N>
		> Base_t;
		typedef typename Base_t::Grid_t	Grid_t;
		typedef typename Grid_t::VecDu	VecDu;
		typedef typename Grid_t::VecDi	VecDi;
	public:
		typedef typename Base_t::PartGrid_t	PartGrid_t;

		MappedPartitionedGrid () : Base_t()
		{}

		MappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base_t(dims, offset, delta)
		{}

		UINT add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			const VecDi& part_pos = this->pos_to_partn(pos);
			const UINT& idx_parent = this->index_parent(part_pos);

			this->list(arr_idx).push_back(part_pos);

			return this->m_grid_parts.data()(idx_parent).add(pos, val, arr_idx);
		}
	};

}
#endif /* PARTITIONEDGRID_HPP_ */
