#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <stdexcept>
#include "MappedGrid.hpp"

namespace felt
{
	template <typename T, UINT D, UINT P, typename Part_t=Grid<T, D> >
	class PartitionedGrid : public Grid<T, D>
	{
	public:
		typedef Part_t	Partition_t;
	protected:
		typedef Grid<T, D>					Grid_t;
		typedef typename Grid_t::VecDu		VecDu;
		typedef typename Grid_t::VecDi		VecDi;

		ArrayMappedGrid<Part_t, D>	m_grid_parts;

	public:
		typedef ArrayMappedGrid<Part_t, D>	PartGrid_t;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

		PartitionedGrid () : Grid_t::Grid(), m_grid_parts()
		{}

		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Grid_t::Grid(), m_grid_parts()
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
			return Grid_t::index(pos, m_grid_parts.dims(), m_grid_parts.offset());
		}


		const VecDi index_parent (const UINT& idx) const
		{
			return Grid_t::index(idx, m_grid_parts.dims(), m_grid_parts.offset());
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
			return Grid_t::dims();
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
				Partition_t& part = m_grid_parts.data()(idx);
				part.dims(udims_child);
			}
		}


		void offset (const VecDi& offset_grid)
		{
			Grid_t::offset(offset_grid);

			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_parts.offset(offset_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Partition_t& part = m_grid_parts.data()(idx);
				const VecDi& pos_part = this->index_parent(idx);
				const VecDi& offset_part = (
					(
						(pos_part.array() - m_grid_parts.offset().array())
						* idims_child.array()
					).matrix() + offset_grid
				);

				part.offset(offset_part);
			}
		}


		const VecDi offset () const
		{
			return Grid_t::offset();
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
				Partition_t& part = m_grid_parts.data()(idx);
				part.fill(val);
			}
		}


		T& get (const VecDi& pos)
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			Partition_t& part = m_grid_parts.data()(idx_parent);
			return part.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			const Partition_t& part = m_grid_parts.data()(idx_parent);
			return part.get(pos);
		}

	};


	template <typename T, UINT D, UINT P, typename Part_t>
		const typename PartitionedGrid<T, D, P, Part_t>::VecDu
			PartitionedGrid<T, D, P, Part_t>::udims_child(P, P, P);
	template <typename T, UINT D, UINT P, typename Part_t>
		const typename PartitionedGrid<T, D, P, Part_t>::VecDi
			PartitionedGrid<T, D, P, Part_t>::idims_child(P, P, P);

}
#endif /* PARTITIONEDGRID_HPP_ */
