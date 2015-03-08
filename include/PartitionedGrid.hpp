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

		VecDu	m_dims_parent;
		VecDi	m_offset_parent;


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


		static const VecDi pos_to_partn (const VecDi& pos)
		{
			return (pos.array() / idims_child.array()).matrix();
		}
		
		
		static const VecDu dims_to_partn (const VecDu& pos)
		{
			return (pos.array() / udims_child.array()).matrix();
		}
		

		const UINT index_parent (const VecDi& pos) const
		{
			return Grid_t::index(pos, m_dims_parent, m_offset_parent);
		}


		const VecDi index_parent (const UINT& idx) const
		{
			return Grid_t::index(idx, m_dims_parent, m_offset_parent);
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
		const VecDu dims (const VecDu& dims_new)
		{
			m_dims_parent = dims_to_partn(dims_new);

			VecDu dims_old = this->m_vec_dims;
			this->m_vec_dims = dims_new;
			VecDu dims_parent = (
				dims_new.array() / udims_child.array()
			).matrix();

			m_grid_parts.dims(dims_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Partition_t& part = m_grid_parts.data()(idx);
				part.dims(udims_child);
			}

			// Return old dimensions.
			return dims_old;
		}


		void offset (const VecDi& vec_offset)
		{
			m_offset_parent = pos_to_partn(vec_offset);
			m_grid_parts.offset(m_offset_parent);
			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				Partition_t& part = m_grid_parts.data()(idx);
				const VecDi& pos_part = this->index_parent(idx);
				const VecDi& offset_part = (
					(pos_part.array() * idims_child.array()).matrix()
					+ (vec_offset.array() / idims_child.array()).matrix()
				);

				part.offset(offset_part);
			}
			Grid_t::offset(vec_offset);
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
