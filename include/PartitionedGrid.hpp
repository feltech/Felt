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
		typedef G	ChildGrid_t;
	public:
		typedef PG	ParentGrid_t;
	protected:
		typedef typename ChildGrid_t::VecDu	VecDu;
		typedef typename ChildGrid_t::VecDi	VecDi;

		ParentGrid_t	m_grid_parts;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

	public:

		PartitionedGrid () : ChildGrid_t(), m_grid_parts()
		{}

		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : ChildGrid_t(), m_grid_parts()
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
			return ChildGrid_t::index(
				pos, m_grid_parts.dims(), m_grid_parts.offset()
			);
		}


		const VecDi index_parent (const UINT& idx) const
		{
			return ChildGrid_t::index(
				idx, m_grid_parts.dims(), m_grid_parts.offset()
			);
		}
		

		ParentGrid_t& parts ()
		{
			return m_grid_parts;
		}


		const ParentGrid_t& parts () const
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
			return ChildGrid_t::dims();
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

			const VecDu& dims_parent = (
				dims_grid.array() / udims_child.array()
			).matrix();

			m_grid_parts.dims(dims_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				ChildGrid_t& part = m_grid_parts.data()(idx);
				part.dims(udims_child);
			}
		}


		void offset (const VecDi& offset_grid)
		{
			ChildGrid_t::offset(offset_grid);

			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_parts.offset(offset_parent);

			for (UINT idx = 0; idx < m_grid_parts.data().size(); idx++)
			{
				ChildGrid_t& part = m_grid_parts.data()(idx);
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
			return ChildGrid_t::offset();
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
				ChildGrid_t& part = m_grid_parts.data()(idx);
				part.fill(val);
			}
		}


		T& get (const VecDi& pos)
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			ChildGrid_t& part = m_grid_parts.data()(idx_parent);
			return part.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			const ChildGrid_t& part = m_grid_parts.data()(idx_parent);
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
	template <typename T, UINT D, UINT N>
	using ChildGridClass = ArrayMappedGrid<T, D, N>;

	template <typename T, UINT D, UINT N>
	using ParentGridClass = Grid<ChildGridClass<T, D, N>, D>;

	template <typename T, UINT D, UINT P, UINT N>
	using BaseClass = PartitionedGrid<
		T, D, P, ChildGridClass<T, D, N>, ParentGridClass<T, D, N>
	>;

	template <typename T, UINT D, UINT P, UINT N>
	class MappedPartitionedGrid : public BaseClass<T, D, P, N>
	{
	protected:
		typedef BaseClass<T, D, P, N>			Base_t;
		typedef typename Base_t::ChildGrid_t	ChildGrid_t;
		typedef typename ChildGrid_t::VecDu		VecDu;
		typedef typename ChildGrid_t::VecDi		VecDi;
	public:
		typedef typename Base_t::ParentGrid_t	ParentGrid_t;


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


		void reset(const T& val, const UINT& arr_idx = 0)
		{
			for (const VecDi& parent_pos : this->list(arr_idx))
			{
				ParentGrid_t& parent = this->parts();
				ChildGrid_t& child = parent(parent_pos);

				for (const VecDi& child_pos : child.list(arr_idx))
					child(child_pos) = val;

				child.list(arr_idx).clear();
			}

			this->list(arr_idx).clear();
		}
	};

}
#endif /* PARTITIONEDGRID_HPP_ */
