#ifndef PARTITIONEDGRID_HPP_
#define PARTITIONEDGRID_HPP_

#include <stdexcept>
#include "MappedGrid.hpp"

namespace felt
{
	template <
		typename T, UINT D, UINT P, class G=Grid<T, D>
	>
	class PartitionedGrid : public G
	{
	protected:
		typedef G	ChildGrid_t;
	public:
		typedef Grid<ChildGrid_t, D>	ParentGrid_t;
	protected:
		typedef typename ChildGrid_t::VecDu	VecDu;
		typedef typename ChildGrid_t::VecDi	VecDi;

		ParentGrid_t	m_grid_parent;

		static const VecDu	udims_child;
		static const VecDi	idims_child;

	public:
		virtual ~PartitionedGrid ()
		{}

		PartitionedGrid () : ChildGrid_t(), m_grid_parent()
		{}


		PartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : ChildGrid_t(), m_grid_parent()
		{
			this->init(dims, offset, delta);
		}


		virtual void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			if ((dims / P) * P != dims)
			{
				throw std::invalid_argument(
					"dims must be divisible by partition size."
				);
			}

			ChildGrid_t::init(dims, offset, delta);
		}


		const VecDi pos_to_partn (const VecDi& pos) const
		{
			return (
				(pos - this->offset()).array() / idims_child.array()
			).matrix() + m_grid_parent.offset();
		}


		const UINT index_parent (const VecDi& pos) const
		{
			return ChildGrid_t::index(
				pos, m_grid_parent.dims(), m_grid_parent.offset()
			);
		}


		const VecDi index_parent (const UINT& idx) const
		{
			return ChildGrid_t::index(
				idx, m_grid_parent.dims(), m_grid_parent.offset()
			);
		}
		

		ParentGrid_t& parent ()
		{
			return m_grid_parent;
		}


		const ParentGrid_t& parent () const
		{
			return m_grid_parent;
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
		virtual void dims (const VecDu& dims_grid)
		{
			this->m_vec_dims = dims_grid;

			const VecDu& dims_parent = (
				dims_grid.array() / udims_child.array()
			).matrix();

			m_grid_parent.dims(dims_parent);

			for (UINT idx = 0; idx < m_grid_parent.data().size(); idx++)
			{
				ChildGrid_t& part = m_grid_parent.data()(idx);
				part.dims(udims_child);
			}
		}


		virtual void offset (const VecDi& offset_grid)
		{
			ChildGrid_t::offset(offset_grid);

			const VecDi& offset_parent = (
				offset_grid.array() / idims_child.array()
			).matrix();

			m_grid_parent.offset(offset_parent);

			for (UINT idx = 0; idx < m_grid_parent.data().size(); idx++)
			{
				ChildGrid_t& part = m_grid_parent.data()(idx);
				const VecDi& pos_part = this->index_parent(idx);
				const VecDi& offset_part = (
					(
						(pos_part - m_grid_parent.offset()).array()
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
			for (UINT idx = 0; idx < m_grid_parent.data().size(); idx++)
			{
				ChildGrid_t& part = m_grid_parent.data()(idx);
				part.fill(val);
			}
		}


		T& get (const VecDi& pos)
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			ChildGrid_t& part = m_grid_parent.data()(idx_parent);
			return part.get(pos);
		}


		const T& get (const VecDi& pos) const
		{
			const UINT& idx_parent = this->index_parent(pos_to_partn(pos));
			const ChildGrid_t& part = m_grid_parent.data()(idx_parent);
			return part.get(pos);
		}

	};


	template <typename T, UINT D, UINT P, class G>
		const typename PartitionedGrid<T, D, P, G>::VecDu
			PartitionedGrid<T, D, P, G>::udims_child(P, P, P);
	template <typename T, UINT D, UINT P, class G>
		const typename PartitionedGrid<T, D, P, G>::VecDi
			PartitionedGrid<T, D, P, G>::idims_child(P, P, P);



	/**
	 * Specialisation for ArrayMappedGrid.
	 */

	template <typename T, UINT D, UINT P, UINT N>
	using MappedPartitionedGridBase = PartitionedGrid<
		T, D, P, ArrayMappedGrid<T, D, N>
	>;

	template <typename T, UINT D, UINT P, UINT N>
	class MappedPartitionedGrid : public MappedPartitionedGridBase<T, D, P, N>
	{
	protected:
		typedef MappedPartitionedGridBase<T, D, P, N>	Base_t;
		typedef typename Base_t::ChildGrid_t			ChildGrid_t;
		typedef typename ChildGrid_t::VecDu				VecDu;
		typedef typename ChildGrid_t::VecDi				VecDi;
	public:
		typedef typename Base_t::ParentGrid_t			ParentGrid_t;


		MappedPartitionedGrid () : Base_t()
		{}


		MappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base_t(dims, offset, delta)
		{}


		UINT add(const VecDi& pos, const T& val, const UINT& arr_idx = 0)
		{
			const VecDi& child_pos = this->pos_to_partn(pos);
			const UINT& idx_child = this->index_parent(child_pos);
			ParentGrid_t& parent = this->parent();

			this->list(arr_idx).push_back(child_pos);

			return parent.data()(idx_child).add(pos, val, arr_idx);
		}


		void reset(const T& val, const UINT& arr_idx = 0)
		{
			ParentGrid_t& parent = this->parent();
			for (const VecDi& pos_child : this->list(arr_idx))
			{
				ChildGrid_t& child = parent(pos_child);

				for (const VecDi& pos : child.list(arr_idx))
					child(pos) = val;

				child.list(arr_idx).clear();
			}

			this->list(arr_idx).clear();
		}
	};


	template <UINT D, UINT P, UINT N>
	using PosMappedPartitionedGridBase = PartitionedGrid<
		Eigen::Matrix<UINT, N, 1>, D, P, PosArrayMappedGrid<D, N>
	>;

	template <UINT D, UINT P, UINT N>
	class PosMappedPartitionedGrid
	: public PosMappedPartitionedGridBase<D, P, N>
	{
	public:
		typedef PosArrayMappedGrid<D, N>			LookupGrid_t;
		typedef typename LookupGrid_t::PosArray		PosArray;
	protected:
		typedef PosMappedPartitionedGridBase<D, P, N>	Base_t;
		typedef typename Base_t::ChildGrid_t			ChildGrid_t;
		typedef typename ChildGrid_t::VecDu				VecDu;
		typedef typename ChildGrid_t::VecDi				VecDi;

		LookupGrid_t									m_grid_lookup;
	public:
		typedef typename Base_t::ParentGrid_t	ParentGrid_t;


		PosMappedPartitionedGrid () : Base_t(), m_grid_lookup()
		{}


		PosMappedPartitionedGrid (
			const VecDu& dims, const VecDi& offset, const FLOAT& delta = 1
		) : Base_t(), m_grid_lookup()
		{
			init(dims, offset, delta);
		}


		void init (
			const VecDu& dims, const VecDi& offset = VecDi::Zero(),
			const FLOAT& delta = 1
		) {
			Base_t::init(dims, offset, delta);
			lookup().init(dims, offset, delta);
		}


		LookupGrid_t& lookup()
		{
			return m_grid_lookup;
		}


		const LookupGrid_t& lookup() const
		{
			return m_grid_lookup;
		}


		UINT add(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_to_partn(pos);
			ParentGrid_t& parent = this->parent();
			LookupGrid_t& lookup = this->lookup();
			ChildGrid_t& child = parent(pos_child);

			lookup.add(pos_child, arr_idx);
			return child.add(pos, arr_idx);
		}


		void reset(const UINT& arr_idx = 0)
		{
			ParentGrid_t& parent = this->parent();
			LookupGrid_t& lookup = this->lookup();
			for (const VecDi& child_pos : lookup.list(arr_idx))
			{
				parent(child_pos).reset(arr_idx);
			}
			lookup.reset(arr_idx);
		}

		void remove(const VecDi& pos, const UINT& arr_idx = 0)
		{
			const VecDi& pos_child = this->pos_to_partn(pos);
			ParentGrid_t& parent = this->parent();
			LookupGrid_t& lookup = this->lookup();
			ChildGrid_t& child = parent(pos_child);
			child.remove(pos, arr_idx);
			if (child.list(arr_idx).size() == 0)
				lookup.remove(pos_child, arr_idx);
		}
	};

}
#endif /* PARTITIONEDGRID_HPP_ */
