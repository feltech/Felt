#ifndef INCLUDE_FELT_POLYGRID_HPP_

#define INCLUDE_FELT_POLYGRID_HPP_

#include <Felt/MultiLookupGrid.hpp>
#include "PartitionedGrid.hpp"
#include "Surface.hpp"
#include "Poly.hpp"

namespace felt
{
	/**
	 * Container for a grid of Poly objects polygonising a spatially partitioned
	 * signed distance grid.
	 *
	 * Each Poly object polygonises a single spatial partition.
	 * 
	 * @tparam D number of dimensions of the surface and thus grid.
	 * @tparam P Poly compatible class for storing the polygonisations of
	 * spatial partitions.
	 */
	template <class Derived>
	class PolyGridBase : public GridBase < PolyGridBase<Derived> >
	{
	public:
		/// Polygonisation of a single surface spatial partition.
		using DerivedType = typename GridTraits<Derived>::ThisType;
		using PolyLeaf = typename GridTraits<Derived>::LeafType;
		static const UINT Dims = GridTraits<Derived>::Dims;

		using VecDu = typename PolyLeaf::VecDu;
		using VecDi = typename PolyLeaf::VecDi;

		using Base = GridBase< PolyGridBase<Derived> >;

		/// Standard 3-layer signed-distance surface (either 2D or 3D).
		using PolySurface = Surface<Dims, 3>;

		/// MultiLookup grid to track partitions containing zero-layer points.
		using PolyChanges = MultiLookupGrid<Dims>;
		
	protected:
		/// MultiLookup grid to track partitions containing zero-layer points.
		PolyChanges		m_grid_changes;

	public:
		~PolyGridBase ()
		{}
		
		/**
		 * Initialise an empty grid of Poly objects.
		 */
		PolyGridBase () : Base(), m_grid_changes()
		{}

		/**
		 * Construct a grid of polygonisations to fit size of given Surface.
		 * 
		 * @param surface
		 */
		PolyGridBase (const PolySurface& surface) : Base(), m_grid_changes()
		{
			this->init(surface);
		}

		/**
		 * Initialise a grid of polygonisations to fit size of given Surface.
		 *
		 * @param surface
		 */
		void init(const PolySurface& surface)
		{
			Base::init(
				surface.isogrid().children().size(), surface.isogrid().children().offset(),
				PolyLeaf()
			);
			m_grid_changes.init(
				surface.isogrid().size(), surface.isogrid().offset()
			);
			for (const VecDi& pos_child : surface.isogrid().children())
				nself->init_child(
					pos_child, surface.isogrid().children().get(pos_child).size(),
					surface.isogrid().children().get(pos_child).offset()
				);
		}

		/**
		 * Initialise a single polygonisation of a spatial partition.
		 *
		 * Override in subclasses for derived Poly classes.
		 *
		 * @param pos_child
		 */
		void init_child(
			const VecDi& pos_child_, const VecDu& size_, const VecDi& offset_
		) {
			// Add a one-element border to account for partition overlap.
			this->get(pos_child_).init(
				size_ + VecDu::Constant(2), offset_ - VecDi::Constant(1)
			);
		}

		/**
		 * Get the grid of tracked changes to the surface vs. last 
		 * polygonisation.
		 * 
		 * @return grid of changes.
		 */
		const PolyChanges& changes () const
		{
			return m_grid_changes;
		}

		/**
		 * Notify of an update to the surface in order to track changes.
		 *
		 * This should be called whenever the surface is updated to ensure
		 * that eventual re-polygonisation only needs to update those spatial
		 * partitions that have actually changed.
		 *
		 * @param surface
		 */
		void notify(const PolySurface& surface)
		{
			for (INT layer_id = surface.LAYER_MIN; layer_id <= surface.LAYER_MAX; layer_id++)
			{
				for (
					const VecDi& pos_child
					: surface.delta().children().list(surface.layer_idx(layer_id))
				) {
					this->notify(surface, pos_child);
				}
			}

			for (const VecDi& pos_child : surface.status_change().children().list())
			{
				m_grid_changes.add(pos_child);
			}
		}

		/**
		 * Notify that a given spatial partition has been updated.
		 *
		 * @param surface
		 * @param pos_child
		 */
		void notify(const PolySurface& surface, const VecDi& pos_child)
		{
			bool is_active = this->get(pos_child).spx().size() > 0;

			if (!is_active)
			{
				is_active = surface.status_change().children().is_active(pos_child);

				for (
					UINT layer_idx = 0; layer_idx < surface.isogrid().NUM_LISTS && !is_active;
					layer_idx++
				) {
					is_active = surface.isogrid().children().is_active(pos_child, layer_idx);
				}
			}

			if (is_active)
				m_grid_changes.add(pos_child);
			else
				m_grid_changes.remove(pos_child);
		}


		/**
		 * (Re-)Polygonise spatial partitions that have been marked as changed.
		 *
		 * @param surface
		 */
		void poly_cubes (const PolySurface& surface)
		{

			#pragma omp parallel for
			for (
				UINT child_idx = 0;
				child_idx < m_grid_changes.list().size(); child_idx++
			) {
				const VecDi& pos_child = (
					m_grid_changes.list()[child_idx]
				);

				PolyLeaf& leaf = this->get(pos_child);
				leaf.reset();

				for (
					INT layer_id = surface.LAYER_MIN;
					layer_id <= surface.LAYER_MAX; layer_id++
				) {
					for (
						const VecDi& pos : surface.layer(pos_child, layer_id)
					) {
						leaf.spx(pos, surface.isogrid());
					}
				}
			}
		}

		/**
		 * Reset change tracking grid ready for next polygonisation.
		 *
		 * This is separated from poly_cubes() so that derived classes can
		 * access the changes() before it is reset (e.g. to update the GPU).
		 */
		void update_end ()
		{
			m_grid_changes.reset();
		}

		/**
		 * Reset all polygonisations and changes.
		 */
		void reset ()
		{
			for (const VecDi& pos_child : *this)
				this->get(pos_child).reset();
			m_grid_changes.reset();
		}

		/**
		 * Reset and polygonise the whole surface.
		 *
		 * @param surface
		 */
		void surf (const PolySurface& surface)
		{
			reset();

			for (
				UINT layer_idx = 0; layer_idx < surface.NUM_LAYERS; layer_idx++
			) {
				for (
					const VecDi& pos_child
					: surface.isogrid().children().list(layer_idx)
				) {
					m_grid_changes.add(pos_child);
				}
			}

			poly_cubes(surface);
		}
	};

	template <UINT D>
	class PolyGrid : public PolyGridBase<PolyGrid<D> >
	{
	public:
		using Base = PolyGridBase<PolyGrid<D> >;
		using Base::PolyGridBase;
	};


	/**
	 * Traits for PolyGridBase.
	 *
	 * Just pass along the traits for the derived class.
	 */
	template <class Derived>
	struct GridTraits< PolyGridBase<Derived> > : GridTraits<Derived> {};


	template <UINT D>
	struct GridTraits< PolyGrid<D> >
	{
		using ThisType = PolyGrid<D>;
		using LeafType = Poly<D>;
		static const UINT Dims = D;
	};

}

#endif /* INCLUDE_FELT_POLYGRID_HPP_ */
