#ifndef INCLUDE_FELT_POLYGRID_HPP_
#define INCLUDE_FELT_POLYGRID_HPP_

#include "MappedGrid.hpp"
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
	template <UINT D, class P=Poly<D> >
	class PolyGrid : public Grid<P, D>
	{
	public:
		/// Polygonisation of a single surface spatial partition.
		typedef P							PolyLeaf;

		typedef typename PolyLeaf::VecDu	VecDu;
		typedef typename PolyLeaf::VecDi	VecDi;

		typedef Grid<PolyLeaf, D>	Base;

		/// Standard 2-layer signed-distance surface (either 2D or 3D).
		typedef Surface<D>	PolySurface;

		/// Lookup grid to track partitions containing zero-layer points.
		typedef LookupPartitionedGrid<D>		PolyChanges;
	protected:
		/// Lookup grid to track partitions containing zero-layer points.
		PolyChanges		m_grid_changes;

	public:
		virtual ~PolyGrid ()
		{}
		
		/**
		 * Initialise an empty grid of Poly objects.
		 */
		PolyGrid () : Base(), m_grid_changes()
		{}

		/**
		 * Construct a grid of polygonisations to fit size of given Surface.
		 * 
		 * @param surface
		 */
		PolyGrid (const PolySurface& surface) : Base(), m_grid_changes()
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
				surface.phi().branch().dims(), surface.phi().branch().offset()
			);
			m_grid_changes.init(
				surface.phi().dims(), surface.phi().offset(),
				surface.phi().child_dims()
			);
			for (const VecDi& pos_child : surface.phi().branch())
				this->init_child(
					pos_child, surface.phi().child(pos_child).dims(),
					surface.phi().child(pos_child).offset()
				);
		}

		/**
		 * Initialise a single polygonisation of a spatial partition.
		 *
		 * Override in subclasses for derived Poly classes.
		 *
		 * @param pos_child
		 */
		virtual void init_child(
			const VecDi& pos_child_, const VecDu& dims_, const VecDi& offset_
		) {
			// Add a one-element border to account for partition overlap.
			this->get(pos_child_).init(
				dims_ + VecDu::Constant(2), offset_ - VecDi::Constant(1)
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
			// Loop over partitions containing active delta phi grid nodes for
			// zero-layer.
			for (
				const VecDi& pos_child
				: surface.dphi().branch().list(surface.layer_idx(0))
			) {
				this->notify(surface, pos_child);
			}

			for (
				const VecDi& pos_child
				: surface.status_change().branch().list(
					PolySurface::StatusChange::layer_idx(0)
				)
			) {
				this->notify(surface, pos_child);
			}
		}

		/**
		 * Notify that a given spatial partition has been updated.
		 *
		 * If there are zero-layer points in the partition then it will be added
		 * to the tracking list for eventual re-polygonisation. If there are no
		 * zero-layer points then the partition will either be added to the
		 * tracking list for eventual deletion, or removed from the tracking
		 * list, as appropriate.
		 *
		 * @param surface
		 * @param pos_child
		 */
		void notify(const PolySurface& surface, const VecDi& pos_child)
		{
			const UINT& zero_layer_idx = surface.layer_idx(0);
			// Zero-layer cuts through spatial partition.
			if (
				surface.phi().branch().lookup().is_active(
					pos_child, zero_layer_idx
				)
			) {
				m_grid_changes.branch().add(pos_child);
			}
			// Zero-layer no longer cuts through spatial partition.
			else if (this->get(pos_child).spx().size())
			{
				m_grid_changes.branch().add(pos_child);
			}
			// Zero-layer doesn't cut through spatial partition and didn't
			// in the previous polygonisation either.
			else
			{
				m_grid_changes.branch().remove(pos_child);
			}
		}


		/**
		 * (Re-)Polygonise spatial partitions that have been marked as changed.
		 *
		 * @param surface
		 */
		void poly_cubes(const PolySurface& surface)
		{
			const UINT& zero_layer_idx = surface.layer_idx(0);
			const UINT branch_size = m_grid_changes.branch().list().size();
			// NOTE: cannot use range-based loop since spatial partition
			// tracking list can be resized in this loop.
			#pragma omp parallel for
			for (UINT child_idx = 0; child_idx < branch_size; child_idx++)
			{
				// pos_child is not a reference since reallocation can occur.
				VecDi pos_child = m_grid_changes.branch().list()[child_idx];
				for (
					const VecDi& pos_centre
					: surface.phi().child(pos_child).list(zero_layer_idx)
				) {
					// Must flag cubes (4x for 2D, 8x for 3D) surrounding the
					// zero-layer point.
					for (const VecDi& pos_offset : PolyLeaf::corners)
					{
						const VecDi& pos_corner = pos_centre - (
							pos_offset - PolyLeaf::SpxGridPosOffset
						);
						if (!m_grid_changes.inside(pos_corner))
							continue;
						m_grid_changes.add_safe(pos_corner);
					}
				}
			}

			// Second pass, since lookaround search above may have added extra
			// spatial partitions.

			#pragma omp parallel for
			for (
				UINT child_idx = 0;
				child_idx < m_grid_changes.branch().list().size(); child_idx++
			) {
				const VecDi& pos_child = (
					m_grid_changes.branch().list()[child_idx]
				);

				// TODO: not thread-safe, but maybe not even needed? Does this
				// really have to be a TrackedGrid?
//				this->add(pos_child);

				PolyLeaf& leaf = this->get(pos_child);
				leaf.reset();
				for (
					const VecDi& pos_cube
					: m_grid_changes.child(pos_child).list()
				) {
					leaf.spx(pos_cube, surface.phi());
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
	};

}

#endif /* INCLUDE_FELT_POLYGRID_HPP_ */
