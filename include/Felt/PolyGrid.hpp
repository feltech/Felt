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
	 * @tparam D the number of dimensions of the surface and thus grid.
	 */
	template <UINT D>
	class PolyGrid : public TrackedGrid<Poly<D>, D>
	{
	public:
		/// Polygonisation of a single surface spatial partition.
		typedef Poly<D>				PolyLeaf;

		typedef typename PolyLeaf::VecDu	VecDu;
		typedef typename PolyLeaf::VecDi	VecDi;

		typedef TrackedGrid<Poly<D>, D>		Base;

		/// Lookup grid to track partitions containing zero-layer points.
		typedef LookupPartitionedGrid<D>		PolyChanges;
		/// Signed-distance surface (either 2D or 3D).
		template <UINT L> using Surface_t = Surface<D, L>;
	protected:
		/// Lookup grid to track partitions containing zero-layer points.
		PolyChanges		m_grid_changes;

	public:
		
		/**
		 * Construct a grid of polygonisations to fit size of given Surface.
		 * 
		 * @param surface
		 */
		template <UINT L>
		PolyGrid (const Surface_t<L>& surface)
		: 	Base(
				surface.phi().branch().dims(), surface.phi().branch().offset()
			), m_grid_changes(
				surface.phi().dims(), surface.phi().offset(),
				surface.phi().child_dims()
			)
		{
			for (const VecDi& pos_child : surface.phi().branch())
			{
				// Add a one-element border to account for partition overlap.
				this->get(pos_child).init(
					surface.phi().child(pos_child).dims() + VecDu::Constant(2),
					surface.phi().child(pos_child).offset() - VecDi::Constant(1)
				);
			}
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
		template <UINT L>
		void notify(const Surface_t<L>& surface)
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
					Surface<D, L>::StatusChange::layer_idx(0)
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
		template <UINT L>
		void notify(const Surface_t<L>& surface, const VecDi& pos_child)
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
		template <UINT L>
		void poly_cubes(const Surface_t<L>& surface)
		{
			const UINT& zero_layer_idx = surface.layer_idx(0);
			for (const VecDi& pos_child : m_grid_changes.branch().list())
			{
				for (
					const VecDi& pos_centre
					: surface.phi().child(pos_child).list(zero_layer_idx)
				) {
					// Must flag cubes (4x for 2D, 8x for 3D) surrounding the
					// zero-layer point.
					for (const VecDi& pos_offset : Poly<D>::corners)
					{
						const VecDi& pos_corner = pos_centre - (
							pos_offset - Poly<D>::SpxGridPosOffset
						);
						if (!m_grid_changes.inside(pos_corner))
							continue;
						m_grid_changes.add(pos_corner);
					}
				}
			}

			// Second pass, since lookaround search above may have added extra
			// spatial partitions.
			for (const VecDi& pos_child : m_grid_changes.branch().list())
			{
				this->add(pos_child);
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
	};

}

#endif /* INCLUDE_FELT_POLYGRID_HPP_ */
