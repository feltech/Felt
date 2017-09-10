#ifndef INCLUDE_FELT_POLYS_HPP_
#define INCLUDE_FELT_POLYS_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/PartitionedMixin.hpp>
#include <Felt/Impl/Mixin/PolyMixin.hpp>
#include <Felt/Impl/Poly.hpp>

namespace Felt
{

/**
 * Polygonisation for a spatially partitioned level set surface.
 *
 * Holds child `Poly::Single` objects that are dynamically created, updated and destroyed as the
 * surface changes.
 *
 * Call `notify` each time the surface is updated to keep track of spatial partitions that need
 * (re)polygonising.
 *
 * Alternatively call `invalidate` to mark the whole isogrid for (re)polygonisation.
 *
 * Call `march` to go through tracked changes, updating the polygonisation of flagged spatial
 * partitions.
 *
 * After each `march`, call `changes` to get the position indices of partitions that were updated.
 */
template <class TSurface>
class Polys : private Impl::Mixin::Partitioned::Children< Polys<TSurface> >
{
private:
	using This = Polys<TSurface>;
	using Traits = Impl::Traits<This>;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<This>;

	/// Surface type.
	using Surface = typename Traits::Surface;
public:
	/// Child grid type.
	using Child = typename Traits::Child;
private:
	/// Isogrid to (partially) polygonise.
	using IsoGrid = typename Traits::IsoGrid;
	/// Spatial partition type this poly will be responsible for.
	using IsoChild = typename IsoGrid::Child;
	/// Lookup grid to track partitions containing zero-layer points.
	using ChangesGrid = Impl::Lookup::SingleListSingleIdx<Impl::Traits<IsoGrid>::t_dims>;

	Surface const* m_psurface;
	std::unique_ptr<ChangesGrid> m_pgrid_update_pending;
	std::unique_ptr<ChangesGrid> m_pgrid_update_done;

public:
	using ChildrenImpl::children;

	Polys(const TSurface& surface_) :
		ChildrenImpl{
			surface_.isogrid().size(), surface_.isogrid().offset(), surface_.isogrid().child_size(),
			Child(surface_.isogrid())
		},
		m_psurface{&surface_},
		m_pgrid_update_pending{
			std::make_unique<ChangesGrid>(
				surface_.isogrid().children().size(), surface_.isogrid().children().offset()
			)
		},
		m_pgrid_update_done{
			std::make_unique<ChangesGrid>(
				surface_.isogrid().children().size(), surface_.isogrid().children().offset()
			)
		}
	{
		// Bind child polygonisation to isogrid child.
		for (
			PosIdx pos_idx_child = 0; pos_idx_child < this->children().data().size();
			pos_idx_child++
		) {
			this->children().get(pos_idx_child).bind(
				surface_.isogrid().children().get(pos_idx_child).lookup()
			);
		}
	}

	/**
	 * Notify of an update to the surface in order to track changes.
	 *
	 * This should be called whenever the surface is updated to ensure that eventual
	 * re-polygonisation only needs to update those spatial partitions that have actually changed.
	 */
	void notify()
	{
		static const TupleIdx num_lists = m_psurface->isogrid().children().lookup().num_lists;

		// Cycle outermost bands of delta update spatial partitions. We have three cases:
		// * Partition is currently polygonised, and since its in the delta grid it needs updating.
		// * Partition is not currently polygonised, but isogrid is tracking it, in which case it
		//   needs polygonising.
		// * Partition is not currently polygonised and isogrid is no longer tracking, so we no
		//   longer need to remember it.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1) {
			for (const PosIdx pos_idx_child : m_psurface->delta(layer_idx))
			{
				bool is_active = this->children().get(pos_idx_child).is_active();

				for (
					TupleIdx layer_idx = 0; !is_active && layer_idx <= num_lists;
					layer_idx += num_lists - 1
				) {
					is_active = m_psurface->isogrid().children().lookup().is_tracked(
						pos_idx_child, layer_idx
					);
				}

				if (is_active)
					m_pgrid_update_pending->track(pos_idx_child);
				else
					m_pgrid_update_pending->untrack(pos_idx_child);
			}
		}

		// Cycle outermost status change lists, where a child may need resetting.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1)
		{
			for (
				const PosIdx pos_idx_child : m_psurface->status_change(layer_idx)
			) {
				if (this->children().get(pos_idx_child).is_active())
					m_pgrid_update_pending->track(pos_idx_child);
			}
		}
	}

	/**
	 * Repolygonise partitions marked as changed since last polygonisation.
	 */
	void march()
	{
		#pragma omp parallel for
		for (ListIdx list_idx = 0; list_idx < m_pgrid_update_pending->list().size(); list_idx++)
		{
			const PosIdx pos_idx_child = m_pgrid_update_pending->list()[list_idx];

			Child& child = this->children().get(pos_idx_child);

			// If isogrid child partition is active, then polygonise it.
			if (m_psurface->isogrid().children().get(pos_idx_child).is_active())
			{
				if (child.is_active())
					child.reset();
				else
					child.activate();
				child.march();

			} else {
				// Otherwise isogrid child partition has become inactive, so destroy the poly child.
				if (child.is_active())
					child.deactivate();
			}
		}

		std::swap(m_pgrid_update_pending, m_pgrid_update_done);
		m_pgrid_update_pending->reset();
	}


	/**
	 * Add all active poly childs and isogrid childs to change tracking for (re)polygonisation.
	 */
	void invalidate()
	{
		static const TupleIdx num_lists = m_psurface->isogrid().children().lookup().num_lists;

		// Remove pending changes, we're about to reconstruct the list.
		m_pgrid_update_pending->reset();
		// Flag curently active Poly::Single childs for re-polygonisation (or deactivation).
		for (const PosIdx pos_idx_child : this->children().lookup().list())
			m_pgrid_update_pending->track(pos_idx_child);

		// Flag active outer-layer isogrid childs for re-polygonisation.
		for (TupleIdx layer_idx = 0; layer_idx <= num_lists; layer_idx += num_lists - 1)
		{
			for (
				const PosIdx pos_idx_child :
				m_psurface->isogrid().children().lookup().list(layer_idx)
			) {
				m_pgrid_update_pending->track(pos_idx_child);
			}
		}
	}

	/**
	 * Get list of partitions that were updated.
	 *
	 * @return list of position indices of partitions that were repolygonised in the last `march`.
	 */
	const PosIdxList& changes() const
	{
		return m_pgrid_update_done->list();
	}
};

namespace Impl
{
/**
 * Traits for Polys.
 *
 * @tparam Surface surface type to polygonise.
 */
template <class TSurface>
struct Traits< Polys<TSurface> >
{
	/// Type of surface to polygonise.
	using Surface = TSurface;
	/// Isogrid type that the surface wraps.
	using IsoGrid = typename Surface::IsoGrid;
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<IsoGrid>::t_dims;
	/// Child poly type to polygonise a single spatial partition.
	using Child = Impl::Poly::Single<IsoGrid>;
	/// Children grid type to store and track active child polys.
	using Children = Impl::Tracked::SingleListSingleIdxByRef<Child, t_dims>;
};
} // Impl.

} // Felt.

#endif /* INCLUDE_FELT_POLYS_HPP_ */
