#ifndef Surface_hpp
#define Surface_hpp
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Partitioned.hpp>
#include <Felt/Impl/Util.hpp>

#include <vector>
#include <functional>
#include <limits>
#include <iostream>
#include <omp.h>
#include <iostream>

#include <boost/archive/binary_iarchive.hpp>
#include <boost/archive/binary_oarchive.hpp>
#include <boost/math/constants/constants.hpp>
#include <boost/math/special_functions/round.hpp>

#ifndef FELT_SURFACE_OMP_MIN_CHUNK_SIZE
/**
 * Minimum number of (active) spatial partitions required before enabling OpenMP loop parallelism.
 *
 * OpenMP loop parallelisation has signficant overhead (a few milliseconds), so the amount of work
 * done by each thread must be enough to warrant this overhead.
 */
#define FELT_SURFACE_OMP_MIN_CHUNK_SIZE 32
#endif

/// Transform args to a char* string.
#define FELT_STR(args) #args
/// The OpenMP parallel loop command.
#define FELT_PARALLEL_FOR_CMD(num, ...)\
	FELT_STR(omp parallel for if(num >= FELT_SURFACE_OMP_MIN_CHUNK_SIZE) __VA_ARGS__)
/// The #pragma containing the OpenMP parallel loop command
#define FELT_PARALLEL_FOR(num, ...) _Pragma(FELT_PARALLEL_FOR_CMD(num, __VA_ARGS__))


namespace Felt
{
/// Layer ID (in -L, ..., +L)
using LayerId = int;

/**
 * Arbitrarily-dimensional sparse-field spatially partitioned level set surface.
 *
 * @tparam D the number of dimensions of the surface.
 * @tparam L the number of narrow band layers surrounding the zero-layer (surface).
 */
template <Dim D, LayerId L>
class Surface
{
private:
	using This = Surface<D, L>;
	/// Dimensions of the underlying isogrid.
	static constexpr Dim s_dims	= D;
	/// Furthest layer from the zero-layer on the inside of the volume.
	static constexpr LayerId s_layer_min	= -L;
	/// Furthest layer from the zero-layer on the outside of the volume.
	static constexpr LayerId s_layer_max	= L;
	/// Value to indicate a "layer" outside of the volume.
	static constexpr LayerId s_outside		= s_layer_max + 1;
	/// Value to indicate a "layer" inside the volume.
	static constexpr LayerId s_inside		= s_layer_min - 1;
	/// Total number of layers.
	static constexpr LayerId s_num_layers = 2*L+1;
	/// A tiny number used for error margin when raycasting.
	static constexpr Distance TINY = 0.00001f;
	/**
	 * A delta isogrid update grid with active (non-zero) grid points tracked.
	 */
	using DeltaIsoGrid = Impl::Partitioned::Tracked::Simple<Distance, D, s_num_layers>;
public:
	/**
	 * A level set embedding isogrid grid, with active grid points (the narrow band) tracked.
	 */
	using IsoGrid = Impl::Partitioned::Tracked::Numeric<Distance, D, s_num_layers>;
private:
	/**
	 * D-dimensional integer vector.
	 */
	using VecDi = Felt::VecDi<D>;
	/**
	 * D-dimensional float vector.
	 */
	using VecDf = Felt::VecDf<D>;
	/**
	 * Grid to track positions that require an update.
	 */
	using AffectedLookupGrid = Impl::Partitioned::Lookup<D, s_num_layers>;

	/**
	 * Grid tracking locations that are to be moved to another narrow band layer.
	 *
	 * The tracking list index encodes the "from" layer and the value in the grid encodes the
	 * "to" layer.
	 */
	using StatusChangeGrid = Impl::Partitioned::Tracked::Simple<LayerId, D, s_num_layers>;

	/**
	 * Structure to store raycast intermediate results.
	 */
	struct ChildHit
	{
		ChildHit(const VecDf& pos_intersect_, PosIdx pos_idx_child_)
		: pos_intersect(pos_intersect_), pos_idx_child(pos_idx_child_)
		{}
		VecDf pos_intersect;
		PosIdx pos_idx_child;
	};
	/**
	 * The main level set embedding isogrid.
	 */
	IsoGrid				m_grid_isogrid;
	/**
	 * The delta isogrid update grid.
	 *
	 * Used to allow for asynchronous updating.
	 */
	DeltaIsoGrid 		m_grid_delta;
	/**
	 * The (spatially partitioned) status change list.
	 *
	 * The list is appended to when a point in the narrow band moves from
	 * one layer to another.
	 */
	StatusChangeGrid	m_grid_status_change;
	/**
	 * Grid for preventing duplicates when doing neighbourhood queries.
	 */
	AffectedLookupGrid 	m_grid_affected;
	AffectedLookupGrid 	m_grid_affected_buffer;

public:
	/// D-dimensional hyperplane type (using Eigen library), for raycasting.
	using Plane = Eigen::Hyperplane<Distance, D>;
	/// D-dimensional parameterised line, for raycasting.
	using Line = Eigen::ParametrizedLine<Distance, D>;

	/**
	 * Store approximate statistics of the number of spatial partitions using up memory.
	 */
	struct Stats
	{
		/// Number of isogrid partitions currently active.
		ListIdx active_isogrid_partitions;
		/// Number of delta isogrid partitions currently active.
		ListIdx active_delta_partitions;
	};

	/**
	 * Get vector representing a raycast miss.
	 */
	static const VecDf ray_miss;

	/**
	 * Construct a level set embedding of size size.
	 *
	 * All points will be marked as outside the surface (i.e. no surface).
	 *
	 * @param size_ size of the isogrid.
	 * @param size_partition_ size of each spatial partition of the isogrid.
	 */
	Surface (const VecDi& size_, const VecDi& size_partition_ = VecDi::Constant(8)) :
		// Configure isogrid embedding, initialising to all outside values.
		m_grid_isogrid(size_, offset(size_), size_partition_, s_outside),
		// Configure delta isogrid embedding, initialising to zero delta.
		m_grid_delta(size_, offset(size_), size_partition_, 0),
		// Configure status change partitioned lists, use "outside" value as convenient "null".
		m_grid_status_change(size_, offset(size_), size_partition_, s_outside),
		// Configure de-dupe grid for neighbourhood queries.
		m_grid_affected(size_, offset(size_), size_partition_),
		m_grid_affected_buffer(size_, offset(size_), size_partition_)
	{}

	/**
	 * Save isogrid to given output stream.
	 *
	 * @param output_stream_ stream to save to.
	 */
	void save(std::ostream& output_stream_) const
	{
		m_grid_isogrid.write(output_stream_);
	}

	/**
	 * Load isogrid from given input stream and construct surface.
	 *
	 * @param input_stream_ stream to load from.
	 *
	 * @return new Surface instance.
	 */
	static This load(std::istream& input_stream_)
	{
		IsoGrid isogrid{IsoGrid::read(input_stream_)};
		return This{std::move(isogrid)};
	}


	/**
	 * Create a single singularity seed point in the isogrid grid.
	 *
	 * NOTE: does not handle overwriting of points currently already on the
	 * surface/in the volume.
	 *
	 * @param pos_centre
	 */
	void seed (const VecDi& pos_centre_)
	{
		// Width of seed.
		const VecDi& vec_width = VecDi::Constant(NodeIdx(s_layer_max));

		// Min and max positions affected by placing seed point.
		const VecDi& pos_min = pos_centre_ - vec_width;
		const VecDi& pos_max = pos_centre_ + vec_width;

		// Get vector size of window formed by pos_min and pos_max.
		const VecDi& pos_window_size = pos_max - pos_min + VecDi::Constant(1); //+1 for zero coord.

		// Calculate number of grid points to be cycled through within window.
		const PosIdx pos_idx_max = PosIdx(pos_window_size.prod());

		// Cycle through each point in window.
		for (PosIdx pos_idx = 0; pos_idx <= pos_idx_max; pos_idx++)
		{
			// Calculate vector position from integer index, using utility function.
			VecDi pos = Felt::index<D>(pos_idx, pos_window_size);
			// Translate position into isogrid grid space.
			pos += pos_min;
			// Calculate vector distance from this position to seed centre.
			const VecDi& vec_dist = pos - pos_centre_;
			// Sum of absolute distance along each axis == city-block distance.
			const Distance dist = Distance(vec_dist.template lpNorm<1>());
			const LayerId layer_id_pos = layer_id(dist);
			// Check distance indicates that this point is within the narrow band.
			if (inside_band(layer_id_pos))
			{
				// Append point to a narrow band layer (if applicable).
				m_grid_isogrid.track(dist, pos, layer_idx(layer_id_pos));
			}
		}
	}

	/**
	 * Perform a full update of the narrow band.
	 *
	 * Lambda function passed will be given the position to process and
	 * a reference to the isogrid grid, and is expected to return delta isogrid to
	 * apply.
	 *
	 * Each spatial partition is processed in parallel.
	 *
	 * @param fn_ (pos, isogrid) -> float
	 */
	template <typename Fn>
	void update(Fn&& fn_)
	{
		update_start();

		// We are iterating over the entire zero-layer, so assume the delta grid should track
		// all active partitions in the main isogrid.
		m_grid_delta.track_children(m_grid_isogrid);

		const PosIdxList& pos_idxs_children = m_grid_isogrid.children().lookup().list(layer_idx(0));

		FELT_PARALLEL_FOR(pos_idxs_children.size(),)
		for (ListIdx list_idx = 0; list_idx < pos_idxs_children.size(); list_idx++)
		{
			const PosIdx pos_idx_child = pos_idxs_children[list_idx];
			typename IsoGrid::Child& isochild = m_grid_isogrid.children().get(pos_idx_child);

			for (const PosIdx pos_idx_leaf : isochild.lookup().list(layer_idx(0)))
			{
				const VecDi& pos_leaf = isochild.index(pos_idx_leaf);
				const Distance dist_delta = callback(
					fn_, pos_leaf, const_cast<const IsoGrid&>(m_grid_isogrid)
				);

				#ifdef FELT_DEBUG_ENABLED
				if (std::abs(dist_delta) > 1.0f)
				{
					std::stringstream strs;
					strs << "Zero layer update value out of bounds: " << Felt::format(pos_leaf)
						<< " with value " << dist_delta;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				#endif

				// Disallow expansion to edge of grid.
				m_grid_delta.children().get(pos_idx_child).track(
					dist_delta, pos_idx_leaf, layer_idx(0)
				);
			}
		}
		update_end_global();
	}

	/**
	 * Perform a bounded update of the narrow band.
	 *
	 * Lambda function passed will be given the position to process and a reference to the isogrid
	 * grid, and is expected to return delta isogrid to apply.
	 *
	 * Each spatial partition is processed in parallel.
	 *
	 * @param pos_leaf_lower_ region selection from.
	 * @param pos_leaf_upper_ region selection to.
	 * @param fn_ (pos, isogrid) -> float.
	 */
	template <typename Fn>
	void update(const VecDi& pos_leaf_lower_, const VecDi& pos_leaf_upper_, Fn&& fn_)
	{
		static const VecDi& one = VecDi::Constant(1);
		// Upper and lower bounds of the grid, inclusive.
		const VecDi& pos_grid_lower = m_grid_isogrid.offset();
		const VecDi& pos_grid_upper =
			m_grid_isogrid.offset() + m_grid_isogrid.size();
		// Child partitions containing upper and lower bounds of grid.
		const VecDi& pos_grid_child_lower = m_grid_isogrid.pos_child(pos_grid_lower);
		const VecDi& pos_grid_child_upper = m_grid_isogrid.pos_child(pos_grid_upper - one);
		// Partition containing lower point of bounding box, bounded by grid.
		const VecDi& pos_child_lower =
			pos_grid_child_lower.cwiseMax(m_grid_isogrid.pos_child(pos_leaf_lower_));
		// Partition containing upper point of bounding box, bounded by grid.
		const VecDi& pos_child_upper =
			pos_grid_child_upper.cwiseMin(m_grid_isogrid.pos_child(pos_leaf_upper_));
		// Size of bounding box at partition level.
		const VecDi& child_bounding_box_size = pos_child_upper - pos_child_lower + one;
		// Upper bound of leaf (1 more than upper point), bounded by grid..
		const VecDi& pos_leaf_upper_bound = pos_grid_upper.cwiseMin(pos_leaf_upper_ + one);
		// Upper index of bounding box.
		const PosIdx child_idx_bound = PosIdx(child_bounding_box_size.prod());
		// Clear previous update.
		update_start();
		// Parallel loop through spatial partitions.
		FELT_PARALLEL_FOR(child_idx_bound)
		for (PosIdx child_idx = 0; child_idx < child_idx_bound; child_idx++)
		{
			// Get spatial partition position.
			const VecDi& pos_child_without_offset = Felt::index<D>(
				child_idx, child_bounding_box_size
			);
			const VecDi& pos_child = pos_child_without_offset + pos_child_lower;
			const PosIdx pos_idx_child = m_grid_isogrid.children().index(pos_child);
			typename IsoGrid::Child& child = m_grid_isogrid.children().get(pos_idx_child);

			// Loop all zero-layer points within this partition.
			for (const PosIdx pos_idx_leaf : child.lookup().list(layer_idx(0)))
			{
				const VecDi& pos_leaf = child.index(pos_idx_leaf);
				// Skip zero-layer points not within finer-grained bounding box.
				if (Felt::inside(pos_leaf, pos_leaf_lower_, pos_leaf_upper_bound))
				{
					const Distance dist_delta = callback(
						fn_, pos_leaf, const_cast<const IsoGrid&>(m_grid_isogrid)
					);

					#ifdef FELT_DEBUG_ENABLED

					if (m_grid_delta.children().get(pos_idx_child).get(pos_idx_leaf) != 0)
					{
						std::stringstream strs;
						strs << "Delta is not zero: " << Felt::format(pos_leaf) << " with delta " <<
							m_grid_delta.children().get(pos_idx_child).get(pos_idx_leaf);
						std::string str = strs.str();
						throw std::domain_error(str);
					}
					if (std::abs(dist_delta) > 1.0f)
					{
						std::stringstream strs;
						strs << "Zero layer update value out of bounds: " << Felt::format(pos_leaf)
						<< " with value " << dist_delta;
						std::string str = strs.str();
						throw std::domain_error(str);
					}

					#endif

					// Update delta isogrid.
					m_grid_delta.track(dist_delta, pos_idx_child, pos_idx_leaf, layer_idx(0));
				}
			}
		}
		// Apply delta to isogrid.
		update_end();
	}

	/**
	 * Reset delta isogrid to zero and clear update lists.
	 *
	 * @snippet test_Surface.cpp Simple delta isogrid update
	 */
	void update_start ()
	{
		m_grid_delta.reset(m_grid_isogrid);
		m_grid_affected.reset(m_grid_isogrid);
		m_grid_affected_buffer.reset(m_grid_isogrid);
		m_grid_status_change.reset(m_grid_isogrid);
	}

	/**
	 * Update zero layer then update distance transform for affected points in each layer.
	 */
	void update_end ()
	{
		// Get points in outer layers that are affected by changes in zero-layer.
		calc_affected();

		m_grid_isogrid.track_children(m_grid_affected);
		m_grid_delta.track_children(m_grid_affected);

		// Update the zero layer, applying delta to isogrid.
		update_zero_layer(&m_grid_affected_buffer);
		AffectedLookupGrid* plookup = &m_grid_affected;
		AffectedLookupGrid* plookup_buffer = &m_grid_affected_buffer;
		bool is_status_changed = update_distance(plookup, plookup_buffer);
		flush_status_change();
		expand_narrow_band();

		while (is_status_changed)
		{
			plookup->reset(m_grid_isogrid);
			std::swap(plookup, plookup_buffer);
			m_grid_status_change.reset(m_grid_affected);
			is_status_changed = update_distance(plookup, plookup_buffer);
			flush_status_change();
			expand_narrow_band();
		}

	}


	/**
	 * Update delta isogrid grid, tracking to tracking list if not already tracked.
	 *
	 * @snippet test_Surface.cpp Simple delta isogrid update
	 *
	 * @param pos position to update.
	 * @param val value to set at given position.
	 */
	void delta (const VecDi& pos_, Distance val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

		const LayerId layer_id_new = layer_id(val_);
		if (layer_id_new != 0 && layer_id_new != 1 && layer_id_new != -1)
		{
			std::stringstream strs;
			strs << "Delta update value out of bounds. Attempted to update position " <<
				Felt::format(pos_) << " by " << val_ << " would give a layer of "
				<< layer_id_new << ", which is too much of a jump";
			std::string str = strs.str();
			throw std::domain_error(str);
		}

		#endif

		m_grid_delta.track(val_, pos_, layer_idx(0));
	}

	/**
	 * Cast a ray to the zero layer.
	 *
	 * @param pos_origin_ originating point of ray
	 * @param dir_ normal vector in direction of ray
	 * @return zero curve hit location or NULL_POS.
	 */
	VecDf ray(const VecDf& pos_origin_, const VecDf& dir_) const
	{

		using ChildHits = std::vector<ChildHit>;

		// If ray is cast from within isogrid grid, first check child grid containing origin point.
		if (m_grid_isogrid.inside(pos_origin_))
		{
			const VecDf& pos_hit = ray(
				pos_origin_, dir_,
				m_grid_isogrid.children().get(
					m_grid_isogrid.pos_child(pos_origin_.template cast<NodeIdx>())
				)
			);
			if (pos_hit != ray_miss)
				return pos_hit;
		}

		// Ray to test against.
		Line line(pos_origin_, dir_);

		// Tracking list for child grids that are hit.
		ChildHits child_hits;

		// Cycle each axis, casting ray to child grid planes marching away from origin.
		for (Dim dim = 0; dim < dir_.size(); dim++)
		{
			// Direction +/-1 along this axis.
			Distance dir_dim = Distance(sgn(dir_(dim)));
			if (dir_dim == 0)
				continue;

			// Get next child plane along this axis.
			Distance pos_plane_dim = round_to_next(
				dim, dir_dim, pos_origin_(dim), m_grid_isogrid.child_size()
			);

			// Construct vector with elements not on this axis at zero.
			VecDf pos_plane = VecDf::Constant(0);
			pos_plane(dim) = pos_plane_dim;

			// If the zero point on this plane is not within the grid, then jump to max/min point
			// on isogrid grid.
			if (!m_grid_isogrid.inside(pos_plane))
			{
				Distance pos_grid_dim;
				// If casting in -'ve direction, get maximum extent.
				if (dir_dim == -1)
				{
					pos_grid_dim = Distance(
						m_grid_isogrid.offset()(dim) + m_grid_isogrid.size()(dim)
					);
					if (pos_plane_dim < pos_grid_dim)
						continue;
				}
				// Else if casting in +'ve direction, get minimum extent.
				else
				{
					pos_grid_dim = Distance(m_grid_isogrid.offset()(dim));
					if (pos_plane_dim > pos_grid_dim)
						continue;
				}
				// Reset plane position to max/min extent as calculated above.
				pos_plane(dim) = pos_grid_dim;
			}

			// Plane normal is opposite to ray direction.
			VecDf normal = VecDf::Constant(0);
			normal(dim) = -dir_dim;

			// Cast ray to plane and track any child grids hit on the way to tracking list.
			// If child size is not a factor of grid size then this first cast could be to outside
			// the grid.  So cannot quit early here and must try next child.
			ray_check_track_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim));

			// Round up/down to next child, in case we started at inexact modulo of child grid size
			// (i.e. when isogrid grid size is not integer multiple of child grid size).
			pos_plane_dim = round_to_next(
				dim, dir_dim, pos_plane(dim), m_grid_isogrid.child_size()
			);
			// If rounding produced a different plane, then cast to that plane, and potentially track
			// child grid to tracking list.
			if (pos_plane_dim != pos_plane(dim))
			{
				pos_plane(dim) = pos_plane_dim;
				if (!ray_check_track_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
					continue;
			}

			// Keep marching along planes, casting ray to each and tracking any candidate child
			// grids to the tracking list.
			const Distance child_size_dim = Distance(m_grid_isogrid.child_size()(dim));
			while (true)
			{
				pos_plane(dim) += dir_dim * child_size_dim;
				if (
					!ray_check_track_child(
						child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)
					)
				) {
					break;
				}
			}
		}

		// untrack children whose raycast hit point lies outside themselves.  Possible if a child
		// has been allowed to be tracked in order for the loop above to continue on to the next
		// child, via dir vs. offset/size check in ray_check_track_child().
//		child_hits.erase(std::untrack_if(
//			child_hits.begin(), child_hits.end(),
//			[this](const ChildHit& child_hit) {
//				return !this->m_grid_isogrid.children()
//					.get(child_hit.pos_child)
//					.inside(child_hit.pos_intersect);
//			}
//		), child_hits.end());

		// Sort candidate child grids in distance order from front to back.
		std::sort(
			child_hits.begin(), child_hits.end(),
			[&pos_origin_](const ChildHit& a, const ChildHit& b) -> bool {
				const VecDf& dist_a = (a.pos_intersect - pos_origin_);
				const VecDf& dist_b = (b.pos_intersect - pos_origin_);
				return dist_a.squaredNorm() < dist_b.squaredNorm();
			}
		);
		// untrack any duplicate child grids from the sorted list (i.e. where ray intersects
		// precisely at the interesction of two or more planes).
		child_hits.erase(std::unique(
			child_hits.begin(), child_hits.end(),
			[](const ChildHit& a, const ChildHit& b) -> bool {
				return a.pos_idx_child == b.pos_idx_child;
			}
		), child_hits.end());

		// For each candidate child, cast ray through until the zero-curve is hit.
		for (const ChildHit& child_hit : child_hits)
		{
			const VecDf& pos_hit = ray(
				child_hit.pos_intersect, dir_,
				m_grid_isogrid.children().get(child_hit.pos_idx_child)
			);

			if (pos_hit != ray_miss)
				return pos_hit;
		}

		return ray_miss;
	}

	/**
	 * Gather statistics about the current state of the surface.
	 *
	 * @return statistics object with values for the number of active spatial partitions.
	 */
	Stats stats() const
	{
		return Stats{
			num_active_partitions(m_grid_isogrid),
			num_active_partitions(m_grid_delta)
		};
	}

	/**
	 * Get list of spatial partitions where iso values were updated in last update.
	 *
	 * @param layer_idx_ index of layer list to get.
	 * @return list of position indices.
	 */
	const PosIdxList& delta(const TupleIdx layer_idx_) const
	{
		return m_grid_delta.children().lookup().list(layer_idx_);
	}

	/**
	 * Get list of spatial partitions where layer status change occurred in last update.
	 *
	 * @param layer_idx_ index of layer list to get.
	 * @return list of position indices.
	 */
	const PosIdxList& status_change(const TupleIdx layer_idx_) const
	{
		return m_grid_status_change.children().lookup().list(layer_idx_);
	}

	/**
	 * Get reference to isogrid grid.
	 *
	 * @return signed distance isogrid embedding the level set surface.
	 */
	const IsoGrid& isogrid () const
	{
		return m_grid_isogrid;
	}

	/**
	 * Get narrow band layer index of ID for indexing into arrays.
	 *
	 * Internally we can of course not have negative array indices, so must convert to an
	 * index first.
	 *
	 * @param id narrow band layer ID.
	 * @return index to use to get given layer in narrow band array.
	 */
	static constexpr TupleIdx layer_idx(const LayerId id_)
	{
		return TupleIdx(id_ + LayerId(s_num_layers) / 2);
	}

	/**
	 * Check if child spatial partition contains zero-layer points.
	 *
	 * @param pos_idx_child position index of child partition.
	 * @return true if zero-curve cuts through partition, false otherwise.
	 */
	bool is_intersected(const PosIdx pos_idx_child)
	{
		return bool(m_grid_isogrid.children().get(pos_idx_child).list(layer_idx(0)).size());
	}

	/**
	 * Get points affected by most recent update, in the form of a lookup grid.
	 *
	 * @return the affected lookup grid.
	 */
	const AffectedLookupGrid& affected() const
	{
		return m_grid_affected;
	}

///////////////////////////////////////////////////////////////////////////////////////////////////
// Begin private
private:

	/**
	 * Explicitly deleted default constructor.
	 */
	Surface () = delete;

	/**
	 * Construct surface from isogrid.
	 *
	 * @param isogrid_
	 */
	Surface (IsoGrid&& isogrid_) :
		// Configure isogrid embedding, initialising to all outside values.
		m_grid_isogrid{std::move(isogrid_)},
		// Configure delta isogrid embedding, initialising to zero delta.
		m_grid_delta{
			m_grid_isogrid.size(), m_grid_isogrid.offset(), m_grid_isogrid.child_size(), 0
		},
		// Configure status change partitioned lists, use "outside" value as convenient "null".
		m_grid_status_change{
			m_grid_isogrid.size(), m_grid_isogrid.offset(), m_grid_isogrid.child_size(), s_outside
		},
		// Configure de-dupe grid for neighbourhood queries.
		m_grid_affected{
			m_grid_isogrid.size(), m_grid_isogrid.offset(), m_grid_isogrid.child_size()
		},
		m_grid_affected_buffer{
			m_grid_isogrid.size(), m_grid_isogrid.offset(), m_grid_isogrid.child_size()
		}
	{}


	/**
	 * Update zero layer then update distance transform for all points in all layers.
	 *
	 * Must do two (or more) passes to clean up singularities and deal with corner cases where
	 * outer layer distances are not quite consistent after layer promotion.
	 */
	void update_end_global ()
	{
		update_zero_layer(&m_grid_affected_buffer);

		bool is_status_changed = update_distance(&m_grid_isogrid, &m_grid_affected_buffer);
		flush_status_change();
		expand_narrow_band();

		if (is_status_changed)
		{
			AffectedLookupGrid* plookup = &m_grid_affected_buffer;
			AffectedLookupGrid* plookup_buffer = &m_grid_affected;

			while (is_status_changed)
			{
				m_grid_status_change.reset(m_grid_isogrid);
				is_status_changed = update_distance(plookup, plookup_buffer);
				flush_status_change();
				expand_narrow_band();
				plookup->reset(m_grid_isogrid);
				std::swap(plookup, plookup_buffer);
			}
		}
	}

	/**
	 * Find all outer layer points who's distance transform is affected by modified zero-layer
	 * points.
	 *
	 * Must do two (or more) passes to clean up singularities and deal with corner cases where
	 * outer layer distances are not quite consistent after layer promotion.
	 *
	 * TODO: several options for optimisation of removing duplicates:
	 * - Use a boolean flag grid to construct a de-duped vector (used here).
	 * - Check std::vector in neighs() using std::find to prevent tracking a duplicate in the
	 *   first place.
	 * - Use std::vector sort, unique, erase.
	 * - Use a std::unordered_set with a suitable hashing function.
	 *
	 * @param apos
	 */
	void calc_affected()
	{
		const TupleIdx layer_idx_zero = layer_idx(0);

		// Loop over delta isogrid modified zero-layer points tracking to
		// tracking grid.

		// Loop spatial partitions of delta for zero-layer.
		for (
			const PosIdx pos_idx_child
			: m_grid_delta.children().lookup().list(layer_idx_zero))
		{
			// Loop leaf grid nodes with spatial partition
			for (
				const PosIdx pos_idx_leaf
				: m_grid_delta.children().get(pos_idx_child).list(layer_idx_zero)
			)
				// track zero-layer point to tracking grid.
				m_grid_affected.track(pos_idx_child, pos_idx_leaf, layer_idx_zero);
		}

		// Arrays to store first and last element in tracking list within each spatial partition of
		// tracking grid.
		std::array<std::vector<ListIdx>, std::size_t(s_num_layers)> aidx_first_neigh;
		std::array<std::vector<ListIdx>, std::size_t(s_num_layers)> aidx_last_neigh;

		// Loop round L times, searching outward for affected outer layer grid nodes.
		for (LayerId udist = 1; udist <= s_layer_max; udist++)
		{
			// Reset the first and last element indices for each spatial partition in each layer.

			for (LayerId layer_id = s_layer_min; layer_id <= s_layer_max; layer_id++)
			{
				const TupleIdx layer_idx = this->layer_idx(layer_id);
				// Get number of spatial partitions for this layer.
				const ListIdx num_childs =
					m_grid_affected.children().lookup().list(layer_idx).size();
				// Resize spatial partition index lists for this layer to
				// to include any newly tracked partitions.
				aidx_last_neigh[ListIdx(layer_idx)].resize(num_childs);
				// Will initialise to zero, so no further work needed for
				// these new indices giving the start of the range.
				aidx_first_neigh[ListIdx(layer_idx)].resize(num_childs);
				// The final index needs to be copied from the current size
				// of each spatial partition, so loop over partitions,
				// copying their size into the respective last index list.
				for (
					ListIdx list_idx_child = 0; list_idx_child < num_childs; list_idx_child++
				) {
					// Get position of this spatial partition in parent
					// lookup grid.
					const PosIdx pos_idx_child =
						m_grid_affected.children().lookup().list(layer_idx)[list_idx_child];
					// Copy number of active grid nodes for this partition
					// into relevant index in the list.
					aidx_last_neigh[ListIdx(layer_idx)][list_idx_child] =
						m_grid_affected.children().get(pos_idx_child).list(layer_idx).size();
				}
			}

			// Loop each layer finding the affected outer layer points
			// for each partition using the start and end points cached
			// above.

			for (LayerId layer_id = s_layer_min; layer_id <= s_layer_max; layer_id++)
			{
				const TupleIdx layer_idx = this->layer_idx(layer_id);

				// Loop over spatial partitions, ignoring newly tracked ones
				// since we're using the cached spatial partition list as
				// the end of the range.
				for (
					ListIdx idx_child = 0;
					idx_child < aidx_first_neigh[ListIdx(layer_idx)].size();
					idx_child++
				) {
					// Get position of this spatial partition in this layer.
					const PosIdx pos_idx_child = (
						m_grid_affected.children().lookup().list(layer_idx)[idx_child]
					);

					// Loop over leaf grid nodes within this spatial partition, using the cached
					// start and end indices, so that newly tracked points are skipped.
					for (
						ListIdx idx_neigh = aidx_first_neigh[ListIdx(layer_idx)][idx_child];
						idx_neigh < aidx_last_neigh[ListIdx(layer_idx)][idx_child];
						idx_neigh++
					) {
						using AffectedChild = typename AffectedLookupGrid::Child;
						const AffectedChild& child =  m_grid_affected.children().get(pos_idx_child);
						// This leaf grid node is the centre to search about.
						const PosIdx pos_idx_centre = child.list(layer_idx)[idx_neigh];
						const VecDi pos_centre = child.index(pos_idx_centre);

						// Cycle neighbours and record them if they are within the narrow band.
						m_grid_isogrid.neighs(
							pos_centre,
							[this](const VecDi& pos_neigh) {
								// Calculate layer of this neighbouring point from the isogrid grid.
								const LayerId layer_id_neigh = this->layer_id(pos_neigh);
								// If the calculated layer lies within the narrow band, then we
								// want to track it.
								if (inside_band(layer_id_neigh))
								{
									// track the neighbour point to the
									// tracking grid. Will reject
									// duplicates.
									m_grid_affected.track(
										pos_neigh, this->layer_idx(layer_id_neigh)
									);
								}
							}
						); // End neighbourhood query.
					} // End for leaf grid node.
				} // End for spatial partition.
			} // End for layer.

			// Now we've found the neighbours of the pre-existing grid
			// points in the tracking list, we want to skip these on the
			// next loop, so set the start index for each partition of
			// each layer to the previous end index.

			for (LayerId layer_id = s_layer_min; layer_id <= s_layer_max; layer_id++)
			{
				const ListIdx layer_idx = ListIdx(this->layer_idx(layer_id));
				// Loop over spatial partitions.
				for (ListIdx idx = 0; idx < aidx_first_neigh[layer_idx].size(); idx++)
					// Set first index in spatial partition's tracking list
					// to be the previous last index, so we start there on
					// next loop around.
					aidx_first_neigh[layer_idx][idx] = aidx_last_neigh[layer_idx][idx];
			}
		}
	} // End calc_affected.

	/**
	 * Loop over the zero-layer of the delta grid and apply it to the isogrid.
	 *
	 * @param plookup_buffer_ grid to store points that need a layer update.
	 */
	void update_zero_layer (AffectedLookupGrid* plookup_buffer_)
	{
		using DeltaChild = typename DeltaIsoGrid::Child;
		using IsoChild = typename IsoGrid::Child;

		const TupleIdx layer_idx_zero = layer_idx(0);
		const PosIdxList& pos_idxs_children = m_grid_delta.children().lookup().list(layer_idx_zero);
		const ListIdx num_childs = pos_idxs_children.size();

		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_buffer_))
		for (ListIdx list_idx_child = 0; list_idx_child < num_childs; list_idx_child++)
		{
			const PosIdx pos_idx_child = pos_idxs_children[list_idx_child];
			DeltaChild& delta_child = m_grid_delta.children().get(pos_idx_child);
			IsoChild& isogrid_child = m_grid_isogrid.children().get(pos_idx_child);

			for (const PosIdx pos_idx_leaf : delta_child.lookup().list(layer_idx_zero))
			{
				const Distance iso_prev = isogrid_child.get(pos_idx_leaf);
				const Distance iso_delta = delta_child.get(pos_idx_leaf);
				const Distance iso_new = iso_prev + iso_delta;
				const LayerId layer_id_new = layer_id(iso_new);

				#ifdef FELT_DEBUG_ENABLED

				const LayerId layer_id_old = layer_id(iso_prev);

				if (layer_id_old != 0)
				{
					std::stringstream strs;
					strs << "Zero layer update attempted at non-zero layer point " <<
						Felt::format(isogrid_child.index(pos_idx_leaf)) << ": " << iso_prev <<
						" + " << iso_delta << " = " << iso_new;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				if (
					std::abs(layer_id_new) != 0
					&& std::abs(layer_id_new) != 1
					&& std::abs(layer_id_new) != -1
				) {
					std::stringstream strs;
					strs << "Zero layer update out of bounds.  Attempting to change value at" <<
						Felt::format(isogrid_child.index(pos_idx_leaf)) << " to " << iso_new <<
						" would give a layer of " << layer_id_new <<
						", which is too much of a jump";
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				// Update value in grid with new signed distance.
				isogrid_child.set(pos_idx_leaf, iso_new);
				// Potentially track to status change, if narrow band layer has changed.
				status_change(pos_idx_child, pos_idx_leaf, 0, layer_id_new, plookup_buffer_);
			}
		}
	}

	/**
	 * Perform distance transform on narrow band layers, from centre working outwards.
	 *
	 * @param plookup_ either isogrid for full update, or affected grid for local updates.
	 * @param plookup_buffer_ buffer to store grid points that need a layer update and may need
	 * additional iterations to converge.
	 */
	template <typename Grid>
	bool update_distance(
		const Grid* plookup_, AffectedLookupGrid* plookup_buffer_
	) {
		bool is_status_changed = false;

		// Update distance transform for inner layers of the narrow band.
		for (LayerId layer_id = -1; layer_id >= s_layer_min; layer_id--)
			is_status_changed |= update_distance(layer_id, -1, plookup_, plookup_buffer_);

		// Update distance transform for outer layers of the narrow band.
		for (LayerId layer_id = 1; layer_id <= s_layer_max; layer_id++)
			is_status_changed |= update_distance(layer_id, 1, plookup_, plookup_buffer_);

		return is_status_changed;
	}

	/**
	 * Update distance transform for all points in given layer.
	 *
	 * @param layer_id_ layer to update
	 * @param side_ side of narrow band (+1 for outside and -1 for inside the volume).
	 * @param lookup_ either isogrid for full update, or affected grid for local updates.
	 * @param plookup_buffer_ grid to store points that need to change layers.
	 * @return true if any point status changed, false otherwise.
	 */
	template <typename Grid>
	bool update_distance(
		const LayerId layer_id_, const LayerId side_,
		const Grid* plookup_, AffectedLookupGrid* plookup_buffer_
	) {
		using IsoChild = typename IsoGrid::Child;
		using DeltaIsoChild = typename DeltaIsoGrid::Child;

		bool is_status_changed = false;

		const TupleIdx layer_idx = this->layer_idx(layer_id_);
		const PosIdxList& pos_idxs_children = plookup_->children().lookup().list(layer_idx);
		const ListIdx num_childs = pos_idxs_children.size();

		// First pass: calculate distance and add to delta isogrid.
		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_, plookup_buffer_))
		for (ListIdx list_idx = 0; list_idx < num_childs; list_idx++)
		{
			// Child spatial partition position
			const ListIdx pos_idx_child = pos_idxs_children[list_idx];

			const PosIdxList& apos_leafs =
				plookup_->children().get(pos_idx_child).list(layer_idx);

			// Calculate distance of every point in this layer to the zero
			// layer, and store in delta isogrid grid.
			// Delta isogrid grid is used to allow for asynchronous updates, that
			// is, to prevent neighbouring points affecting the distance
			// transform.
			for (const PosIdx pos_idx_leaf : apos_leafs)
			{
				// Distance from this position to zero layer.
				const Distance dist = distance(pos_idx_child, pos_idx_leaf, side_);

				#ifdef FELT_DEBUG_ENABLED

				const LayerId layer_id_new = layer_id(dist);

				if (
					layer_id_new != layer_id_ &&
					layer_id_new != layer_id_ + 1 &&
					layer_id_new != layer_id_ - 1
				) {
					const VecDi& pos =
						m_grid_isogrid.children().get(pos_idx_child).index(pos_idx_leaf);
					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds." << std::endl <<
						str_neighs(pos) << " distance of " << dist <<
						", which is too much of a jump";

					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				// Update delta isogrid grid.
				m_grid_delta.children().get(pos_idx_child).track(dist, pos_idx_leaf, layer_idx);
			}
		}

		// Second pass: apply distance to isogrid and update status change lists.
		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_, plookup_buffer_))
		for (ListIdx list_idx = 0; list_idx < num_childs; list_idx++)
		{
			const PosIdx pos_idx_child = pos_idxs_children[list_idx];

			IsoChild& grid_isogrid_child = m_grid_isogrid.children().get(pos_idx_child);
			DeltaIsoChild& grid_delta_child = m_grid_delta.children().get(pos_idx_child);

			const PosIdxList& pos_idxs_leafs =
				plookup_->children().get(pos_idx_child).list(layer_idx);

			// Update distance in isogrid from delta isogrid and append any points that
			// move out of their layer to a status change list.
			for (const PosIdx pos_idx_leaf : pos_idxs_leafs)
			{
				// Distance calculated above.
				const Distance dist = grid_delta_child.get(pos_idx_leaf);
				const LayerId layer_id_new = layer_id(dist);

				#ifdef FELT_DEBUG_ENABLED

				if (
					layer_id_new != layer_id_ &&
					layer_id_new != layer_id_ + 1 &&
					layer_id_new != layer_id_ - 1
				) {
					const VecDi& pos =
						m_grid_isogrid.children().get(pos_idx_child).index(pos_idx_leaf);
					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds. Attempting to" <<
						" move " << Felt::format(pos) << " in layer " << layer_id_ << " to a" <<
						" distance of " << dist << " would result in a layer of " <<
						layer_id_new << ", which is too much of a jump" ;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				grid_isogrid_child.set(pos_idx_leaf, dist);
				is_status_changed |= status_change(
					pos_idx_child, pos_idx_leaf, layer_id_, layer_id_new, plookup_buffer_
				);
			}

		}

		return is_status_changed;

	} // End update_distance.

	/**
	 * Potentially track a point to the status change list to eventually be moved from one layer to
	 * another.
	 *
	 * @param pos_ position in grid to check.
	 * @param layer_id_from_ layer moving from.
	 * @param layer_id_to_ layer moving to.
	 * @param plookup_buffer_ grid to store points that need status change, for potential additional
	 * processing.
	 * @return true if status change is needed for this position, false otherwise.
	 */
	bool status_change (
		const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_, const LayerId layer_id_from_,
		const LayerId layer_id_to_, AffectedLookupGrid* plookup_buffer_
	) {
		using StatusChangeChild = typename StatusChangeGrid::Child;

		if (layer_id_from_ == layer_id_to_)
			return false;

		#ifdef FELT_DEBUG_ENABLED
		m_grid_isogrid.children().assert_pos_idx_bounds(pos_idx_child_, "status_change child: ");
		m_grid_isogrid.children().get(pos_idx_child_).assert_pos_idx_bounds(
			pos_idx_leaf_, "status_change leaf: "
		);
		#endif

		StatusChangeChild& child = m_grid_status_change.children().get(pos_idx_child_);
		LayerId layer_id_to = child.get(pos_idx_leaf_);

		// If the position is already marked for status change, this must be a subsequent loop
		// around `converge_distance`, and this leaf position is "jumping" more than one layer.
		if (layer_id_to != s_outside)
		{
			child.set(pos_idx_leaf_, layer_id_to_);
		}
		else
		{
			m_grid_status_change.track(
				layer_id_to_, pos_idx_child_, pos_idx_leaf_, layer_idx(layer_id_from_)
			);
		}

		// Keep a record of points that undergo a status change, since these will have to have
		// their distance calculated again, and again, until they no longer status change.  This
		// is required for collapsing regions where the surface should disappear because the
		// zero-layer is gone.
		if (inside_band(layer_id_to_))
			plookup_buffer_->track(pos_idx_child_, pos_idx_leaf_, layer_idx(layer_id_to_));

		return true;
	}

	/**
	 * Loop through the status change lists moving the points from one layer to another.
	 */
	void flush_status_change ()
	{
		using StatusChangeChild = typename StatusChangeGrid::Child;

		for (LayerId layer_id_from = s_layer_min; layer_id_from <= s_layer_max; layer_id_from++)
		{
			const TupleIdx layer_idx_from = layer_idx(layer_id_from);

			const PosIdxList& pos_idxs_children =
				m_grid_status_change.children().lookup().list(layer_idx_from);
			const ListIdx num_childs = pos_idxs_children.size();

			FELT_PARALLEL_FOR(num_childs, firstprivate(layer_id_from))
			for (
				ListIdx list_idx_child = 0; list_idx_child < num_childs;
				list_idx_child++
			) {
				const ListIdx pos_idx_child = pos_idxs_children[list_idx_child];
				const StatusChangeChild& child =
					m_grid_status_change.children().get(pos_idx_child);

				for (const PosIdx pos_idx_leaf : child.lookup().list(layer_idx_from))
				{
					const LayerId layer_id_to = child.get(pos_idx_leaf);
					const TupleIdx layer_idx_to = layer_idx(layer_id_to);

					if (inside_band(layer_id_to))
					{
						#ifdef FELT_DEBUG_ENABLED

						if (child.lookup().list(layer_idx_from).size() == 0)
						{
							std::stringstream strs;
							strs << "Layer empty when attempting to move " <<
								Felt::format(child.index(pos_idx_leaf)) <<
								" from layer " << layer_id_from << " to layer " << layer_id_to <<
								" in partition " <<
								Felt::format(m_grid_isogrid.children().index(pos_idx_child)) <<
								" = " << Felt::format(child.offset()) << "-" <<
								Felt::format(child.offset() + child.size());
							std::string str = strs.str();
							throw std::domain_error(str);
						}

						#endif

						m_grid_isogrid.retrack(
							pos_idx_child, pos_idx_leaf, layer_idx_from, layer_idx_to
						);
					}
					else
					{
						// Remove from tracking, potentially deactivating child and setting it's
						// background value (distance) to value of target layer id.
						m_grid_isogrid.untrack(
							Distance(layer_id_to), pos_idx_child, pos_idx_leaf, layer_idx_from
						);
					}
				}
			}
		}
	}

	/**
	 * Add new points to the narrow band when expanding/contracting.
	 */
	void expand_narrow_band()
	{
		using StatusChangeChild = typename StatusChangeGrid::Child;

		// Cycle innermost layer and outermost layer.
		for (LayerId layer_id : std::array<LayerId, 2>{{s_layer_min, s_layer_max}})
		{
			const TupleIdx layer_idx = this->layer_idx(layer_id);

			const PosIdxList& apos_children =
				m_grid_status_change.children().lookup().list(layer_idx);

			const LayerId side = sgn(layer_id);

			// TODO: not thread-safe - neighbouring nodes at edge of child grid could be in another
			// child, so `get` and `track` calls below are unsafe. Perhaps worth a lock?
//			#pragma omp parallel for
			for (
				ListIdx list_idx_child = 0; list_idx_child < apos_children.size();
				list_idx_child++
			) {
				const PosIdx pos_idx_child = apos_children[list_idx_child];
				const StatusChangeChild& child = m_grid_status_change.children().get(pos_idx_child);

				for (const PosIdx pos_idx : child.lookup().list(layer_idx))
				{
					// If not expanding/contracting, then nothing to do here.
					if (child.get(pos_idx) != s_layer_max * side - side)
						continue;

					const VecDi& pos = child.index(pos_idx);

					// Cycle over neighbours of this outer layer point.
					m_grid_isogrid.neighs(
						pos,
						[this, layer_idx, side, layer_id, &pos](const VecDi& pos_neigh_) {
							if (not m_grid_isogrid.inside(pos_neigh_))
								return;

							Distance distance_neigh = m_grid_isogrid.get(pos_neigh_);
							const LayerId layer_id_from = this->layer_id(distance_neigh);

							// Only track if neighbouring point is not already within the narrow
							// band.
							if (inside_band(layer_id_from))
							{
								#ifdef FELT_DEBUG_ENABLED

								const ListIdx lookup_idx = m_grid_isogrid.children().get(
										m_grid_isogrid.pos_child(pos_neigh_)
									).lookup().get(pos_neigh_);

								if (lookup_idx == null_idx)
								{
									std::stringstream sstr;
									sstr << "pos not tracked but should be: " <<
										str_pos(pos_neigh_);
									std::string str = sstr.str();
									throw std::domain_error(str);
								}

								#endif

								return;
							}

							// Calculate updated distance of this neighbour to the zero curve.
							distance_neigh = distance(pos_neigh_, distance_neigh, Distance(side));

							#ifdef FELT_DEBUG_ENABLED

							const LayerId layer_id_to = this->layer_id(distance_neigh);

							if (layer_id_to != layer_id)
							{
								std::stringstream strs;
								strs << "Neighbour is further away than expected." << std::endl <<
									"pos:" << std::endl <<
									"  " << str_pos(pos) << std::endl
									<< "Neigh:" << std::endl <<
									"  " << str_pos(pos_neigh_) << std::endl <<
									"Calculated distance " << distance_neigh <<
									" would give a layer of " << layer_id_to <<
									" when we expect a layer of " << layer_id;
								std::string str = strs.str();
								throw std::domain_error(str);
							}

							if (layer_id_to != s_layer_min && layer_id_to != s_layer_max)
							{
								std::stringstream strs;
								strs << "Attempting to track " << Felt::format(pos_neigh_) <<
									" to the narrow band but the distance is " << distance_neigh <<
									" which would give a layer of " << layer_id_to;
								std::string str = strs.str();
								throw std::domain_error(str);
							}

							#endif

							this->m_grid_isogrid.track(distance_neigh, pos_neigh_, layer_idx);
						}
					);
				} // End for pos_idx.
			}
		}
	}

	/**
	 * Calculate city-block distance from position to zero curve.
	 *
	 * @param pos_ target position.
	 * @param side_ side of narrow band (+/-1).
	 * @return city-block distance from pos_ to the zero curve.
	 */
	Distance distance (
		const PosIdx pos_idx_child_, const PosIdx pos_idx_leaf_, const LayerId side_
	) const {

		const typename IsoGrid::Child& child = m_grid_isogrid.children().get(pos_idx_child_);
		// Position vector from position index.
		VecDi pos = child.index(pos_idx_leaf_);
		// Current distance recorded at this point.
		const Distance dist = child.get(pos_idx_leaf_);
		// Direction away from surface along surface normal in distance units.
		const Distance dir = Distance(side_);

		return distance(std::move(pos), dist, dir);
	}

	/**
	 * Calculate city-block distance from position to zero curve.
	 *
	 * @param pos_ position who's distance to zero-curve we wish to find.
	 * @param dist_ current signed distance at this position.
	 * @param dir_ direction away from surface along surface normal.
	 * @return
	 */
	Distance distance (VecDi pos_, Distance dist_, const Distance dir_) const
	{
		#ifdef FELT_DEBUG_ENABLED
		const VecDi pos_original(pos_);
		#endif
		// Transform to unsigned distance.
		dist_ *= dir_;

		// Get neighbouring point that is next closest to the zero-layer.
		m_grid_isogrid.neighs(
			std::move(pos_),
			[this, &dist_, dir_](const VecDi& pos_neigh_) {
				if (not m_grid_isogrid.inside(pos_neigh_))
					return;
				const Distance dist_neigh = m_grid_isogrid.get(pos_neigh_);
				// Check absolute value of this neighbour is less than nearest point.
				// Multiplying by `dir` has two effects:
				// - It has the same effect as abs() for points on the same side of the band.
				// - It ensures points on the opposite side are negative, so that < comparison
				//   prefers those points, which is good because we're interested in the neighbour
				//   in the *direction* of the zero-curve.
				if (dist_neigh*dir_ < dist_) {
					dist_ = dist_neigh*dir_;
				}
			}
		);

		// This point's distance is then the distance of the closest neighbour +/-1, depending
		// which side of the band we are looking at. So first transform back into signed distance
		// then add +/-1.
		const Distance dist_neigh = dist_*dir_;
		dist_ = dist_neigh + dir_;

		#ifdef FELT_DEBUG_ENABLED
		const LayerId layer_id_pos = layer_id(pos_original);
		const LayerId layer_id_neigh = layer_id(dist_neigh);

		if (
			std::abs(layer_id_pos) < std::abs(layer_id_neigh) &&
			sgn(layer_id_pos) == sgn(layer_id_neigh)
		) {
			std::stringstream sstr;
			sstr << "Neighbour closest to zero curve is further away than source position: " <<
				Felt::format(pos_original) << " at " << m_grid_isogrid.get(pos_original) <<
				" is closer than " << Felt::format(pos_) << " at " << m_grid_isogrid.get(pos_) <<
				" but should not be";
			std::string str = sstr.str();
			throw std::domain_error(str);
		}

		#endif

		return dist_;
	}

	/**
	 * Cast a ray to the zero layer within a given child grid.
	 *
	 * @param pos_origin
	 * @param dir
	 * @return NULL_POS if no hit, otherwise interpolated position on zero
	 * curve.
	 */
	VecDf ray(
		VecDf pos_sample, const VecDf& dir, const typename IsoGrid::Child& child
	) const {
		using Line = Eigen::ParametrizedLine<Distance, D>;

		const Line line_leaf(pos_sample, dir);
		Distance t_leaf = 0;

//		std::cerr << "Child: (" << child.offset().transpose() << ") + (" <<
//			child.size().transpose() << ")" << std::endl;

		while (child.inside(pos_sample))
		{
			const LayerId layer_id = this->layer_id(pos_sample);

//			std::cerr << layer_id << std::endl;

			if (abs(layer_id) == 0)
			{
				VecDf normal = m_grid_isogrid.grad(pos_sample);

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
				if (std::isnan(normal.normalized()[0]))
				{
					std::stringstream strs;
					strs << "ERROR: raycast isosurface gradient normal is NaN " <<
						"when normalising " << Felt::format(normal);
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				#endif

				normal.normalize();

//				std::cerr << "Normal: " << normal.transpose() << std::endl;

				if (normal.dot(dir) < 0)
				{
					static constexpr int MAX_CONVERGE_STEPS = 100;
					int num_converge_steps = 0;
					Distance dist = 0;
					for (; num_converge_steps < MAX_CONVERGE_STEPS; num_converge_steps++)
					{
						dist = m_grid_isogrid.interp(pos_sample);

						pos_sample -= normal*dist;

						if (!m_grid_isogrid.inside(pos_sample))
							return ray_miss;

						if (std::abs(dist) <= TINY || normal.dot(dir) >= 0)
							break;

						normal = m_grid_isogrid.grad(pos_sample);
						normal.normalize();
					}

					#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
					if (num_converge_steps == MAX_CONVERGE_STEPS)
					{
//						std::cerr << "WARNING: raycast failed to get to distance < "
//							<< TINY << ".\n"
//							<< "dir(" << Felt::format(dir) << ");"
//							<< " normal(" << Felt::format(normal) << ");"
//							<< " sample(" << Felt::format(pos_sample) << ")"
//							<< " at dist " << dist;
					}
					#endif

					return pos_sample;
				}
			}

			t_leaf += 0.5f;

			pos_sample = line_leaf.pointAt(t_leaf);
		} // End while inside child grid.

		return ray_miss;
	}

	/**
	 * Cast ray to plane, get child at that point, and track to list if it contains zero-curve.
	 *
	 * @param child_hits
	 * @param line
	 * @param plane
	 * @return
	 */
	bool ray_check_track_child(
		std::vector<ChildHit>& child_hits, const Line& line, const Plane& plane
	) const {
		const VecDf& pos_intersect =
			line.intersectionPoint(plane) + line.direction() * Distance(TINY);

		const VecDi& size = m_grid_isogrid.size();
		const VecDi& offset = m_grid_isogrid.offset();
		const VecDf& dir = line.direction();

		for (Dim i = 0; i < dir.size(); i++)
		{
			if (
				(dir(i) > 0 && pos_intersect(i) > Distance(size(i))) ||
				(dir(i) < 0 && pos_intersect(i) < Distance(offset(i)))
			)
				return false;
		}

		if (!m_grid_isogrid.inside(pos_intersect))
			return true;

		const VecDi& pos_floor = pos_intersect.array().floor().matrix().template cast<NodeIdx>();
		const PosIdx pos_idx_child = m_grid_isogrid.pos_idx_child(pos_floor);

		if (
			layer(pos_idx_child, 0).size() ||
			layer(pos_idx_child, 1).size() || layer(pos_idx_child, -1).size()
		) {
			child_hits.push_back(ChildHit(pos_intersect, pos_idx_child));
		}
		return true;
	}

	/**
	 * Along a given dimension at given position, round up or down to border of next partition.
	 *
	 * @param dim dimension e.g. 0 for x, 1 for y.
	 * @param dir direction along given dimension, eiether -1 or 1.
	 * @param pos position to round from
	 * @param part_size size of partition in space to round to.
	 * @return position of plane
	 */
	Distance round_to_next(
		const Dim dim, const Distance dir, const Distance pos, const VecDi part_size
	) const {
		// Real-valued child pos translated to [0, 2*childsize) space.
		Distance pos_plane_dim =
			(pos - Distance(m_grid_isogrid.offset()(dim))) / Distance(part_size(dim));

		// Round to next child en route in [0, 2*childsize) space.
		pos_plane_dim = (dir == -1.0f) ?
			std::floor(pos_plane_dim) : std::ceil(pos_plane_dim);
		// Scale back to isogrid grid in [0, 2*fullsize) space.
		pos_plane_dim *= Distance(m_grid_isogrid.child_size()(dim));
		// Translate back to isogrid grid in [-fullsize, fullsize) space.
		pos_plane_dim += Distance(m_grid_isogrid.offset()(dim));

		return pos_plane_dim;
	}


#ifdef FELT_DEBUG_ENABLED

	/**
	 * Stringify a position vector, including information about the isogrid at that point.
	 *
	 * Includes the position vector, spatial partition position, isogrid value, narrow band layer,
	 * and index in the narrow band tracking list (both spatial partition level and leaf level).
	 *
	 * @param pos_ position to stringify
	 * @return string representation of point and the isogrid state at that point.
	 */
	std::string str_pos(const VecDi& pos_) const
	{
		const Distance dist_pos = m_grid_isogrid.get(pos_);
		const LayerId layer_id_pos = layer_id(pos_);
		const VecDi& pos_child = m_grid_isogrid.pos_child(pos_);
		const typename IsoGrid::Child child = m_grid_isogrid.children().get(pos_child);
		const VecDi& pos_child_lower = child.offset();
		const VecDi& pos_child_upper = child.offset() + child.size();
		const Felt::Tuple<ListIdx, s_num_layers>& list_idxs_child =
			m_grid_isogrid.children().lookup().get(pos_child);
		const TupleIdx list_id_pos = layer_idx(layer_id_pos);
		const ListIdx list_idx_pos = child.lookup().get(pos_);

		std::stringstream sstr;
		sstr << Felt::format(pos_) << " ∈ P(" <<
			Felt::format(pos_child) << ") = [" <<
			Felt::format(pos_child_lower) << "," << Felt::format(pos_child_upper) << "] @ " <<
			dist_pos << " ∈ L(" << layer_id_pos << ") @ " << Felt::format(list_idxs_child) <<
			"[" << list_id_pos << "][" << list_idx_pos << "]";
		return sstr.str();
	}


	/**
	 * Call str_pos on all neighbours of given position.
	 *
	 * @param pos_ position to report along with it's neighbours.
	 * @return string representation of the point and its neighbours, with trackitional isogrid info.
	 */
	std::string str_neighs(const VecDi& pos_) const
	{
		std::stringstream sstr;
		sstr << str_pos(pos_) << std::endl << "in:" << std::endl;
		m_grid_isogrid.neighs(pos_, [this, &sstr](const VecDi& pos_neigh) {
			sstr << "    " << str_pos(pos_neigh) << std::endl;
		});
		return sstr.str();
	}
#endif

	/**
	 * Get reference to a single layer of the narrow band at a given
	 * spatial partition.
	 *
	 * @param pos_child_idx_ location of child spatial partition.
	 * @param layer_id_ narrow band layer id.
	 * @return list of positions in given narrow band layer within given child partition.
	 */
	const PosIdxList& layer (const PosIdx pos_child_idx_, const LayerId layer_id_) const
	{
		return m_grid_isogrid.children().get(pos_child_idx_).lookup().list(layer_idx(layer_id_));
	}

	/**
	 * Get narrow band layer id of location in isogrid grid.
	 *
	 * Calculates layer ID based on value in grid (i.e. is the layer the location *should*
	 * be in).
	 *
	 * The position vector type is templated, so that interpolation overloads can be used
	 * if a float vector (e.g. Vec3f) is provided.
	 *
	 * @param pos_ position vector of location
	 * @return integer layer ID that given location should belong to.
	 */
	template <class Pos>
	LayerId layer_id(const Pos& pos_) const
	{
		const LayerId layer_id_pos = layer_id(m_grid_isogrid.get(pos_));

		return layer_id_pos;
	}

	/**
	 * Get narrow band layer ID of given value.
	 *
	 * Rounds value to nearest integer, but with an epsilon to prefer rounding up, to keep
	 * consistent when we have floating point rounding errors.
	 *
	 * @param val value to round to give narrow band layer ID.
	 * @return layer ID that given value should belong to
	 */
	LayerId layer_id(const Distance val_) const
	{
		return boost::math::iround(val_ + std::numeric_limits<Distance>::epsilon());
	}

	/**
	 * Test whether a given value lies (or should lie) within the narrow band or not.
	 *
	 * @param val float or int value to test
	 * @return true if inside the narrow band, false otherwise.
	 */
	template <typename Val>
	bool inside_band (const Val& val) const
	{
		return LayerId(std::abs(val)) <= s_layer_max;
	}

	/**
	 * Offset of isogrid from given size
	 *
	 * Simply minus half the size.
	 */
	static constexpr VecDi offset(const VecDi& size_)
	{
		return -1 * size_ / 2;
	}

	/**
	 * Caclculate an approximation of the number of active spatial partitions for a given grid.
	 *
	 * Is an approximation, since it takes the maximum number of partitions across all
	 * layers, rather than doing a more expensive union.
	 *
	 * @tparam TGrid type of grid to query
	 * @param grid _the grid to query
	 *
	 * @return the max of the number of spatial partitions in each layer.
	 */
	template <class TGrid>
	ListIdx num_active_partitions(const TGrid& grid_) const
	{
		ListIdx num = 0;
		for (TupleIdx layer_idx = 0; layer_idx < s_num_layers; layer_idx++)
		{
			num = std::max(num, grid_.children().lookup().list(layer_idx).size());
		}

		return num;
	}

	/**
	 * Callback wrapper to call an arbitrary functor with a given pos and grid.
	 *
	 * @tparam Fn deduced functor type.
	 * @param fn_  functor.
	 * @param pos_ position in grid.
	 * @param grid_ grid itself.
	 *
	 * @return a Distance value returned by the functor.
	 */
	template <typename Fn>
	std::enable_if_t<
		std::is_same<std::result_of_t<Fn(const VecDi&, const IsoGrid&)>, Distance>::value,
		Distance
	>
	callback(Fn&& fn_, const VecDi& pos_, const IsoGrid& grid_)
	{
		return fn_(pos_, grid_);
	}

	/**
	 * Callback wrapper to call an arbitrary functor with a given pos.
	 *
	 * @tparam Fn deduced functor type.
	 * @param fn_  functor.
	 * @param pos_ position in grid.
	 * @param grid_ grid itself (unused).
	 *
	 * @return a Distance value returned by the functor.
	 */
	template <typename Fn>
	std::enable_if_t<
		std::is_same<std::result_of_t<Fn(const VecDi&)>, Distance>::value,
		Distance
	>
	callback(Fn&& fn_, const VecDi& pos_, const IsoGrid&)
	{
		return fn_(pos_);
	}

	/**
	 * Callback wrapper to call an arbitrary functor.
	 *
	 * @tparam Fn deduced functor type.
	 * @param fn_  functor (unused).
	 * @param pos_ position in grid (unused).
	 * @param grid_ grid itself (unused).
	 *
	 * @return a Distance value returned by the functor.
	 */
	template <typename Fn>
	std::enable_if_t<
		std::is_same<std::result_of_t<Fn()>, Distance>::value,
		Distance
	>
	callback(Fn&& fn_, const VecDi&, const IsoGrid&)
	{
		return fn_();
	}
};

template <Dim D, LayerId L>
const typename Surface<D, L>::VecDf Surface<D, L>::ray_miss =
	VecDf::Constant(std::numeric_limits<Distance>::max());

}

#endif
