#ifndef Impl_Surface_hpp
#define Impl_Surface_hpp

#include <vector>
#include <functional>
#include <limits>
#include <set>
#include <iostream>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/constants/constants.hpp>
#include <eigen3/Eigen/Dense>
#include <omp.h>

#include <Felt/Impl/Partitioned.hpp>
#include <Felt/Impl/Util.hpp>

#ifndef Felt_SURFACE_OMP_CHUNK_SIZE
/**
 * Minimum number of (active) spatial partitions required before enabling OpenMP loop parallelism.
 *
 * OpenMP loop parallelisation has signficant overhead (a few milliseconds), so the amount of work
 * done by each thread must be enough to warrant this overhead.
 */
#define Felt_SURFACE_OMP_MIN_CHUNK_SIZE 32
#endif

/// Transform args to a char* string.
#define FELT_STR(args) #args
/// The OpenMP parallel loop command.
#define FELT_PARALLEL_FOR_CMD(num, args)\
	FELT_STR(omp parallel for if(num >= 32) args)
/// The #pragma containing the OpenMP parallel loop command
#define FELT_PARALLEL_FOR(num, args) _Pragma(FELT_PARALLEL_FOR_CMD(num, args))


namespace Felt
{

/**
 * A n-dimensional sparse-field spatially partitioned level set.
 *
 * @tparam D the number of dimensions of the surface.
 * @tparam L the number of narrow band layers surrounding the zero-level
 * surface.
 */
template <Dim D, UINT L>
class Surface
{
public:
	/// Layer ID (in -L, ..., +L)
	using LayerId = INT;
	/// Furthest layer from the zero-layer on the inside of the volume.
	static constexpr LayerId LAYER_MIN	= -LayerId(L);
	/// Furthest layer from the zero-layer on the outside of the volume.
	static constexpr LayerId LAYER_MAX	= LayerId(L);
	/// Total number of layers.
	static constexpr TupleIdx NUM_LAYERS = 2*L+1;
	/// A tiny number used for error margin when raycasting.
	static constexpr FLOAT TINY = 0.00001f;

	using ThisType = Surface<D, L>;
	/**
	 * A delta isogrid update grid with active (non-zero) grid points tracked.
	 */
	using DeltaIsoGrid = Impl::Partitioned::Tracked::Simple<FLOAT, D, NUM_LAYERS>;
	/**
	 * A level set embedding isogrid grid, with active grid points (the narrow
	 * band) tracked.
	 */
	using IsoGrid = Impl::Partitioned::Tracked::Numeric<FLOAT, D, NUM_LAYERS>;
	/**
	 * D-dimensional unsigned int vector.
	 */
	using VecDu = Felt::VecDu<D>;
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
	using AffectedLookupGrid = Impl::Partitioned::Lookup<D, NUM_LAYERS>;

	/// D-dimensional hyperplane type (using Eigen library), for raycasting.
	using Plane = Eigen::Hyperplane<FLOAT, D>;
	/// D-dimensional parameterised line, for raycasting.
	using Line = Eigen::ParametrizedLine<FLOAT, D>;
	/**
	 * Grid tracking locations that are to be moved to another narrow band layer.
	 *
	 * The tracking list index encodes the "from" layer and the value in the grid encodes the
	 * "to" layer.
	 */
	using StatusChangeGrid = Impl::Partitioned::Tracked::Simple<INT, D, NUM_LAYERS>;


private:
	/**
	 * Structure to store raycast intermediate results.
	 */
	struct ChildHit
	{
		ChildHit(const VecDf& pos_intersect_, const VecDi& pos_child_)
		: pos_intersect(pos_intersect_), pos_child(pos_child_)
		{}
		VecDf pos_intersect;
		VecDi pos_child;
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

	/**
	 * Explicitly defined default constructor.
	 */
	Surface () = delete;

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
		m_grid_isogrid(size_, offset(size_), size_partition_, LAYER_MAX+1),
		// Configure delta isogrid embedding, initialising to zero delta.
		m_grid_delta(size_, offset(size_), size_partition_, 0),
		// Configure status change partitioned lists.
		m_grid_status_change(size_, offset(size_), size_partition_, LAYER_MAX+1),
		// Configure de-dupe grid for neighbourhood queries.
		m_grid_affected(size_, offset(size_), size_partition_),
		m_grid_affected_buffer(size_, offset(size_), size_partition_)
	{}

	/**
	 * Get reference to isogrid grid.
	 *
	 * @return signed distance isogrid embedding the level set surface.
	 */
	IsoGrid& isogrid ()
	{
		return m_grid_isogrid;
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
	 * Get grid of affected narrow band points used during localised update mode.
	 *
	 * @return ephemeral "affected" grid.
	 */
	const AffectedLookupGrid& affected ()
	{
		return m_grid_affected;
	}

	/**
	 * Get the status change grid that flags when a point is moving between narrow band layers
	 *
	 * @return status change grid.
	 */
	const StatusChangeGrid& status_change () const
	{
		return m_grid_status_change;
	}

	/**
	 * Get reference to delta grid of isogrid updates.
	 *
	 * @return grid of delta values to be applied to isogrid.
	 */
	DeltaIsoGrid& delta ()
	{
		return m_grid_delta;
	}

	/**
	 * @copydoc Surface::delta()
	 *
	 * const version.
	 */
	const DeltaIsoGrid& delta () const
	{
		return m_grid_delta;
	}

	/**
	 * Get reference to the active partitions of a given layer of the narrow.
	 *
	 * @return list of positions of active child partitions in the main isogrid.
	 */
	const PosArray& parts () const
	{
		return m_grid_isogrid.children().lookup().list(layer_idx(0));
	}

	/**
	 * Get reference to the zero layer of the narrow band at a given spatial partition.
	 *
	 * @param pos_child_idx_ location of child spatial partition.
	 * @return list of positions in given narrow band layer within given child partition.
	 */
	const PosArray& layer (const PosIdx pos_child_idx_) const
	{
		return m_grid_isogrid.children().get(pos_child_idx_).lookup().list(layer_idx(0));
	}

	/**
	 * Get reference to a single layer of the narrow band at a given
	 * spatial partition.
	 *
	 * @param pos_child_idx_ location of child spatial partition.
	 * @param layer_id_ narrow band layer id.
	 * @return list of positions in given narrow band layer within given child partition.
	 */
	const PosArray& layer (const PosIdx pos_child_idx_, const LayerId layer_id_) const
	{
		return m_grid_isogrid.children().get(pos_child_idx_).lookup().list(layer_idx(layer_id_));
	}

	/**
	 * Update delta isogrid grid, tracking to tracking list if not already tracked.
	 *
	 * @snippet test_Surface.cpp Simple delta isogrid update
	 *
	 * @param pos position to update.
	 * @param val value to set at given position.
	 */
	void delta (const VecDi& pos_, FLOAT val_)
	{
		#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

		const INT newlayer_id = this->layer_id(val_);
		if (newlayer_id != 0 && newlayer_id != 1 && newlayer_id != -1)
		{
			std::stringstream strs;
			strs << "Delta update value out of bounds. Attempted to update position " <<
				Felt::format(pos_) << " by " << val_ << " would give a layer of "
				<< newlayer_id << ", which is too much of a jump";
			std::string str = strs.str();
			throw std::domain_error(str);
		}

		#endif

		m_grid_delta.track(pos_, val_, layer_idx(0));
	}

	/**
	 * Create a single singularity seed point in the isogrid grid.
	 *
	 * NOTE: does not handle overwriting of points currently already on the
	 * surface/in the volume.
	 *
	 * @param pos_centre
	 */
	void seed (const VecDi& pos_centre)
	{
		// Width of seed.
		const VecDi& vec_width = VecDi::Constant(INT(LAYER_MAX));

		// Min and max positions affected by placing seed point.
		const VecDi& pos_min = pos_centre - vec_width;
		const VecDi& pos_max = pos_centre + vec_width;

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
			const VecDi& vec_dist = pos - pos_centre;
			// Sum of absolute distance along each axis == city-block distance.
			const FLOAT dist = FLOAT(vec_dist.template lpNorm<1>());
			// Check distance indicates that this point is within the narrow band.
			if (inside_band(layer_id(dist)))
			{
				// Append point to a narrow band layer (if applicable).
				layer_track(pos, dist);
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
	void update(std::function<FLOAT(const VecDi&, const IsoGrid&)> fn_)
	{
		const PosArray& pos_children = parts();
		update_start();
		FELT_PARALLEL_FOR(pos_children.size(),)
		for (UINT part_idx = 0; part_idx < pos_children.size(); part_idx++)
		{
			const VecDi& pos_part = pos_children[part_idx];
			for (const VecDi& pos : this->layer(pos_part))
				this->delta(pos, fn_(pos, m_grid_isogrid));
		}
		update_end();
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
	void update(
		const VecDi& pos_leaf_lower_, const VecDi& pos_leaf_upper_,
		std::function<FLOAT(const VecDi&, const IsoGrid&)> fn_
	) {
		static const VecDi& one = VecDi::Constant(1);
		static const VecDi& two = VecDi::Constant(2);
		// Upper and lower bounds of the grid, inclusive.
		const VecDi& pos_grid_lower = m_grid_isogrid.offset();
		const VecDi& pos_grid_upper =
			m_grid_isogrid.offset() + m_grid_isogrid.size().template cast<INT>();
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
		const VecDu& child_bounding_box_size =
			(pos_child_upper - pos_child_lower + one).template cast<UINT>();
		// Upper bound of leaf (1 more than upper point), bounded by grid..
		const VecDi& pos_leaf_upper_bound = pos_grid_upper.cwiseMin(pos_leaf_upper_ + one);
		// Upper index of bounding box.
		const UINT child_idx_bound = child_bounding_box_size.prod();
		// Clear previous update.
		update_start();
		// Parallel loop through spatial partitions.
		FELT_PARALLEL_FOR(child_idx_bound,)
		for (UINT child_idx = 0; child_idx < child_idx_bound; child_idx++)
		{
			// Get spatial partition position.
			const VecDi& pos_child_without_offset = Felt::index(
				child_idx, child_bounding_box_size
			);
			const VecDi& pos_child = pos_child_without_offset + pos_child_lower;
			// Loop all zero-layer points within this partition.
			for (const VecDi& pos : layer(pos_child))
			{
				// Skip zero-layer points not within finer-grained bounding box.
				if (IsoGrid::inside(pos, pos_leaf_lower_, pos_leaf_upper_bound))
				{
					const FLOAT amt = fn_(pos, m_grid_isogrid);

					#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

					if (m_grid_delta.get(pos) != 0)
					{
						std::stringstream strs;
						strs << "Delta is not zero: " << Felt::format(pos) << " with delta " <<
							m_grid_delta.get(pos);
						std::string str = strs.str();
						throw std::domain_error(str);
					}
					if (std::abs(amt) > 1.0f)
					{
						std::stringstream strs;
						strs << "Zero layer update value out of bounds: " << Felt::format(pos)
						<< " with value " << amt;
						std::string str = strs.str();
						throw std::domain_error(str);
					}

					#endif

					// Update delta isogrid.
					delta(pos, amt);
				}
			}
		}
		// Apply delta to isogrid.
		update_end_local();
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
					m_grid_isogrid.pos_child(pos_origin_.template cast<INT>())
				)
			);
			if (pos_hit != NULL_POS<FLOAT>())
				return pos_hit;
		}

		// Ray to test against.
		Line line(pos_origin_, dir_);

		// Tracking list for child grids that are hit.
		ChildHits child_hits;

		// Cycle each axis, casting ray to child grid planes marching away from origin.
		for (UINT dim = 0; dim < dir_.size(); dim++)
		{
			// Direction +/-1 along this axis.
			FLOAT dir_dim = sgn(dir_(dim));
			if (dir_dim == 0)
				continue;

			// Get next child plane along this axis.
			FLOAT pos_plane_dim = round_to_next(
				dim, dir_dim, pos_origin_(dim), m_grid_isogrid.child_size()
			);

			// Construct vector with elements not on this axis at zero.
			VecDf pos_plane = VecDf::Constant(0);
			pos_plane(dim) = pos_plane_dim;

			// If the zero point on this plane is not within the grid, then jump to max/min point
			// on isogrid grid.
			if (!m_grid_isogrid.inside(pos_plane))
			{
				FLOAT pos_grid_dim;
				// If casting in -'ve direction, get maximum extent.
				if (dir_dim == -1)
				{
					pos_grid_dim = m_grid_isogrid.offset()(dim) + m_grid_isogrid.size()(dim);
					if (pos_plane_dim < pos_grid_dim)
						continue;
				}
				// Else if casting in +'ve direction, get minimum extent.
				else
				{
					pos_grid_dim = m_grid_isogrid.offset()(dim);
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
			!ray_check_track_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim));

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
			const FLOAT child_size_dim = m_grid_isogrid.child_size()(dim);
			while (true)
			{
				pos_plane(dim) += dir_dim * child_size_dim;
				if (!ray_check_track_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
					break;
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
				return a.pos_child == b.pos_child;
			}
		), child_hits.end());

		// For each candidate child, cast ray through until the zero-curve is hit.
		for (const ChildHit& child_hit : child_hits)
		{
			const VecDf& pos_hit = this->ray(
				child_hit.pos_intersect, dir_,
				m_grid_isogrid.children().get(child_hit.pos_child)
			);

			if (pos_hit != NULL_POS<FLOAT>())
				return pos_hit;
		}

		return NULL_POS<FLOAT>();
	}

	/**
	 * Get null position vector for given template typename.
	 *
	 * TODO: c++14 should support variable templates, which is a better solution,
	 * but gcc doesn't yet.
	 *
	 * @return D-dimensional vector with each element set to numeric_limits<T>::max.
	 */
	template <typename T>
	static constexpr Felt::VecDT<T, D> NULL_POS()
	{
		return Felt::VecDT<T, D>::Constant(std::numeric_limits<T>::max());
	}

	/**
	 * Reset delta isogrid to zero and clear update lists.
	 *
	 * @snippet test_Surface.cpp Simple delta isogrid update
	 */
	void update_start ()
	{
		m_grid_delta.reset_all(m_grid_isogrid);
		m_grid_affected.reset_all(m_grid_isogrid);
		m_grid_affected_buffer.reset_all(m_grid_isogrid);
		m_grid_status_change.reset_all(m_grid_isogrid);
	}

	/**
	 * Update zero layer then update distance transform for all
	 * points in all layers.
	 */
	void update_end ()
	{
		update_zero_layer(&m_grid_affected_buffer);

		track_children_delta(m_grid_isogrid);

		if (update_distance(&m_grid_isogrid, &m_grid_affected_buffer))
			converge_distance(&m_grid_affected_buffer, &m_grid_affected);

		flush_status_change();

		expand_narrow_band();
	}

	/**
	 * Update zero layer then update distance transform for affected points in each layer.
	 */
	void update_end_local ()
	{
		// Get points in outer layers that are affected by changes in zero-layer.
		calc_affected();

		// Update the zero layer, applying delta to isogrid.
		update_zero_layer(&m_grid_affected_buffer);

		track_children_delta(m_grid_affected);

		converge_distance(&m_grid_affected, &m_grid_affected_buffer);

		flush_status_change();

		expand_narrow_band();
	}

#ifndef _TESTING
private:
#endif

	/**
	 * Find all outer layer points who's distance transform is affected by modified zero-layer
	 * points.
	 *
	 * @snippet test_Surface.cpp Calculate affected outer layers for localised narrow band updates
	 *
	 * TODO: several options for optimisation of removing duplicates:
	 * - Use a boolean flag grid to construct a de-duped vector (used here).
	 * - Check std::vector in Grid::neighs() using std::find to prevent tracking a duplicate in the
	 *   first place.
	 * - Use std::vector sort, unique, erase.
	 * - Use a std::unordered_set with a suitable hashing function.
	 *
	 * @param apos
	 */
	void calc_affected()
	{
		const UINT layer_idx_zero = layer_idx(0);

		// Loop over delta isogrid modified zero-layer points tracking to
		// tracking grid.

		// Loop spatial partitions of delta for zero-layer.
		for (
			const VecDi& pos_child
			: m_grid_delta.children().list(layer_idx_zero))
		{
			// Loop leaf grid nodes with spatial partition
			for (
				const VecDi& pos_leaf
				: m_grid_delta.children().get(pos_child).list(layer_idx_zero)
			)
				// track zero-layer point to tracking grid.
				m_grid_affected.track(pos_leaf, layer_idx_zero);
		}

		// Arrays to store first and last element in tracking list within each spatial partition of
		// tracking grid.
		std::array<std::vector<UINT>, NUM_LAYERS> aidx_first_neigh;
		std::array<std::vector<UINT>, NUM_LAYERS> aidx_last_neigh;

		// Loop round L times, searching outward for affected outer layer grid nodes.
		for (UINT udist = 1; udist <= LAYER_MAX; udist++)
		{
			// Reset the first and last element indices for each spatial partition in each layer.

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT layer_idx = this->layer_idx(layer_id);
				// Get number of spatial partitions for this layer.
				const UINT num_childs = (
					m_grid_affected.children().list(layer_idx).size()
				);
				// Resize spatial partition index lists for this layer to
				// to include any newly tracked partitions.
				aidx_last_neigh[layer_idx].resize(num_childs);
				// Will initialise to zero, so no further work needed for
				// these new indices giving the start of the range.
				aidx_first_neigh[layer_idx].resize(num_childs);
				// The final index needs to be copied from the current size
				// of each spatial partition, so loop over partitions,
				// copying their size into the respective last index list.
				for (
					UINT idx_child = 0; idx_child < num_childs; idx_child++
				) {
					// Get position of this spatial partition in parent
					// lookup grid.
					const PosIdx pos_child_idx = m_grid_affected.children().list(layer_idx)[idx_child];
					// Copy number of active grid nodes for this partition
					// into relevant index in the list.
					aidx_last_neigh[layer_idx][idx_child] = (
						m_grid_affected.children().get(pos_child_idx).list(layer_idx).size()
					);
				}
			}

			// Loop each layer finding the affected outer layer points
			// for each partition using the start and end points cached
			// above.

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT layer_idx = this->layer_idx(layer_id);

				// Loop over spatial partitions, ignoring newly tracked ones
				// since we're using the cached spatial partition list as
				// the end of the range.
				for (
					UINT idx_child = 0;
					idx_child < aidx_first_neigh[layer_idx].size();
					idx_child++
				) {
					// Get position of this spatial partition in this layer (not by-reference since
					// list will be modified).
					const VecDi pos_child = (
						m_grid_affected.children().list(layer_idx)[idx_child]
					);

					// Loop over leaf grid nodes within this spatial
					// partition, using the cached start and end indices,
					// so that newly tracked points are skipped.
					for (
						UINT idx_neigh = aidx_first_neigh[layer_idx][idx_child];
						idx_neigh < aidx_last_neigh[layer_idx][idx_child];
						idx_neigh++
					) {
						// Get list of active leaf grid nodes in this
						// spatial partition.
						const PosArray& apos_neigh = (
							m_grid_affected.children().get(pos_child).list(layer_idx)
						);
						// This leaf grid node is the centre to search about.
						const VecDi& pos_centre = apos_neigh[idx_neigh];

						// Use utility method from Grid to get neighbouring
						// grid nodes and call a lambda to track the point
						// to the appropriate tracking list.

						m_grid_isogrid.neighs(
							pos_centre,
							[this](const VecDi& pos_neigh) {
								// Calculate layer of this neighbouring
								// point from the isogrid grid.
								const INT layer_id_neigh = this->layer_id(pos_neigh);
								// If the calculated layer lies within the
								// narrow band, then we want to track it.
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

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT layer_idx = this->layer_idx(layer_id);
				// Loop over spatial partitions.
				for (UINT idx = 0; idx < aidx_first_neigh[layer_idx].size(); idx++)
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
		using DeltaChild = typename DeltaIsoGrid::ChildType;
		using IsoChild = typename IsoGrid::ChildType;

		const UINT layer_idx_zero = layer_idx(0);
		const PosArray& apos_children =  m_grid_delta.children().list(layer_idx_zero);

		// Performance optimisation to bulk track spatial partitions to be tracked in isogrid.
		m_grid_isogrid.track_children(apos_children, layer_idx_zero);

		const UINT num_childs = apos_children.size();

		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_buffer_))
		for (UINT idx_child = 0; idx_child < num_childs; idx_child++)
		{
			const VecDi& pos_child = apos_children[idx_child];
			DeltaChild& delta_child = m_grid_delta.children().get(pos_child);
			IsoChild& isogrid_child = m_grid_isogrid.children().get(pos_child);

			for (const PosIdx pos_idx : delta_child.list(layer_idx_zero))
			{
				const FLOAT fisogrid = isogrid_child.get(pos_idx);
				const FLOAT fdelta = delta_child.get(pos_idx);
				const FLOAT fval = fisogrid + fdelta;
				const INT layer_id_new = layer_id(fval);

				#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

				const INT layer_id_old = layer_id(fisogrid);

				if (layer_id_old != 0)
				{
					std::stringstream strs;
					strs << "Zero layer update attempted at non-zero layer point " <<
						Felt::format(pos_idx) << ": " << fisogrid << " + " << fdelta << " = " << fval;
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
						Felt::format(pos_idx) << " to " << fval << " would give a layer of " <<
						layer_id_new << ", which is too much of a jump";
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				// Update value in grid with new signed distance.
				isogrid_child(pos_idx) = fval;
				// Potentially track to status change, if narrow band layer has changed.
				this->status_change(pos_idx, 0, layer_id_new, plookup_buffer_);
			}
		}
	}

	/**
	 * Performance optimisation to bulk track spatial partitions to be tracked in delta grid.
	 *
	 * This will activate the child if it is deactivated (i.e. allocate memory for lazy grid).
	 */
	template <class GridType>
	void track_children_delta(const GridType& grid)
	{
		for (INT layer_id_child = LAYER_MIN; layer_id_child <= LAYER_MAX; layer_id_child++)
		{
			// Assume all of zero layer already tracked.
			if (layer_id_child == 0)
				continue;

			const UINT layer_idx_child = layer_idx(layer_id_child);
			m_grid_delta.track_children(grid.children().list(layer_idx_child), layer_idx_child);
		}
	}

	/**
	 * Repeatedly update the distance in the affected grid, until no more status changes are made.
	 *
	 * @param plookup initial lookup grid supplying points to check.
	 * @param plookup_buffer double-buffer grid for async iterative convergence.
	 */
	void converge_distance(AffectedLookupGrid* plookup, AffectedLookupGrid* plookup_buffer)
	{
		while (update_distance(plookup, plookup_buffer))
		{
			plookup->reset_all(m_grid_isogrid);
			std::swap(plookup, plookup_buffer);
		}
	}

	/**
	 * Perform distance transform on narrow band layers, from centre working outwards.
	 *
	 * @param plookup_ either isogrid for full update, or affected grid for local updates.
	 * @param plookup_buffer_ buffer to store grid points that need a layer update and may need
	 * trackitional iterations to converge.
	 */
	template <typename GridType>
	bool update_distance(
		const GridType* plookup_, AffectedLookupGrid* plookup_buffer_
	) {
		bool is_status_changed = false;

		// Update distance transform for inner layers of the narrow band.
		for (INT layer_id = -1; layer_id >= LAYER_MIN; layer_id--)
			is_status_changed |= update_distance(layer_id, -1, plookup_, plookup_buffer_);

		// Update distance transform for outer layers of the narrow band.
		for (INT layer_id = 1; layer_id <= LAYER_MAX; layer_id++)
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
	template <typename GridType>
	bool update_distance(
		const INT layer_id_, const INT side_,
		const GridType* plookup_, AffectedLookupGrid* plookup_buffer_
	) {
		using IsoChild = typename IsoGrid::ChildType;
		using DeltaIsoChild = typename DeltaIsoGrid::ChildType;

		bool is_status_changed = false;

		const UINT layer_idx = this->layer_idx(layer_id_);
		const PosArray& apos_children = plookup_->children().list(layer_idx);
		const UINT num_childs = apos_children.size();

		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_, plookup_buffer_))
		for (UINT pos_idx = 0; pos_idx < num_childs; pos_idx++)
		{
			// Child spatial partition position
			const VecDi& pos_child = apos_children[pos_idx];

			DeltaIsoChild& grid_delta_child = m_grid_delta.children().get(pos_child);

			const PosArray& apos_leafs = plookup_->children().get(pos_child).list(layer_idx);

			// Calculate distance of every point in this layer to the zero
			// layer, and store in delta isogrid grid.
			// Delta isogrid grid is used to allow for asynchronous updates, that
			// is, to prevent neighbouring points affecting the distance
			// transform.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance from this position to zero layer.
				const FLOAT dist = distance(pos, side_);

				#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

				const INT layer_id_new = layer_id(dist);

				if (
					layer_id_new != layer_id_ &&
					layer_id_new != layer_id_ + 1 &&
					layer_id_new != layer_id_ - 1
				) {
					const VecDi& pos_neigh = next_closest(pos, side_);

					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds." << std::endl <<
						str_neighs(pos) << "Chose " << Felt::format(pos_neigh) <<
						", giving distance of " << dist << ", which is too much of a jump";

					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				// Update delta isogrid grid.
				grid_delta_child.track(pos, dist, layer_idx);
			}
		}

		FELT_PARALLEL_FOR(num_childs, firstprivate(plookup_, plookup_buffer_))
		for (UINT pos_idx = 0; pos_idx < num_childs; pos_idx++)
		{
			const VecDi& pos_child = apos_children[pos_idx];

			IsoChild& grid_isogrid_child = m_grid_isogrid.children().get(pos_child);
			DeltaIsoChild& grid_delta_child = m_grid_delta.children().get(pos_child);

			const PosArray& apos_leafs = plookup_->children().get(pos_child).list(layer_idx);

			// Update distance in isogrid from delta isogrid and append any points that
			// move out of their layer to a status change list.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance calculated above.
				const FLOAT dist = grid_delta_child(pos);
				const INT layer_id_new = layer_id(dist);

				#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

				if (
					layer_id_new != layer_id_ &&
					layer_id_new != layer_id_ + 1 &&
					layer_id_new != layer_id_ - 1
				) {
					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds. Attempting to" <<
						" move " << Felt::format(pos) << " in layer " << layer_id_ << " to a" <<
						" distance of " << dist << " would result in a layer of " <<
						layer_id_new << ", which is too much of a jump" ;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				#endif

				grid_isogrid_child(pos) = dist;
				is_status_changed |= status_change(pos, layer_id_, layer_id_new, plookup_buffer_);
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
	 * @param plookup_buffer_ grid to store points that need status change, for potential trackitional
	 * processing.
	 * @return true if status change is needed for this position, false otherwise.
	 */
	bool status_change (
		const VecDi& pos_, const INT layer_id_from_, const INT layer_id_to_,
		AffectedLookupGrid* plookup_buffer_
	) {
		if (layer_id_from_ == layer_id_to_)
			return false;

		#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)
		m_grid_isogrid.assert_pos_bounds(pos_, "status_change: ");
		#endif

		INT& layer_id_to = m_grid_status_change.get(pos_);

		if (layer_id_to != m_grid_status_change.background())
		{
			layer_id_to = layer_id_to_;
		}
		else
		{
			m_grid_status_change.track(
				pos_, layer_id_to_, layer_idx(layer_id_from_)
			);
		}

		if (inside_band(layer_id_to_))
			plookup_buffer_->track(pos_, layer_idx(layer_id_to_));

		return true;
	}

	/**
	 * Loop through the status change lists moving the points from one layer to another.
	 */
	void flush_status_change ()
	{
		using StatusChangeChild = typename StatusChangeGrid::ChildType;
		using StatusChangeChildrenList = typename StatusChangeGrid::Lookup::PosArray;

		for (INT layer_id_from = LAYER_MIN; layer_id_from <= LAYER_MAX; layer_id_from++)
		{
			const UINT list_idx = this->layer_idx(layer_id_from);

			const StatusChangeChildrenList& pos_children =
				m_grid_status_change.children().list(list_idx);

			const UINT num_childs = pos_children.size();

			FELT_PARALLEL_FOR(num_childs, firstprivate(layer_id_from))
			for (UINT idx_child = 0; idx_child < num_childs; idx_child++)
			{
				const VecDi& pos_child = pos_children[idx_child];
				const StatusChangeChild& child = m_grid_status_change.children().get(pos_child);

				for (const PosIdx pos_idx : child.list(list_idx))
				{
					const INT layer_id_to = child.get(pos_idx);
					this->layer_move(pos_idx, layer_id_from, layer_id_to);
				}
			}
		}
	}

	/**
	 * track new points to the narrow band when expanding/contracting.
	 */
	void expand_narrow_band()
	{
		using StatusChangeChild = typename StatusChangeGrid::ChildType;
		using StatusChangeChildrenList = typename StatusChangeGrid::Lookup::PosArray;

		// Cycle innermost layer and outermost layer.
		for (
			INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id += NUM_LAYERS-1
		) {
			const UINT list_idx = this->layer_idx(layer_id);

			const StatusChangeChildrenList& apos_children =
				m_grid_status_change.children().list(list_idx);

			const INT side = sgn(layer_id);

			// TODO: it appears this is (no longer) thread-safe...
//			#pragma omp parallel for
			for (UINT idx_child = 0; idx_child < apos_children.size(); idx_child++)
			{
				const PosIdx pos_child_idx = apos_children[idx_child];
				const StatusChangeChild& child = m_grid_status_change.children().get(pos_child_idx);

				for (const PosIdx pos_idx : child.list(list_idx))
				{
					// If not expanding/contracting, then nothing to do here.
					if (child.get(pos_idx) != LAYER_MAX * side - side)
						continue;

					// Cycle over neighbours of this outer layer point.
					m_grid_isogrid.neighs(
						pos_idx,
						[this, list_idx, side, layer_id, pos_idx](const VecDi& pos_neigh) {
							const INT layer_id_from = this->layer_id(pos_neigh);

							// Only track if neighbouring point is not already within the narrow band.
							if (inside_band(layer_id_from))
							{
								#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

								const UINT lookup_idx = m_grid_isogrid.children().get(
										m_grid_isogrid.pos_child(pos_neigh)
									).lookup().get(pos_neigh);

								if (lookup_idx == IsoGrid::ChildType::Lookup::NULL_IDX)
								{
									std::stringstream sstr;
									sstr << "pos_idx not tracked but should be: " << str_pos(pos_neigh);
									std::string str = sstr.str();
									throw std::domain_error(str);
								}

								#endif

								return;
							}

							// Calculate distance of this neighbour to the zero curve.
							const FLOAT distance_neigh = distance(pos_neigh, side);

							#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

							const INT layer_id_to = this->layer_id(distance_neigh);

							if (layer_id_to != layer_id)
							{
								std::stringstream strs;
								strs << "Neighbour is further away than expected." << std::endl <<
									"pos_idx:" << std::endl <<
									"  " << str_pos(pos_idx) << std::endl
									<< "Neigh:" << std::endl <<
									"  " << str_pos(pos_neigh) << std::endl <<
									"Calculated distance " << distance_neigh <<
									" would give a layer of " << layer_id_to <<
									" when we expect a layer of " << layer_id;
								std::string str = strs.str();
								throw std::domain_error(str);
							}

							if (layer_id_to != LAYER_MIN && layer_id_to != LAYER_MAX)
							{
								std::stringstream strs;
								strs << "Attempting to track " << Felt::format(pos_neigh) <<
									" to the narrow band but the distance is " << distance_neigh <<
									" which would give a layer of " << layer_id_to;
								std::string str = strs.str();
								throw std::domain_error(str);
							}

							#endif

							// Use thread-safe update and track function, since neighbouring point
							// could be in another spatial partition.
							this->m_grid_isogrid.track(pos_neigh, distance_neigh, list_idx);
						}
					);
				} // End for pos_idx.
			}
		}
	}

	/**
	 * Append position to the narrow band using given value, if value is within the band.
	 *
	 * @param pos_ position to potentially track to narrow band.
	 * @param val_ value to set position to.
	 */
	void layer_track (const VecDi& pos_ ,const FLOAT val_)
	{
		const INT layer_id_pos = layer_id(val_);
		// Do nothing if position is outside narrow band.
		if (!inside_band(layer_id_pos))
			return;
		m_grid_isogrid.track(val_, pos_, layer_idx(layer_id_pos));
	}

	/**
	 * untrack a point from one layer and move it into another.
	 *
	 * Will only modify underlying grid value to if removing from a layer without tracking
	 * to another (i.e. pos_ is now outside the narrow band) and will set to a background value
	 * of layer_id_to_).
	 *
	 * @param pos_ position to move between layers.
	 * @param layer_id_from_ layer moving from.
	 * @param layer_id_to_ layer moving to.
	 */
	void layer_move(
		const VecDi& pos_, const INT layer_id_from_, const INT layer_id_to_
	) {
		const bool is_from_inside = inside_band(layer_id_from_);
		const bool is_to_inside = inside_band(layer_id_to_);

		if (is_from_inside && is_to_inside)
		{

			#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

			const VecDi& pos_child = m_grid_isogrid.pos_child(pos_);
			const TupleIdx list_idx_from = layer_idx(layer_id_from_);
			const TupleIdx list_idx_to = layer_idx(layer_id_to_);
			const typename IsoGrid::ChildType& child = m_grid_isogrid.children().get(pos_child);

			if (child.list(list_idx_from).size() == 0)
			{
				std::stringstream strs;
				strs << "Layer empty when attempting to move " << Felt::format(pos_) <<
					" from layer " << layer_id_from_ << " to layer " << layer_id_to_ <<
					" in partition " << Felt::format(pos_child) << " = " <<
					Felt::format(child.offset()) << "-" <<
					Felt::format(child.offset() + child.size().template cast<INT>());
				std::string str = strs.str();
				throw std::domain_error(str);
			}

			#endif

			m_grid_isogrid.retrack(
				pos_, layer_idx(layer_id_from_), layer_idx(layer_id_to_)
			);
		}
		else if (is_from_inside)
		{
			m_grid_isogrid.untrack(
				pos_, layer_idx(layer_id_from_), layer_id_to_
			);
		}
		else if (is_to_inside)
		{
			m_grid_isogrid.track(pos_, layer_idx(layer_id_to_));
		}
		#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)
		else
		{
			std::stringstream strs;
			strs << "Invalid layer move: attempting to move " << Felt::format(pos_) <<
				" from layer " << layer_id_from_ << " to layer " << layer_id_to_;
			std::string str = strs.str();
			throw std::domain_error(str);
		}
		#endif
	}

	/**
	 * Calculate city-block distance from position to zero curve.
	 *
	 * @param pos_ target position.
	 * @param side_ side of narrow band (+/-1).
	 * @return city-block distance from pos_ to the zero curve.
	 */
	FLOAT distance (const VecDi& pos_, const FLOAT side_) const
	{
		// Get neighbouring point that is next closest to the zero-layer.
		const VecDi pos_closest = next_closest(pos_, side_);
		const FLOAT val_closest = m_grid_isogrid(pos_closest);
		// This point's distance is then the distance of the closest neighbour +/-1, depending
		// which side of the band we are looking at.
		const FLOAT dist = val_closest + side_;


		#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

		const INT layer_id_pos = layer_id(pos_);
		const INT layer_id_neigh = layer_id(val_closest);

		if (
			std::abs(layer_id_pos) < std::abs(layer_id_neigh) &&
			sgn(layer_id_pos) == sgn(layer_id_neigh)
		) {
			std::stringstream sstr;
			sstr << "Neighbour closest to zero curve is further away than source position: " <<
				Felt::format(pos_) << " at " << m_grid_isogrid.get(pos_) << " is closer than " <<
				Felt::format(pos_closest) << " at " << m_grid_isogrid.get(pos_closest) <<
				" but should not be";
			std::string str = sstr.str();
			throw std::domain_error(str);
		}

		#endif

		return dist;
	}

	/**
	 * Get neighbouring position in isogrid grid that is closest to zero-curve.
	 *
	 * @param pos_ position in grid to get neighbour for.
	 * @param side_ side of narrow band (+/-1).
	 * @return neighbour of pos_ closest to the zero-curve, or pos_ itself if already on zero curve.
	 */
	VecDi next_closest (const VecDi& pos_, const FLOAT side_) const
	{
		// Trivially return if this is already a zero-layer point.
		if (this->layer_id(pos_) == 0)
			return pos_;

		VecDi pos_nearest(pos_);
		FLOAT val_nearest = m_grid_isogrid(pos_)*side_;

		m_grid_isogrid.neighs(
			pos_,
			[this, &pos_nearest, &val_nearest, side_](const VecDi& pos_neigh_) {
				const FLOAT val_neigh = m_grid_isogrid.get(pos_neigh_);
				// Check absolute value of this neighbour is less than nearest point.
				// Multiplying by side_ has two effects:
				// - It has the same effect as abs() for points on the same side of the band.
				// - It ensures points on the opposite side are negative, so that < comparison
				//   prefers those points, which is good because we're interested in the neighbour
				//   in the *direction* of the zero-curve.
				if (val_neigh*side_ < val_nearest) {
					pos_nearest = pos_neigh_;
					val_nearest = val_neigh*side_;
				}
			}
		);

		return pos_nearest;
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
		VecDf pos_sample, const VecDf& dir, const typename IsoGrid::ChildType& child
	) const {
		using Line = Eigen::ParametrizedLine<FLOAT, D>;

		const Line line_leaf(pos_sample, dir);
		FLOAT t_leaf = 0;

//		std::cerr << "Child: (" << child.offset().transpose() << ") + (" <<
//			child.size().transpose() << ")" << std::endl;

		while (child.inside(pos_sample))
		{
			const INT layer_id = this->layer_id(pos_sample);

//			std::cerr << layer_id << std::endl;

			if (abs(layer_id) == 0)
			{
				VecDf normal = m_grid_isogrid.grad(pos_sample);

				#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)
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
					static const UINT MAX_CONVERGE_STEPS = 100;
					UINT num_converge_steps = 0;
					FLOAT dist = 0;
					for (; num_converge_steps < MAX_CONVERGE_STEPS; num_converge_steps++)
					{
						dist = m_grid_isogrid.interp(pos_sample);

						pos_sample -= normal*dist;

						if (!m_grid_isogrid.inside(pos_sample))
							return NULL_POS<FLOAT>();

						if (std::abs(dist) <= TINY || normal.dot(dir) >= 0)
							break;

						normal = m_grid_isogrid.grad(pos_sample);
						normal.normalize();
					}

					#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)
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

		return NULL_POS<FLOAT>();
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
		const VecDf& pos_intersect = (
			line.intersectionPoint(plane) + line.direction() * FLOAT(TINY)
		);

		const VecDu& size = m_grid_isogrid.size();
		const VecDi& offset = m_grid_isogrid.offset();
		const VecDf& dir = line.direction();

		for (UINT i = 0; i < dir.size(); i++)
		{
			if (
				(dir(i) > 0 && pos_intersect(i) > size(i)) ||
				(dir(i) < 0 && pos_intersect(i) < offset(i))
			)
				return false;
		}

		if (!m_grid_isogrid.inside(pos_intersect))
			return true;

		const VecDi& pos_floor = floor(pos_intersect);
		const VecDi& pos_child = this->m_grid_isogrid.pos_child(pos_floor);

		if (
			this->layer(pos_child, 0).size() ||
			this->layer(pos_child, 1).size() || this->layer(pos_child, -1).size()
		) {
			child_hits.push_back(ChildHit(pos_intersect, pos_child));
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
	FLOAT round_to_next(
		const UINT dim, const FLOAT dir, const FLOAT pos, const VecDu part_size
	) const {
		// Real-valued child pos translated to [0, 2*childsize) space.
		FLOAT pos_plane_dim = (
			FLOAT(pos - m_grid_isogrid.offset()(dim)) / part_size(dim)
		);
		// Round to next child en route in [0, 2*childsize) space.
		pos_plane_dim = (dir == -1) ?
			std::floor(pos_plane_dim) : std::ceil(pos_plane_dim);
		// Scale back to isogrid grid in [0, 2*fullsize) space.
		pos_plane_dim *= m_grid_isogrid.child_size()(dim);
		// Translate back to isogrid grid in [-fullsize, fullsize) space.
		pos_plane_dim += m_grid_isogrid.offset()(dim);

		return pos_plane_dim;
	}


#if defined(Felt_EXCEPTIONS) || !defined(NDEBUG)

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
		const FLOAT dist_pos = m_grid_isogrid.get(pos_);
		const FLOAT layer_id_pos = layer_id(pos_);
		const VecDi& pos_child = m_grid_isogrid.pos_child(pos_);
		const typename IsoGrid::ChildType child = m_grid_isogrid.children().get(pos_child);
		const VecDi& pos_child_lower = child.offset();
		const VecDi& pos_child_upper =
			child.offset() + child.size().template cast<INT>();
		const Felt::VecDu<NUM_LAYERS>& list_idxs_child =
			m_grid_isogrid.children().lookup().get(pos_child);
		const UINT list_id_pos = layer_idx(layer_id_pos);
		const TupleIdx list_idx_pos = m_grid_isogrid.children().get(pos_child).lookup().get(pos_);

		std::stringstream sstr;
		sstr << Felt::format(pos_) << " ∈ P(" << Felt::format(pos_child) << ") = [" <<
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
	template <class PosType>
	LayerId layer_id(const PosType& pos_) const
	{
		const LayerId layer_id_pos = layer_id(m_grid_isogrid(pos_));

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
	LayerId layer_id(const FLOAT val_) const
	{
		// Round to value+epsilon, to catch cases of precisely +/-0.5.
		return boost::math::iround(
			val_ + std::numeric_limits<FLOAT>::epsilon()
		);
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
		return TupleIdx(id_ + LayerId(NUM_LAYERS) / 2);
	}

	/**
	 * Test whether a given value lies (or should lie) within the narrow band or not.
	 *
	 * @param val float or int value to test
	 * @return true if inside the narrow band, false otherwise.
	 */
	template <typename ValType>
	bool inside_band (const ValType& val) const
	{
		return (UINT)std::abs(val) <= LAYER_MAX;
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
};
}
#endif