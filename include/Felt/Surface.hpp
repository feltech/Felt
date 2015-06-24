#ifndef Surface_hpp
#define Surface_hpp

#include <vector>
#include <functional>
#include <limits>
#include <boost/math/special_functions/round.hpp>
#include <eigen3/Eigen/Dense>
#include <omp.h>

#include "PartitionedGrid.hpp"


namespace felt {

	/**
	 * A n-dimensional sparse-field spatially partitioned level set.
	 *
	 * @tparam D the number of dimensions of the surface.
	 * @tparam L the number of narrow band layers surrounding the zero-level
	 * surface.
	 * @tparam P the size of a spatial partition.
	 */
	template <UINT D, UINT L=2>
	class Surface
	{
	public:
		/**
		 * A delta phi update grid with active (non-zero) grid points tracked.
		 */
		typedef TrackedPartitionedGrid<FLOAT, D, 2*L+1>		DeltaPhiGrid;
		/**
		 * A level set embedding phi grid, with active grid points (the narrow
		 * band) tracked.
		 */
		typedef SharedTrackedPartitionedGrid<FLOAT, D, 2*L+1>	PhiGrid;

		static const INT LAYER_MIN	= -L;
		static const INT LAYER_MAX	= L;

	protected:

		/**
		 * A resizable array of D-dimensional grid positions.
		 */
		typedef typename PhiGrid::PosArray	PosArray;
		/**
		 * D-dimensional unsigned int vector.
		 */
		typedef typename PhiGrid::VecDu	VecDu;
		/**
		 * D-dimensional integer vector.
		 */
		typedef typename PhiGrid::VecDi VecDi;
		/**
		 * D-dimensional float vector.
		 */
		typedef typename PhiGrid::VecDf VecDf;

		/**
		 * Storage class for flagging a point in the grid to be moved from one
		 * narrow band layer to another.
		 */
		struct StatusChange
		{
			StatusChange(
				const VecDi& pos_, const INT& from_layer_, const INT& to_layer_
			) : pos(pos_), from_layer(from_layer_), to_layer(to_layer_){}

			VecDi pos;
			INT from_layer;
			INT to_layer;
		};

		/**
		 * A spatially partitioned array of StatusChange objects.
		 */
		typedef PartitionedArray<StatusChange, D>	StatusChangeGrid;

		typedef LookupPartitionedGrid<D, 1>			AffectedLookupGrid;


	protected:

		/**
		 * The minimum usable position in the grid for the surface (zero-layer).
		 *
		 * Additional space is required at the border of the phi embedding for
		 * the outer layers, so this is the effective minimum position a point
		 * on the surface can occupy.
		 */
		VecDi m_pos_min;
		/**
		 * The maximum usable position in the grid for the surface (zero-layer).
		 *
		 * Additional space is required at the border of the phi embedding for
		 * the outer layers, so this is the effective maximum position a point
		 * on the surface can occupy.
		 */
		VecDi m_pos_max;

		/**
		 * Grid for preventing duplicates when doing neighbourhood queries.
		 */
		AffectedLookupGrid 	m_grid_affected;

		/**
		 * The main level set embedding phi grid.
		 *
		 * Named after the greek letter often used to represent the level set
		 * function.
		 */
		PhiGrid				m_grid_phi;

		/**
		 * The delta phi update grid.
		 *
		 * Used to allow for asynchronous updating.
		 */
		DeltaPhiGrid 		m_grid_dphi;

		/**
		 * The (spatially partitioned) status change list.
		 *
		 * The list is appended to when a point in the narrow band moves from
		 * one layer to another.
		 */
		StatusChangeGrid	m_grid_status_change;


	public:
//		/**
//		 * Default constructor initialising a zero-dimensional embedding.
//		 */
//		Surface ()
//		:	m_grid_phi(),
//			m_grid_status_change(),
//			m_grid_dphi()
//		{}

		/**
		 * Construct a level set embedding of size dims.
		 *
		 * All points will be marked as outside the surface (i.e. no surface).
		 *
		 * @param dims
		 * @param uborder
		 */
		Surface (
			const VecDu& dims,
			const VecDu& dims_partition = VecDu::Constant(DEFAULT_PARTITION)
		)
		:	m_grid_phi(),
			m_grid_status_change(),
			m_grid_dphi()
		{
			this->dims(dims, dims_partition);
		}

		/**
		 * Shorthand for accessing the phi (level set) grid at a given position.
		 *
		 * @param pos
		 * @return
		 */
		FLOAT& operator() (const VecDi& pos)
		{
			return this->phi()(pos);
		}

		/**
		 * Shorthand for accessing the phi (level set) grid (const version).
		 *
		 * @param pos
		 * @return
		 */
		const FLOAT& operator() (const VecDi& pos) const
		{
			return this->phi()(pos);
		}

		/**
		 * Initialise level set embedding with given dimensions.
		 *
		 * Initialises the various lookup grids and (indirectly) the level set
		 * sparse field layers. Also calculates and stores the spatial limits
		 * of the grid accounting for the narrow band space required.
		 *
		 * @param udims
		 */
		void dims (const VecDu& udims, const VecDu& dims_partition)
		{
			const VecDi idims = udims.template cast<INT>();
			const VecDi offset = -idims/2;

			m_grid_phi.init(udims, offset, dims_partition);

			// Configure delta phi embedding.
			m_grid_dphi.init(udims, offset, dims_partition);
			// Configure status change partitioned lists.
			m_grid_status_change.init(udims, offset, dims_partition);

			// Configure de-dupe grid for neighbourhood queries.
			m_grid_affected.init(udims, offset, dims_partition);

			// Store min and max usable positions in phi embedding.
			this->pos_min(
				VecDi::Constant(L + 1) + m_grid_phi.offset()
			);
			this->pos_max(
				(idims - VecDi::Constant(L + 1))
				+ m_grid_phi.offset() - VecDi::Constant(1)
			);
			// Fill phi grid with 'outside' value.
			m_grid_phi.fill(L+1);
			// Initialise delta phi to zero.
			m_grid_dphi.fill(0);
		}

		/**
		 * Get minimum usable position in phi grid.
		 *
		 * @return
		 */
		const VecDi& pos_min () const
		{
			return m_pos_min;
		}

		/**
		 * Get maximum usable position in phi grid.
		 *
		 * @return
		 */
		const VecDi& pos_max () const
		{
			return m_pos_max;
		}


		/**
		 * Set minimum usable position in phi grid.
		 *
		 * @param pos
		 */
		void pos_min (const VecDi& pos)
		{
			m_pos_min = pos;
		}

		/**
		 * Set maximum usable position in phi grid.
		 *
		 * @param pos
		 */
		void pos_max (const VecDi& pos)
		{
			m_pos_max = pos;
		}


		/**
		 * Get reference to phi grid.
		 *
		 * @return
		 */
		PhiGrid& phi ()
		{
			return m_grid_phi;
		}

		/**
		 * Get reference to phi grid.
		 *
		 * @return
		 */
		const PhiGrid& phi () const
		{
			return m_grid_phi;
		}

		/**
		 * Shorthand for accessing phi grid at given position (const
		 * version).
		 *
		 * @param pos
		 * @return
		 */
		const FLOAT& phi (const VecDi& pos) const
		{
			return m_grid_phi(pos);
		}

		/**
		 * Shorthand for accessing phi grid at given position.
		 *
		 * @param pos
		 * @return
		 */
		FLOAT& phi (const VecDi& pos)
		{
			return m_grid_phi(pos);
		}

		/**
		 * Test whether a given value lies within the narrow band or not.
		 *
		 * @param val
		 * @return
		 */
		template <typename ValType>
		bool inside_band (const ValType& val)
		{
			return (UINT)std::abs(val) <= L;
		}

		/**
		 * Update phi grid point at pos by val.
		 *
		 * Checks if the layer should change and if so adds to the appropriate
		 * status change list.
		 *
		 * Also expands the outer layers as the surface expands/contracts.
		 *
		 * TODO: because of the outer layer expansion code, this function is
		 * not, in general, thread safe.  Must move outer layer expansion to a
		 * separate routine.
		 *
		 * @param pos
		 * @param val
		 * @param layerID
		 */
		void phi (const VecDi& pos, const FLOAT& val, const INT& layerID = 0)
		{
			PhiGrid& phi = this->phi();
			const INT newLayerID = this->layerID(val);
			phi(pos) = val;

			if (newLayerID == layerID)
				return;

			this->status_change(pos, layerID, newLayerID);

			// If outside point moving inward, must create new outside
			// points.
			if (!(std::abs(layerID) == L && std::abs(newLayerID) == L-1))
				return;

			// Get neighbouring points.
			PosArray neighs;
			phi.neighs(pos, neighs);
			// Get which side of the zero-layer this point lies on.
			const INT side = sgn(newLayerID);

			for (
				UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++
			) {
				const VecDi& pos_neigh = neighs[neighIdx];
				const INT fromLayerID = this->layerID(pos_neigh);
				// If neighbouring point is not already within the
				// narrow band.
				if (!this->inside_band(fromLayerID))
				{
					// Get distance of this new point to the zero layer.
					const FLOAT dist_neigh = this->distance(
						pos_neigh, side
					);
					// Set distance in phi grid.
					phi(pos_neigh) = dist_neigh;
					// Add to status change list to be added to the
					// outer layer.
					this->status_change(pos_neigh, fromLayerID, side*L);
				}
			}
		}

		/**
		 * Add a point to the status change list to eventually be moved
		 * from one layer to another.
		 *
		 * @param pos
		 * @param fromLayerID
		 * @param toLayerID
		 */
		void status_change (
			const VecDi& pos, const INT& fromLayerID, const INT& toLayerID
		) {
			m_grid_status_change.add(
				pos, StatusChange(pos, fromLayerID, toLayerID)
			);
		}

		/**
		 * Loop through the status change lists moving the referenced points
		 * from one layer to another.
		 */
		void status_change ()
		{
			for (
				const VecDi& pos_child : m_grid_status_change.branch().list()
			) {
				for (
					const StatusChange& change
					: m_grid_status_change.child(pos_child)
				) {
					this->layer_move(
						change.pos, change.from_layer, change.to_layer
					);
				}
			}
		}

		/**
		 * Get reference to delta phi grid.
		 *
		 * @return
		 */
		DeltaPhiGrid& dphi ()
		{
			return m_grid_dphi;
		}

		/**
		 * Get reference to delta phi grid (const version).
		 *
		 * @return
		 */
		const DeltaPhiGrid& dphi () const
		{
			return m_grid_dphi;
		}

		/**
		 * Shorthand for access to the delta phi grid at given position.
		 *
		 * @param pos
		 * @return
		 */
		FLOAT& dphi (const VecDi& pos)
		{
			return m_grid_dphi(pos);
		}

		/**
		 * Shorthand for access to the delta phi grid at given position (const
		 * version)
		 *
		 * @param pos
		 * @return
		 */
		const FLOAT& dphi (const VecDi& pos) const
		{
			return m_grid_dphi(pos);
		}

		/**
		 * Update delta phi grid and append point to change list for
		 * given thread.
		 * @param pos
		 * @param val
		 */
		void dphi (
			const VecDi& pos, FLOAT val, const INT& layerID = 0
		) {
			// If this is the zero-layer, then ensure we cannot leave the grid
			// boundary.
			if (layerID == 0)
			{
				const VecDi& pos_min = this->pos_min();
				const VecDi& pos_max = this->pos_max();
				// Cycle each axis.
				for (UINT d = 0; d < D; d++)
					// Check if pos lies at the max bound of this axis.
					if (pos_min(d) == pos(d) || pos_max(d) == pos(d))
					{
						// Get phi at this point.
						const FLOAT fphi = this->phi(pos);
						// Max value that will not be rounded and thus trigger
						// a layer_move.
						const FLOAT val_max = -0.5 +
							(std::numeric_limits<FLOAT>::epsilon() * 2);
						// Clamp the value of delta phi.
						val = std::max(val_max - fphi, val);
						break;
					}
			}

			this->dphi().add(pos, val, this->layerIdx(layerID));
		}


		/**
		 * Get reference to the active partitions of a given layer of the narrow
		 * band.
		 *
		 * @param layerID
		 * @return
		 */
		PosArray& parts (const INT& id = 0)
		{
			return m_grid_phi.branch().list(id+L);
		}

		/**
		 * Get reference to a single layer of the narrow band at a given
		 * spatial partition.
		 *
		 * @param pos_child
		 * @param id
		 * @return
		 */
		PosArray& layer (const VecDi& pos_child, const INT& id = 0)
		{
			return m_grid_phi.child(pos_child).list(id+L);
		}

		/**
		 * Get reference to a single layer of the narrow band at a given
		 * spatial partition (const version).
		 *
		 * @param pos_child
		 * @param id
		 * @return
		 */
		const PosArray& layer (
			const VecDi& pos_child, const INT& id = 0
		) const
		{
			return m_grid_phi.child(pos_child).list(id+L);
		}

		/**
		 * Get a container providing an iterator to the grid positions in a
		 * given layer of the narrow band.
		 *
		 * Spatial partitioning complicates matters, so a special data structure
		 * is required for iterating over all the positions in a layer in
		 * sequence.
		 *
		 * @param layerID
		 * @return
		 */
		const LeafsContainer<PhiGrid> layer(const UINT& layerID) const
		{
			return m_grid_phi.leafs(this->layerIdx(layerID));
		}

		/**
		 * Append position to a layer of the narrow band.
		 *
		 * @param id
		 * @param pos
		 */
		void layer_add (const VecDi& pos, const INT& layerID)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;
			m_grid_phi.add(pos, this->layerIdx(layerID));
		}

		/**
		 * Append position to a layer of the narrow band.
		 *
		 * Skip the lookup in the phi grid by assuming passed val is the phi
		 * grid value of this point.
		 *
		 * @param pos
		 * @param val
		 */
		void layer_add (const FLOAT& val, const VecDi& pos)
		{
			const INT& id = this->layerID(val);
			this->layer_add(pos, id);
		}

		/**
		 * Remove a position from a given layer of the narrow band.
		 *
		 * Does not modify underlying grid value, just the layer list.
		 *
		 * @param pos
		 * @param layerID
		 */
		void layer_remove(const VecDi& pos, const INT& layerID)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;

			m_grid_phi.remove(pos, this->layerIdx(layerID));
		}

		/**
		 * Remove a point from one layer and move it into another.
		 *
		 * Simply adjusts the lists, does not modify the underlying grid values.
		 *
		 * @param pos
		 * @param fromLayerID
		 * @param toLayerID
		 */
		void layer_move(
			const VecDi& pos, const INT& fromLayerID, const INT& toLayerID
		) {
			this->layer_remove(pos, fromLayerID);
			this->layer_add(pos, toLayerID);
		}

		/**
		 * Get narrow band layer id of location in phi grid.
		 * @param pos
		 * @return
		 */
		INT layerID(const VecDi& pos) const
		{
			return this->layerID(this->phi(pos));
		}

		/**
		 * Get narrow band layer id of value.
		 * @param val
		 * @return
		 */
		INT layerID(const FLOAT& val) const
		{
			// Round to value+epsilon, to catch cases of precisely +/-0.5.
			return boost::math::round(
				val + std::numeric_limits<FLOAT>::epsilon()
			);
		}

		/**
		 * Get narrow band layer index of ID.
		 *
		 * @param id
		 * @return
		 */
		UINT layerIdx(const INT& id) const
		{
			return id+L;
		}

		/**
		 * Create a single singularity seed point in the phi grid.
		 * NOTE: does not handle overwriting of points currently already on the
		 * surface/in the volume.
		 *
		 * @param pos_centre
		 */
		void seed (const VecDi& pos_centre)
		{
			PhiGrid& phi = this->phi();
			const VecDu dims = phi.dims();

			// Width of seed.
			const VecDi vec_width = VecDi::Constant(L);

			// Min and max positions affected by placing seed point.
			const VecDi pos_min = pos_centre - vec_width;
			const VecDi pos_max = pos_centre + vec_width;

			// Get vector size of window formed by pos_min and pos_max.
			const VecDu pos_size = (
				pos_max - pos_min + VecDi::Constant(1)
			).template cast<UINT>(); //+1 for zero coord.

			// Calculate number of grid points to be cycled through within
			// window.
			UINT size = 1;
			for (INT axis = 0; axis < dims.size(); axis++)
				size *= pos_size(axis);

			// Cycle through each point in window.
			for (UINT u_pos = 0; u_pos <= size; u_pos++)
			{
				// Calculate vector position from integer index,
				// using Felt::Grid utility function, index().
				VecDi pos = PhiGrid::index(u_pos, pos_size);
				// Translate position into phi grid space.
				pos += pos_min;
				// Calculate vector distance from this position to seed centre.
				const VecDi vec_dist = pos - pos_centre;
				// Sum of absolute distance along each axis == city-block
				// distance.
				FLOAT f_dist = (FLOAT)vec_dist.template lpNorm<1>();
				// Check distance indicates that this point is within the
				// narrow band.
				if ((UINT)std::abs(this->layerID(f_dist)) <= L)
				{
					// Set distance as value in phi grid.
					phi(pos) = f_dist;
					// Append point to a narrow band layer (if applicable).
					this->layer_add(f_dist, pos);
				}
			}
		}

		/**
		 * Get neighbouring position in phi grid that is closest to
		 * zero-curve.
		 *
		 * @param pos
		 * @param side
		 * @return
		 */
		VecDi next_closest (const VecDi& pos, const FLOAT& side) const
		{
			// Trivially return if this is already a zero-layer point.
			if (this->layerID(pos) == 0)
				return pos;

			const PhiGrid& phi = this->phi();

			// Get all neighbours of this point.
			PosArray neighs;
			phi.neighs(pos, neighs);

			VecDi pos_nearest = VecDi(pos);
			FLOAT val_nearest = phi(pos)*side;
			// Cycle neighbours finding one that is closest to zero-layer.
			for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
			{
				const VecDi pos_neigh = neighs[neighIdx];
				const FLOAT val_neigh = phi(pos_neigh);
				// Check absolute value of this neighbour is less than nearest
				// point.
				// NOTE: cannot simply use abs() because during an update,
				// points close to the zero-curve may end up detecting points
				// on the other side.  By multiplying by the side (+/-1) value,
				// we ensure points on the opposite side of the band will
				// always be considered further away from the zero-layer than
				// points on the same side.
				if (val_neigh*side < val_nearest) {
					pos_nearest = pos_neigh;
					val_nearest = val_neigh*side;
				}
			}

			// TODO: lovely elegant solution using gradient vector doesn't work
			// in all cases because pos may not yet be initialised (this
			// next_closest() function is used to initialise it), so gradient
			// can be erroneous.
//			const VecDf vec_dir = phi.grad(pos) * dir;
//			const UINT axis_best = ublas::index_norm_inf(vec_dir);
//			VecDi pos_nearest = VecDi(pos);
//			const FLOAT val_best = vec_dir(axis_best);
//			if (fabs(val_best) > 0)
//				pos_nearest(axis_best) += sgn(val_best);

			return pos_nearest;
		}

		/**
		 * Get neighbouring position in phi grid that is closest to the
		 * zero-curve.
		 *
		 * @param pos
		 * @return
		 */
		VecDi next_closest (const VecDi& pos) const
		{
			const Grid<FLOAT,D>& phi = this->phi();
			const FLOAT val_centre = phi(pos);
			// Direction multiplier for gradient toward zero-curve.
			const FLOAT side = sgn(val_centre);

			return this->next_closest(pos, side);
		}

		/**
		 * Reset delta phi to zero and clear update lists.
		 */
		void update_start ()
		{
			for (UINT layerIdx = 0; layerIdx < 2*L+1; layerIdx++)
				this->dphi().reset(0, layerIdx);

			m_grid_status_change.reset();
			m_grid_affected.reset();
		}


		/**
		 * Apply delta phi to phi along the zero layer.
		 */
		void update_zero_layer ()
		{
			const typename DeltaPhiGrid::BranchGrid&
			branch = this->dphi().branch();

			const UINT& layerIdx = this->layerIdx(0);

//			#pragma omp parallel for
			for (
				UINT idx_child = 0; idx_child < branch.list(layerIdx).size();
				idx_child++
			) {
				const VecDi& pos_child = branch.list(layerIdx)[idx_child];
				for (const VecDi& pos : branch(pos_child).list(layerIdx))
				{
					const FLOAT& fphi = this->phi(pos);
					const FLOAT& fdphi = this->dphi(pos);
					const FLOAT& fval = fphi + fdphi;
					this->phi(pos, fval);
				}
			}
		}

		/**
		 * Update zero layer then update distance transform for all
		 * points in all layers.
		 */
		void update_end ()
		{
			this->update_zero_layer();

			// Update distance transform for inner layers of the narrow band.
			for (INT layerID = -1; layerID >= -(INT)L; layerID--)
				this->update_distance(layerID, -1);
			// Update distance transform for outer layers of the narrow band.
			for (INT layerID = 1; layerID <= (INT)L; layerID++)
				this->update_distance(layerID, 1);

			this->status_change();
		}

		/**
		 * Update zero layer then update distance transform for affected
		 * points in each layer.
		 */
		void update_end_local ()
		{
			// Get points in outer layers that are affected by changes in
			// zero-layer.
			PosArray aAffected[2*L+1];
			this->affected(aAffected);

			// Update the zero layer, applying delta to phi.
			this->update_zero_layer();

			// Update distance transform for inner layers of the narrow band.
			for (INT layerID = -1; layerID >= -(INT)L; layerID--)
			{
				PosArray& apos = aAffected[this->layerIdx(layerID)];
				this->update_distance(layerID, -1, apos);
			}
			// Update distance transform for outer layers of the narrow band.
			for (INT layerID = 1; layerID <= (INT)L; layerID++)
			{
				PosArray& apos = aAffected[this->layerIdx(layerID)];
				this->update_distance(layerID, 1, apos);
			}
			this->status_change();
		}

		/**
		 * Update distance transform for all points in given layer.
		 *
		 * @param layerID
		 * @param side
		 */
		void update_distance(const INT& layerID, const INT& side)
		{
			const UINT& layerIdx = this->layerIdx(layerID);
			for (
				const VecDi& pos_child
				: m_grid_phi.branch().list(layerIdx)
			) {
				this->update_distance(
					layerID, side,
					m_grid_phi.child(pos_child).list(layerIdx)
				);
			}
		}

		/**
		 * Update distance transform for points in layer layerID given
		 * in alayer.
		 *
		 * @param layerID
		 * @param side
		 * @param alayer
		 */
		void update_distance(
			const INT& layerID, const INT& side, PosArray& alayer
		) {
			DeltaPhiGrid& dphi = this->dphi();
			const UINT& size = alayer.size();

			// Calculate distance of every point in this layer to the zero
			// layer, and store in delta phi grid.
			// Delta phi grid is used to allow for asynchronous updates, that
			// is, to prevent neighbouring points affecting the distance
			// transform.
//#pragma omp parallel for
			for (UINT pos_idx = 0; pos_idx < alayer.size(); pos_idx++)
			{
				// Current position along this layer.
				const VecDi& pos = alayer[pos_idx];
				// Distance from this position to zero layer.
				const FLOAT dist = this->distance(pos, side);
				// Update delta phi grid.
				this->dphi(pos, dist, layerID);
			}

			// Update distance in phi from delta phi and append any points that
			// move out of their layer to a status change list.

// TODO: cannot parallelise, since phi() can create new layer items for outer
// layers.
// Should split outer layer expansion to separate process.
//#pragma omp parallel for
			for (UINT pos_idx = 0; pos_idx < alayer.size(); pos_idx++)
			{
				// Current position along this layer.
				const VecDi& pos = alayer[pos_idx];
				// Distance calculated above.
				const FLOAT& dist = dphi(pos);
				// Update phi grid. Note that '=' is used here, rather than
				// '+=' as with the zero layer.
				this->phi(pos, dist, layerID);
			} // End for pos_idx.

		}


		/**
		 * Calculate city-block distance from position to zero curve.
		 *
		 * @param pos
		 * @param side
		 * @return
		 */
		FLOAT distance (const VecDi& pos, const FLOAT& side) const
		{
			const PhiGrid& phi = this->phi();
			// Get neighbouring point that is next closest to the zero-layer.
			const VecDi pos_closest = this->next_closest(pos, side);
			const FLOAT val_closest = phi(pos_closest);
			// This point's distance is then the distance of the closest
			// neighbour +/-1, depending which side of the band we are looking
			// at.
			const FLOAT dist = val_closest + side;
			return dist;
		}


		/**
		 * Find all outer layer points who's distance transform is
		 * affected by modified zero-layer points.
		 *
		 * TODO: several options for optimisation of removing duplicates:
		 * - Use a boolean flag grid to construct a de-duped vector (used
		 * here).
		 * - Check std::vector in Grid::neighs() using std::find to prevent
		 * adding a duplicate in the first place.
		 * - Use std::vector sort, unique, erase.
		 * - Use a std::unordered_set with a suitable hashing function.
		 *
		 * @param apos
		 */
		void affected(PosArray* apos)
		{
			// Reference to phi grid.
			typedef typename DeltaPhiGrid::BranchGrid DeltaPhiBranch;
			typedef typename DeltaPhiGrid::ChildGrid DeltaPhiChild;
			typedef typename AffectedLookupGrid::PosArray PosArray;

			// Add modified zero-layer points to tracking list.
			const UINT& layerIdx = this->layerIdx(0);
			for (const VecDi& pos_child : this->dphi().branch().list(layerIdx))
			{
				const PosArray& aposdphi = (
					this->dphi().child(pos_child).list(layerIdx)
				);
				for (const VecDi& pos_leaf : aposdphi)
					m_grid_affected.add(pos_leaf);
			}

			std::vector<UINT> aidx_first_neigh;
			std::vector<UINT> aidx_last_neigh;

			aidx_first_neigh.resize(m_grid_affected.branch().list().size());
			aidx_last_neigh.resize(m_grid_affected.branch().list().size());

			for (UINT udist = 1; udist <= L; udist++)
			{
				const UINT size_prev = aidx_first_neigh.size();
				const UINT size_now = m_grid_affected.branch().list().size();

				aidx_first_neigh.resize(size_now);
				aidx_last_neigh.resize(size_now);

				for (
					UINT idx_child = size_prev; idx_child < size_now;
					idx_child++
				)
					aidx_first_neigh[idx_child] = 0;

				for (UINT idx_child = 0; idx_child < size_now; idx_child++)
				{
					const PosArray& apos_child = m_grid_affected.branch().list();
					const VecDi& pos_child = apos_child[idx_child];
					aidx_last_neigh[idx_child] = (
						m_grid_affected.child(pos_child).list().size()
					);
				}

				for (UINT idx_child = 0; idx_child < size_now; idx_child++)
				{
					// Not by reference since tracking list can be resized.
					const VecDi pos_child = (
						m_grid_affected.branch().list()[idx_child]
					);

					for (
						UINT idx_neigh = aidx_first_neigh[idx_child];
						idx_neigh < aidx_last_neigh[idx_child]; idx_neigh++
					) {
						const PosArray& apos_neigh = (
							m_grid_affected.child(pos_child).list()
						);
						const VecDi& pos_neigh = apos_neigh[idx_neigh];

						this->phi().neighs(pos_neigh, &m_grid_affected);
					}
				}
				for (UINT idx = 0; idx < size_now; idx++)
					aidx_first_neigh[idx] = aidx_last_neigh[idx];
			}


//			// Cycle through neighbours out to a distance of L.
//			UINT idx_first_neigh = 0;
//			UINT idx_last_neigh = 0;
//			for (UINT udist = 1; udist <= L; udist++)
//			{
//				idx_last_neigh = aneighs.size();
//				// Cycle current subset of neighbours.
//				for (
//					UINT idx_neigh = idx_first_neigh;
//					idx_neigh < idx_last_neigh; idx_neigh++
//				) {
//					// NOTE: cannot get pos by reference, must make a copy,
//					// because phi.neighs() can cause a reallocation, making
//					// the reference invalid.
//					const VecDi pos_neigh = aneighs[idx_neigh];
//					// Append new neighbours to neighbour vector.
//					phi.neighs(pos_neigh, aneighs, m_grid_flag);
//				}
//				idx_first_neigh = idx_last_neigh;
//			}

			// Cycle de-duped set of neighbours.
			for (VecDi pos_neigh : m_grid_affected.leafs())
			{
				const INT layer_neigh = this->layerID(pos_neigh);
				// Ensure this point lies in an outer layer.
				if (this->inside_band(layer_neigh))
				{
					// Append the point index to the output list at the
					// appropriate layer index.
					apos[L + layer_neigh].push_back(pos_neigh);
				}
			}
		}

#ifndef _TESTING
	protected:
#endif

	};
}
#endif
