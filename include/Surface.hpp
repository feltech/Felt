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

	template <UINT D, UINT L=2, UINT P=2>
	class Surface
	{
	protected:
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;

		struct StatusChange
		{
			StatusChange(
				const VecDi& pos_, const INT& from_layer_, const INT& to_layer_
			) : pos(pos_), from_layer(from_layer_), to_layer(to_layer_){}

			VecDi pos;
			INT from_layer;
			INT to_layer;
		};

		typedef PartitionedArray<StatusChange, D, P>	StatusChangeGrid;

	public:

		typedef TrackedPartitionedGrid<FLOAT, D, P, 2*L+1>	DeltaPhiGrid;
		typedef SharedTrackedPartitionedGrid<FLOAT, D, P, 2*L+1>	PhiGrid;

#ifndef _TESTING
	protected:
#else
	public:
#endif

		VecDi m_pos_min;
		VecDi m_pos_max;

		Grid<bool,D> m_grid_flag;

		PhiGrid				m_grid_phi_parent;

		DeltaPhiGrid 		m_grid_dphi;

		StatusChangeGrid	m_grid_status_change;


	public:
		Surface ()
		:	m_grid_phi_parent(),
			m_grid_status_change()
		{}


		Surface (const VecDu& dims, const UINT& uborder = 0)
		:	m_grid_phi_parent(),
			m_grid_status_change()
		{
			this->dims(dims, uborder);
		}


		FLOAT operator() (const VecDi& pos)
		{
			return this->phi()(pos);
		}


		/**
		 * @brief Set dimensions and store limits adjusted for narrow band
		 * space.
		 * @param vec_dims
		 */
		void dims (const VecDu& udims, const UINT& uborder = 0)
		{
			const VecDi idims = udims.template cast<INT>();
			const VecDi offset = -idims/2;

			m_grid_phi_parent.init(udims, offset);

			// Configure delta phi embedding.
			m_grid_dphi.init(udims, offset);
			// Configure status change partitioned lists.
			m_grid_status_change.init(udims, offset);

			// Configure boolean flag grid.
			m_grid_flag.dims(udims);
			m_grid_flag.offset(offset);

			// Store min and max usable positions in phi embedding.
			this->pos_min(
				VecDi::Constant(L + uborder + 1) + m_grid_phi_parent.offset()
			);
			this->pos_max(
				(idims - VecDi::Constant(L + uborder + 1))
				+ m_grid_phi_parent.offset() - VecDi::Constant(1)
			);
			// Fill phi grid with 'outside' value.
			m_grid_phi_parent.fill(L+1);
			// Initialise delta phi to zero.
			m_grid_dphi.fill(0);
			// Initialise flag grid to false.
			m_grid_flag.fill(false);
		}

		/**
		 * @brief Get minimum usable position in phi grid.
		 * @return
		 */
		const VecDi& pos_min () const
		{
			return m_pos_min;
		}

		/**
		 * @brief Get maximum usable position in phi grid.
		 * @return
		 */
		const VecDi& pos_max () const
		{
			return m_pos_max;
		}


		/**
		 * @brief Set minimum usable position in phi grid.
		 * @param pos
		 */
		void pos_min (const VecDi& pos)
		{
			m_pos_min = pos;
		}

		/**
		 * @brief Set maximum usable position in phi grid.
		 * @param pos
		 */
		void pos_max (const VecDi& pos)
		{
			m_pos_max = pos;
		}


		/**
		 * @brief Get reference to phi grid.
		 * @return
		 */
		PhiGrid& phi ()
		{
			return m_grid_phi_parent;
		}

		/**
		 * @brief Get reference to phi grid.
		 * @return
		 */
		const PhiGrid& phi () const
		{
			return m_grid_phi_parent;
		}


		const FLOAT& phi (const VecDi& pos) const
		{
			return m_grid_phi_parent(pos);
		}

		FLOAT& phi (const VecDi& pos)
		{
			return m_grid_phi_parent(pos);
		}

		/**
		 * @brief Test whether a given value lies within the narrow band or not.
		 * @param val
		 * @return
		 */
		template <typename ValType>
		bool inside_band (const ValType& val)
		{
			return (UINT)std::abs(val) <= L;
		}


		/**
		 * @brief Update phi grid point at pos by val.
		 * Checks if the layer should change and adds to the appropriate status
		 * change list.
		 * Also expands the outer layers as the surface expands/contracts.
		 * TODO: because of the outer layer expansion code, this function is
		 * not, in general, thread safe.  Must move outer layer expansion to a
		 * separate routine.
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
			std::vector<VecDi> neighs;
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
		 * @brief Add a point to the status change list to eventually be moved
		 * from one layer to another.
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


		void status_change ()
		{
			for (
				const VecDi& pos_child : m_grid_status_change.branch().list()
			) {
				for (
					const StatusChange& change
					: m_grid_status_change.child(pos_child)
				) {

					this->layer_remove(
						change.pos, change.from_layer
					);
				}
			}
			for (
				const VecDi& pos_child : m_grid_status_change.branch().list()
			) {
				for (
					const StatusChange& change
					: m_grid_status_change.child(pos_child)
				) {

					this->layer_add(
						change.pos, change.to_layer
					);
				}
			}
		}


		/**
		 * @brief Get reference to delta phi grid.
		 * @return
		 */
		DeltaPhiGrid& dphi ()
		{
			return m_grid_dphi;
		}

		const DeltaPhiGrid& dphi () const
		{
			return m_grid_dphi;
		}


		FLOAT& dphi (const VecDi& pos)
		{
			return m_grid_dphi(pos);
		}

		const FLOAT& dphi (const VecDi& pos) const
		{
			return m_grid_dphi(pos);
		}


		/**
		 * @brief Update delta phi grid and append point to change list for
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
		 * @brief Get reference to a single layer of the narrow band.
		 * @param id
		 * @return
		 */
		typename PhiGrid::PosArray& layer (
			const VecDi& pos_child, const INT& id = 0
		) {
			return m_grid_phi_parent.child(pos_child).list(id+L);
		}

		/**
		 * @brief Get reference to a single layer of the narrow band.
		 * @param id
		 * @return
		 */
		const typename PhiGrid::PosArray& layer (
			const VecDi& pos_child, const INT& id = 0
		) const
		{
			return m_grid_phi_parent.child(pos_child).list(id+L);
		}


		const LeafsContainer<PhiGrid> layer(const UINT& layerID) const
		{
			return m_grid_phi_parent.leafs(this->layerIdx(layerID));
		}

		/**
		 * @brief Append position to a layer of the narrow band.
		 * @param id
		 * @param pos
		 */
		void layer_add (const VecDi& pos, const INT& layerID)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;
			m_grid_phi_parent.add(pos, this->layerIdx(layerID));
		}

		/**
		 * @brief Append position to a layer of the narrow band.
		 * Skip the lookup in the phi grid by assuming passed val is the phi
		 * grid value of this point.
		 * @param pos
		 * @param val
		 */
		void layer_add (const FLOAT& val, const VecDi& pos)
		{
			const INT& id = this->layerID(val);
			this->layer_add(pos, id);
		}


		void layer_remove(const VecDi& pos, const INT& layerID)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;

			m_grid_phi_parent.remove(pos, this->layerIdx(layerID));
		}


		void layer_move(
			const VecDi& pos, const INT& fromLayerID, const INT& toLayerID
		) {
			this->layer_remove(pos, fromLayerID);
			this->layer_add(pos, toLayerID);
		}

		/**
		 * @brief Get narrow band layer id of location in phi grid.
		 * @param pos
		 * @return
		 */
		INT layerID(const VecDi& pos) const
		{
			return this->layerID(this->phi(pos));
		}
		/**
		 * @brief Get narrow band layer id of value.
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
		 * @brief Get narrow band layer index of ID.
		 * @param val
		 * @return
		 */
		UINT layerIdx(const INT& id) const
		{
			return id+L;
		}

		/**
		 * @brief Create a single singularity seed point in the phi grid.
		 * NOTE: does not handle overwriting of points currently already on the
		 * surface/in the volume.
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
		 * @brief Get neighbouring position in phi grid that is closest to
		 * zero-curve.
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
			std::vector<VecDi> neighs;
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
		 * @brief Get neighbouring position in phi grid that is closest to
		 * zero-curve.
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
		 * @brief Reset delta phi to zero and clear update lists.
		 */
		void update_start ()
		{
			for (UINT layerIdx = 0; layerIdx < 2*L+1; layerIdx++)
				this->dphi().reset(0, layerIdx);

			m_grid_status_change.reset();
		}


		/**
		 * @brief Apply delta phi to phi along the zero layer.
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
		 * @brief Update zero layer then update distance transform for all
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
		 * @brief Update zero layer then update distance transform for affected
		 * points in each layer.
		 */
		void update_end_local ()
		{
			typedef typename PhiGrid::PosArray PosArray;
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
		 * @brief Update distance transform for all points in given layer
		 * @param layerID
		 * @param side
		 */
		void update_distance(const INT& layerID, const INT& side)
		{
			const UINT& layerIdx = this->layerIdx(layerID);
			for (
				const VecDi& pos_child
				: m_grid_phi_parent.branch().list(layerIdx)
			) {
				this->update_distance(
					layerID, side,
					m_grid_phi_parent.child(pos_child).list(layerIdx)
				);
			}
		}

		/**
		 * @brief Update distance transform for points in layer layerID given
		 * in alayer.
		 * @param layerID
		 * @param side
		 * @param alayer
		 */
		void update_distance(
			const INT& layerID, const INT& side,
			typename PhiGrid::PosArray& alayer
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
		 * @brief Calculate city-block distance from position to zero curve.
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
		 * @brief Find all outer layer points who's distance transform is
		 * affected by modified zero-layer points.
		 * TODO: several options for opmisation of removing duplicates:
		 * - Use a boolean flag grid to construct a de-duped vector (used
		 * here).
		 * - Check std::vector in Grid::neighs() using std::find to prevent
		 * adding a duplicate in the first place.
		 * - Use std::vector sort, unique, erase.
		 * - Use a std::unordered_set with a suitable hashing function.
		 * @param apos
		 */
		void affected(typename PhiGrid::PosArray* apos)
		{
			// Reference to phi grid.
			const PhiGrid& phi = this->phi();
			const DeltaPhiGrid& dphi = this->dphi();

			// Vector of all neighbours of all modified phi points.
			std::vector<VecDi> aneighs;

			typedef typename DeltaPhiGrid::BranchGrid DeltaPhiBranch;
			typedef typename DeltaPhiGrid::ChildGrid DeltaPhiChild;
			typedef typename DeltaPhiGrid::PosArray PosArray;
			// Loop through dphi lists, copying to neighbour list.
			const DeltaPhiBranch& branch = dphi.branch();
			const UINT& layerIdx = this->layerIdx(0);
			for (const VecDi& pos_child : branch.list(layerIdx))
			{
				const DeltaPhiChild& child = branch(pos_child);
				const PosArray& aposdphi = child.list(layerIdx);
				aneighs.insert(aneighs.end(), aposdphi.begin(), aposdphi.end());
			}

			// Loop again through dphi/neighbours, setting flag to true to
			// avoid re-visiting.
			for (UINT upos = 0; upos < aneighs.size(); upos++)
			{
				const VecDi& pos = aneighs[upos];
				m_grid_flag(pos) = true;
			}

			// Cycle through neighbours out to a distance of L.
			UINT idx_first_neigh = 0;
			UINT idx_last_neigh = 0;
			for (UINT udist = 1; udist <= L; udist++)
			{
				idx_last_neigh = aneighs.size();
				// Cycle current subset of neighbours.
				for (
					UINT idx_neigh = idx_first_neigh;
					idx_neigh < idx_last_neigh; idx_neigh++
				) {
					// NOTE: cannot get pos by reference, must make a copy,
					// because phi.neighs() can cause a reallocation, making
					// the reference invalid.
					const VecDi pos_neigh = aneighs[idx_neigh];
					// Append new neighbours to neighbour vector.
					phi.neighs(pos_neigh, aneighs, m_grid_flag);
				}
				idx_first_neigh = idx_last_neigh;
			}

			// Cycle de-duped set of neighbours.
			for (auto pos_neigh : aneighs)
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

			// Reset flag grid to all false.
			for (auto pos_neigh : aneighs)
			{
				m_grid_flag(pos_neigh) = false;
			}
		}

#ifndef _TESTING
	protected:
#endif

	};
}
#endif
