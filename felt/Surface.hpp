#include <vector>
#include <functional>
#include <limits>
#include <boost/math/special_functions/round.hpp>
#include <eigen3/Eigen/Dense>
#include <omp.h>

#include "Grid.hpp"


namespace felt {

	template <UINT D, UINT L=2>
	class Surface {
		typedef Eigen::Matrix<UINT, D, 1> VecDu;
		typedef Eigen::Matrix<INT, D, 1> VecDi;
		typedef Eigen::Matrix<FLOAT, D, 1> VecDf;

#ifndef _TESTING
	protected:
#else
	public:
#endif

		VecDu m_dims;
		VecDi m_offset;
		VecDi m_pos_min;
		VecDi m_pos_max;

		Grid<FLOAT,D> m_grid_phi;
		Grid<FLOAT,D> m_grid_dphi;
		Grid<UINT,D> m_grid_idx;

		std::vector<VecDi> m_layers[2*L+1];

		UINT m_uThreads;
		std::vector<std::vector<std::vector<VecDi> > > m_omp_adphi;
		std::vector<std::vector<VecDi> > m_omp_aStatusChangePos;
		std::vector<std::vector<INT> > m_omp_aStatusChangeFromLayer;
		std::vector<std::vector<INT> > m_omp_aStatusChangeToLayer;

	public:
		Surface ()
		:	m_grid_phi(),
			m_grid_idx()
		{
			this->init();
		}

		Surface (const VecDu& vec_dims)
		:	m_grid_phi(),
			m_grid_idx()
		{
			this->init();
			this->dims(vec_dims);
		}

		void init (const UINT uThreads = 0)
		{
			for (UINT layerIdx = 0; layerIdx < 2*L+1; layerIdx++)
				m_layers[layerIdx].reserve(100);

			this->num_threads(uThreads);
		}


		FLOAT operator() (const VecDi& pos)
		{
			return this->phi()(pos);
		}

		VecDi operator[] (const UINT& index)
		{
			return this->layer(0)[index];
		}


		/**
		 * @brief Set number of threads to use for OpenMP parallelisation.
		 *
		 * Felt requires data structures to store (mainly) lists of points to process (one for
		 * each narrow band layer).  These lists can be spread across multiple threads, with a
		 * separate list for each thread, so we must create them here to be used throughout.
		 *
		 * @param uThreads
		 */
		void num_threads(UINT uThreads) {
			if (uThreads == 0)
				uThreads = omp_get_max_threads();
			m_uThreads = uThreads;
			m_omp_adphi.resize(m_uThreads);
			m_omp_aStatusChangePos.resize(m_uThreads);
			m_omp_aStatusChangeFromLayer.resize(m_uThreads);
			m_omp_aStatusChangeToLayer.resize(m_uThreads);
			for (UINT threadIdx = 0; threadIdx < m_omp_adphi.size(); threadIdx++)
			{
				m_omp_adphi[threadIdx] = std::vector<std::vector<VecDi> >(2*L+1);
				m_omp_aStatusChangePos[threadIdx] = std::vector<VecDi>();
				m_omp_aStatusChangeFromLayer[threadIdx] = std::vector<INT>();
				m_omp_aStatusChangeToLayer[threadIdx] = std::vector<INT>();
			}
		}

		/**
		 * @brief Get number of threads to be used by Felt.
		 * @return
		 */
		UINT num_threads() const {
			return m_uThreads;
		}

		/**
		 * @brief Set dimensions and store limits adjusted for narrow band space.
		 * @param vec_dims
		 */
		void dims (const VecDu& vec_dims)
		{
			Grid<FLOAT,D>& phi = this->phi();
			Grid<FLOAT,D>& dphi = this->dphi();
			Grid<UINT,D>& idx = this->idx();

			const VecDi veci_dims = vec_dims.template cast<INT>();

			// Configure phi embedding.
			phi.dims(vec_dims);
			phi.offset(-VecDi(veci_dims)/2);
			// Configure delta phi embedding.
			dphi.dims(vec_dims);
			dphi.offset(-VecDi(veci_dims)/2);
			// Configure layer index spatial lookup.
			idx.dims(vec_dims);
			idx.offset(-VecDi(veci_dims)/2);
			// Store min and max usable positions in phi embedding.
			this->pos_min(VecDi::Constant(L) + phi.offset());
			this->pos_max((veci_dims - VecDi::Constant(L)) + phi.offset() - VecDi::Constant(1));
			// Fill phi grid with 'outside' value.
			phi.fill(L+1);
			idx.fill(this->null_idx());
			dphi.fill(0);
		}

		/**
		 * @brief Get dimensions.
		 * @return
		 */
		const VecDu& dims () const
		{
			return m_grid_phi.dims();
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
		Grid<FLOAT,D>& phi ()
		{
			return m_grid_phi;
		}

		/**
		 * @brief Get reference to phi grid.
		 * @return
		 */
		const Grid<FLOAT,D>& phi () const
		{
			return m_grid_phi;
		}

		/**
		 * @brief Get reference to phi grid.
		 * @return
		 */
		FLOAT phi (const VecDi& pos) const
		{
			return this->phi()(pos);
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


		void phi (const VecDi& pos, const FLOAT& val, const INT& layerID = 0)
		{
			Grid<FLOAT,D>& phi = this->phi();
			const INT newLayerID = this->layerID(val);
			phi(pos) = val;

			if (newLayerID != layerID)
			{
				this->status_change(pos, layerID, newLayerID);

				// If outside point moving inward, must create new outside points.
				if (std::abs(layerID) == L && std::abs(newLayerID) == L-1)
				{
					// Get neighbouring points.
					std::vector<VecDi> neighs;
					phi.neighs(pos, neighs);
					// Get which side of the zero-layer this point lies on.
					const INT side = sgn(newLayerID);

					for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
					{
						const VecDi& pos_neigh = neighs[neighIdx];
						const INT fromLayerID = this->layerID(pos_neigh);
						// If neighbouring point is not already within the narrow band.
						if (!this->inside_band(fromLayerID))
						{
							// Get distance of this new point to the zero layer.
							const FLOAT dist_neigh = this->distance(pos_neigh, side);
							// Set distance in phi grid.
							phi(pos_neigh) = dist_neigh;
							// Add to status change list to be added to the outer layer.
							this->status_change(pos_neigh, fromLayerID, side*L);
						}
					}
				}
			}
		}

		/**
		 * @brief Add a point to the status change list to eventually be moved from one layer to another.
		 * @param pos
		 * @param fromLayerID
		 * @param toLayerID
		 */
		void status_change (const VecDi& pos, const INT& fromLayerID, const INT& toLayerID)
		{
			m_omp_aStatusChangePos[omp_get_thread_num()].push_back(pos);
			m_omp_aStatusChangeFromLayer[omp_get_thread_num()].push_back(fromLayerID);
			m_omp_aStatusChangeToLayer[omp_get_thread_num()].push_back(toLayerID);
		}


		void status_change ()
		{
			for (UINT threadIdx = 0; threadIdx < m_omp_aStatusChangePos.size(); threadIdx++)
			{
				for (UINT posIdx = 0; posIdx < m_omp_aStatusChangePos[threadIdx].size(); posIdx++)
				{
					const VecDi& pos = m_omp_aStatusChangePos[threadIdx][posIdx];
					const INT& fromLayerID = m_omp_aStatusChangeFromLayer[threadIdx][posIdx];
					const INT& toLayerID = m_omp_aStatusChangeToLayer[threadIdx][posIdx];

					this->layer_move(pos, fromLayerID, toLayerID);
				}
			}
		}


		/**
		 * @brief Get reference to delta phi grid.
		 * @return
		 */
		Grid<FLOAT,D>& dphi ()
		{
			return m_grid_dphi;
		}

		/**
		 * @brief Get reference to delta phi grid.
		 * @return
		 */
		const Grid<FLOAT,D>& dphi () const
		{
			return m_grid_dphi;
		}

		/**
		 * @brief Get reference to delta phi array.
		 *
		 * Contains positions of non-zero points in dPhi, split across threads.
		 *
		 * @param idx - the thread index
		 * @return
		 */
		std::vector<VecDi>& dphi(const INT& threadIdx, const INT& layerID = 0)
		{
			return m_omp_adphi[threadIdx][layerID+L];
		}

		/**
		 * @brief Update delta phi grid and append point to change list for current thread.
		 * @param pos
		 * @param val
		 */
		void dphi (const UINT& uPos, const FLOAT& val, const INT& layerID = 0)
		{
			this->dphi(layer(layerID)[uPos], val, omp_get_thread_num(), layerID);
		}
		/**
		 * @brief Update delta phi grid and append point to change list for current thread.
		 * @param pos
		 * @param val
		 */
		void dphi (const VecDi& pos, const FLOAT& val, const INT& layerID = 0)
		{
			this->dphi(pos, val, omp_get_thread_num(), layerID);
		}
		/**
		 * @brief Update delta phi grid and append point to change list for given thread.
		 * @param pos
		 * @param val
		 */
		void dphi (const VecDi& pos, FLOAT val, const INT& threadIdx, const INT& layerID)
		{
			Grid<FLOAT,D>& dphi = this->dphi();
			std::vector<VecDi>& adphi = this->dphi(threadIdx, layerID);

			// If this is the zero-layer, then ensure we cannot leave the grid boundary.
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
						// Max value that will not be rounded and thus trigger a layer_move.
						const FLOAT val_max = 0.5 - (std::numeric_limits<FLOAT>::epsilon());
						// Clamp the value of delta phi.
						val = std::max(-val_max - fphi, val);
						break;
					}
			}

			dphi(pos) = val;
			adphi.push_back(pos);
		}
		/**
		 * @brief Get reference to index lookup grid.
		 * @return
		 */
		Grid<UINT,D>& idx ()
		{
			return m_grid_idx;
		}
		/**
		 * @brief Get reference to index lookup grid.
		 * @return
		 */
		const Grid<UINT,D>& idx () const
		{
			return m_grid_idx;
		}

		/**
		 * @brief Get default (null) narrow band index lookup.
		 * Used to indicate that a grid point is outside the narrow band.
		 * @return
		 */
		UINT null_idx ()
		{
			return std::numeric_limits<UINT>::max();
		}


		/**
		 * @brief Get reference to a single layer of the narrow band.
		 * @param id
		 * @return
		 */
		std::vector<VecDi>& layer (const INT& id)
		{
			return m_layers[id+L];
		}


		/**
		 * @brief Get reference to the zero-layer of the narrow band.
		 * @param id
		 * @return
		 */
		std::vector<VecDi>& layer()
		{
			return m_layers[L];
		}

		/**
		 * @brief Get reference to the zero-layer of the narrow band.
		 * @param id
		 * @return
		 */
		const std::vector<VecDi>& layer() const
		{
			return m_layers[L];
		}

		/**
		 * @brief Shortcut to zero-layer vector begin() iterator.
		 * @return
		 */
		typename std::vector<VecDi>::iterator begin()
		{
			return layer().begin();
		}

		/**
		 * @brief Shortcut to zero-layer vector end() iterator.
		 * @return
		 */
		typename std::vector<VecDi>::iterator end()
		{
			return layer().end();
		}

		/**
		 * @brief Get size, in voxels, of the surface. That is, the size of the zero-layer.
		 * @return size of zero-layer.
		 */
		UINT size()
		{
			return layer().size();
		}


		void each(std::function<void (const VecDi)> func)
		{
			std::for_each(this->begin(), this->end(), func);
		}



		/**
		 * @brief Append position to a layer of the narrow band.
		 * @param id
		 * @param pos
		 */
		void layer_add (const INT& layerID, const VecDi& pos)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;
			std::vector<VecDi>& layer = this->layer(layerID);
			Grid<UINT,D>& idx = this->idx();

			layer.push_back(pos);
			idx(pos) = layer.size()-1;
		}


		/**
		 * @brief Append position to a layer of the narrow band.
		 * Skip the lookup in the phi grid by assuming passed val is the phi grid value of this point.
		 * @param pos
		 * @param val
		 */
		void layer_add (const VecDi& pos, const FLOAT& val)
		{
			const INT id = this->layerID(val);
			this->layer_add(id, pos);
		}

		void layer_remove(const VecDi& pos, const INT& layerID)
		{
			// Do nothing if position is outside narrow band.
			if (!this->inside_band(layerID))
				return;

			std::vector<VecDi>& layer = this->layer(layerID);
			Grid<UINT,D>& grid_idx = this->idx();
			// Get index of point in layer.
			const UINT idx = grid_idx(pos);
			// Reset index lookup to null value.
			grid_idx(pos) = this->null_idx();
			const UINT layer_size = layer.size();
			if (layer_size > 1)
			{
				// Duplicate last element into this index.
				VecDi& posLast = layer[layer_size-1];
				grid_idx(posLast) = idx;
				layer[idx] = posLast;
			}
			// Pop last element.
			layer.pop_back();
		}

		void layer_move(const VecDi& pos, const INT& fromLayerID, const INT& toLayerID)
		{
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
			const Grid<FLOAT,D>& phi = this->phi();
			return this->layerID(phi(pos));
		}
		/**
		 * @brief Get narrow band layer id of value.
		 * @param val
		 * @return
		 */
		INT layerID(const FLOAT& val) const
		{
			return boost::math::round(val);
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
		 * NOTE: does not handle overwriting of points currently already on the surface/in the volume.
		 * @param pos_centre
		 */
		void seed (const VecDi& pos_centre)
		{
			Grid<FLOAT,D>& phi = this->phi();
			const VecDu dims = phi.dims();

			// Width of seed.
			const VecDi vec_width = VecDi::Constant(L);

			// Min and max positions affected by placing seed point.
			const VecDi pos_min = pos_centre - vec_width;
			const VecDi pos_max = pos_centre + vec_width;

			// Get vector size of window formed by pos_min and pos_max.
			const VecDu pos_size = (pos_max - pos_min + VecDi::Constant(1)).template cast<UINT>(); //+1 for zero coord.
			// Calculate number of grid points to be cycled through within window.
			UINT size = 1;
			for (INT axis = 0; axis < dims.size(); axis++)
				size *= pos_size(axis);

			// Cycle through each point in window.
			for (UINT u_pos = 0; u_pos <= size; u_pos++)
			{
				// Calculate vector position from integer index,
				// using Felt::Grid utility function, index().
				VecDi pos = Grid<FLOAT,D>::index(u_pos, pos_size);
				// Translate position into phi grid space.
				pos += pos_min;
				// Calculate vector distance from this position to seed centre.
				const VecDi vec_dist = pos - pos_centre;
				// Sum of absolute distance along each axis == city-block distance.
				FLOAT f_dist = (FLOAT)vec_dist.template lpNorm<1>();
				// Check distance indicates that this point is within the narrow band.
				if ((UINT)std::abs(this->layerID(f_dist)) <= L)
				{
					// Set distance as value in phi grid.
					phi(pos) = f_dist;
					// Append point to a narrow band layer (if applicable).
					this->layer_add(pos, f_dist);
				}
			}
		}

		/**
		 * @brief Get neighbouring position in phi grid that is closest to zero-curve.
		 * @param pos
		 * @param side
		 * @return
		 */
		VecDi next_closest (const VecDi& pos, const FLOAT& side) const
		{
			// Trivially return if this is already a zero-layer point.
			if (this->layerID(pos) == 0)
				return pos;

			const Grid<FLOAT,D>& phi = this->phi();

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
				// Check absolute value of this neighbour is less than nearest point.
				// NOTE: cannot simply use abs() because during an update, points close to the zero-curve
				// may end up detecting points on the other side.  By multiplying by the side (+/-1) value,
				// we ensure points on the opposite side of the band will always be considered further
				// away from the zero-layer than points on the same side.
				if (val_neigh*side < val_nearest) {
					pos_nearest = pos_neigh;
					val_nearest = val_neigh*side;
				}
			}

			// TODO: lovely elegant solution using gradient vector doesn't work in all cases because pos may not yet
			// be initialised (this next_closest() function is used to initialise it), so gradient can be erroneous.
//			const VecDf vec_dir = phi.grad(pos) * dir;
//			const UINT axis_best = ublas::index_norm_inf(vec_dir);
//			VecDi pos_nearest = VecDi(pos);
//			const FLOAT val_best = vec_dir(axis_best);
//			if (fabs(val_best) > 0)
//				pos_nearest(axis_best) += sgn(val_best);

			return pos_nearest;
		}

		/**
		 * @brief Get neighbouring position in phi grid that is closest to zero-curve.
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
			Grid<FLOAT,D>& dphi = this->dphi();
			for (UINT threadIdx = 0; threadIdx < this->num_threads(); threadIdx++)
			{
				std::vector<VecDi>& apos = this->dphi(threadIdx);
				// Reset delta phi to zero.
				for (UINT pos_idx = 0; pos_idx < apos.size(); pos_idx++)
					dphi(apos[pos_idx]) = 0;
				// Clear position list.
				apos.clear();
				// Clear status change lists.
				m_omp_aStatusChangePos[threadIdx].clear();
				m_omp_aStatusChangeFromLayer[threadIdx].clear();
				m_omp_aStatusChangeToLayer[threadIdx].clear();
			}
		}

		/**
		 * @brief Apply delta phi to phi, update distance transforms, and update narrow band.
		 */
		void update_end ()
		{
			Grid<FLOAT,D>& phi = this->phi();
			Grid<FLOAT,D>& dphi = this->dphi();
			// Loop through thread-localised update points along zero layer.
#pragma omp parallel for
			for (UINT threadIdx = 0; threadIdx < this->num_threads(); threadIdx++)
			{
				// Get zero-layer update points for this thread.
				std::vector<VecDi>& apos = this->dphi(threadIdx);
				// Loop through positions, updating phi by delta phi.
				for (UINT pos_idx = 0; pos_idx < apos.size(); pos_idx++)
				{
					const VecDi& pos = apos[pos_idx];
					const FLOAT fphi = phi(pos);
					const FLOAT fdphi = dphi(pos);
					const FLOAT fval = fphi + fdphi;
					this->phi(pos, fval);
				}
			}
			// Update distance transform for inner layers of the narrow band.
			for (INT layerID = -1; layerID >= -(INT)L; layerID--)
				this->update_distance(layerID, -1);
			// Update distance transform for outer layers of the narrow band.
			for (INT layerID = 1; layerID <= (INT)L; layerID++)
				this->update_distance(layerID, 1);

			this->status_change();
		}


		void update_distance(const INT& layerID, const INT& side)
		{
//			Grid<FLOAT,D>& phi = this->phi();
			Grid<FLOAT,D>& dphi = this->dphi();
			std::vector<VecDi>& alayer = this->layer(layerID);

			// Calculate distance of every point in this layer to the zero layer,
			// and store in delta phi grid.
			// Delta phi grid is used to allow for asynchronous updates, that is,
			// to prevent neighbouring points affecting the distance transform.
#pragma omp parallel for
			for (UINT pos_idx = 0; pos_idx < alayer.size(); pos_idx++)
			{
				// Current position along this layer.
				const VecDi& pos = alayer[pos_idx];
				// Distance from this position to zero layer.
				const FLOAT dist = this->distance(pos, side);
				// Update delta phi grid.
				this->dphi(pos, dist, layerID);
			}

			// Update distance in phi from delta phi and append any points that move
			// out of their layer to a status change list.

// TODO: cannot parallelise, since phi() can create new layer items for outer layers.
// Should split outer layer expansion to separate process.
//#pragma omp parallel for
			for (UINT pos_idx = 0; pos_idx < alayer.size(); pos_idx++)
			{
				// Current position along this layer.
				const VecDi& pos = alayer[pos_idx];
				// Distance calculated above.
				const FLOAT& dist = dphi(pos);
				// Update phi grid. Note that '=' is used here, rather than '+=' as with the zero layer.
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
			const Grid<FLOAT,D>& phi = this->phi();
			// Get neighbouring point that is next closest to the zero-layer.
			const VecDi pos_closest = this->next_closest(pos, side);
			const FLOAT val_closest = phi(pos_closest);
			// This point's distance is then the distance of the closest neighbour +/-1,
			// depending which side of the band we are looking at.
			const FLOAT dist = val_closest + side;
			return dist;
		}

#ifndef _TESTING
	protected:
#endif

	};
}
