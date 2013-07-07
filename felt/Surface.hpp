#include <vector>
#include <limits>
#include "Grid.hpp"
#include <boost/math/special_functions/round.hpp>
#include <boost/numeric/ublas/assignment.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <omp.h>

namespace felt {

	template <UINT D, UINT L=2>
	class Surface {
		typedef ublas::vector<UINT,ublas::bounded_array<UINT,D> > VecDu;
		typedef ublas::vector<INT,ublas::bounded_array<INT,D> > VecDi;
		typedef ublas::vector<FLOAT,ublas::bounded_array<FLOAT,D> > VecDf;

		typedef ublas::scalar_vector<INT,ublas::bounded_array<INT,D> > ScalarDi;

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

			// Configure phi embedding.
			phi.dims(vec_dims);
			phi.offset(-VecDi(vec_dims)/2);
			// Configure delta phi embedding.
			dphi.dims(vec_dims);
			dphi.offset(-VecDi(vec_dims)/2);
			// Configure layer index spatial lookup.
			idx.dims(vec_dims);
			idx.offset(-VecDi(vec_dims)/2);
			// Store min and max usable positions in phi embedding.
			this->pos_min(ScalarDi(D,L) + phi.offset());
			this->pos_max((phi.dims() - ScalarDi(D,L)) + phi.offset() - ScalarDi(D,1));
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
					const INT dir = -sgn(newLayerID);
					std::vector<VecDi> neighs;
					phi.neighs(pos, neighs);

					for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
					{
						const VecDi& pos_neigh = neighs[neighIdx];
						const INT fromLayerID = this->layerID(pos_neigh);
						if ((UINT)std::abs(fromLayerID) > L)
						{
							const FLOAT dist_neigh = this->distance(pos_neigh, dir);
							phi(pos_neigh) = dist_neigh;
							this->status_change(pos_neigh, fromLayerID, -dir*L);
						}
					}
				}
			}
		}


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
		void dphi (const VecDi& pos, const FLOAT& val, const INT& threadIdx, const INT& layerID)
		{
			Grid<FLOAT,D>& dphi = this->dphi();
			std::vector<VecDi>& adphi = this->dphi(threadIdx, layerID);
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
		std::vector<VecDi>& layer ()
		{
			return m_layers[L];
		}


		/**
		 * @brief Append position to a layer of the narrow band.
		 * @param id
		 * @param pos
		 */
		void layer_add (const INT& id, const VecDi& pos)
		{
			if (!this->phi().inside(pos))
				throw "Adding point that is outside of grid!";
			// Do nothing if position is outside narrow band.
			if (std::fabs(id) > L)
				return;
			std::vector<VecDi>& layer = this->layer(id);
			Grid<UINT,D>& idx = this->idx();

			layer.push_back(pos);
			idx(pos) = layer.size()-1;
		}


		/**
		 * @brief Append position to a layer of the narrow band.
		 * Calculates layer ID from phi value at given position.
		 * @param pos
		 */
		void layer_add (const VecDi& pos)
		{
			const INT id = this->layerID(pos);
			this->layer_add(id, pos);
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
			if ((UINT)std::abs(fromLayerID) <= L)
				this->layer_remove(pos, fromLayerID);
			if ((UINT)std::abs(toLayerID) <= L)
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
		 * @param pos_centre
		 */
		void seed (const VecDi& pos_centre)
		{
			Grid<FLOAT,D>& phi = this->phi();
			const VecDu dims = phi.dims();

			// Width of seed.
			const ScalarDi vec_width(D, L);

			// Min and max positions affected by placing seed point.
			const VecDi pos_min = pos_centre - vec_width;
			const VecDi pos_max = pos_centre + vec_width;

			// Get vector size of window formed by pos_min and pos_max.
			const VecDi pos_size = pos_max - pos_min + ScalarDi(D,1); //+1 for zero coord.
			// Calculate number of grid points to be cycled through within window.
			UINT size = 1;
			for (UINT axis = 0; axis < dims.size(); axis++)
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
				FLOAT f_dist = (FLOAT)ublas::norm_1(vec_dist);
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
		 * @param dir
		 * @return
		 */
		VecDi next_closest (const VecDi& pos, const FLOAT& dir) const
		{
			const Grid<FLOAT,D>& phi = this->phi();
			if (this->layerID(pos) == 0)
				return pos;

			std::vector<VecDi> neighs;
			phi.neighs(pos, neighs);
			VecDi pos_nearest = VecDi(pos);
			FLOAT val_nearest = 1000;
			for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
			{
				const VecDi pos_neigh = neighs[neighIdx];
				const FLOAT val_neigh = phi(pos_neigh);
				if (val_neigh*-dir < val_nearest) {
					pos_nearest = pos_neigh;
					val_nearest = val_neigh*-dir;
				}
			}

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
			const FLOAT dir = val_centre > 0 ? -1 : 1;

			return this->next_closest(pos, dir);
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
				this->update_distance(layerID, 1);
			// Update distance transform for outer layers of the narrow band.
			for (INT layerID = 1; layerID <= (INT)L; layerID++)
				this->update_distance(layerID, -1);

			this->status_change();
		}


		void update_distance(const INT& layerID, const INT& dir)
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
				const FLOAT dist = this->distance(pos, dir);
				// Update delta phi grid.
				this->dphi(pos, dist, layerID);
			}

			// Update distance in phi from delta phi and append any points that move
			// out of their layer to a status change list.

// NOTE: cannot parallelise, since phi() can create new layer items for outer layers.
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
		 * @param dir
		 * @return
		 */
		FLOAT distance (const VecDi& pos, const FLOAT& dir) const
		{
			const Grid<FLOAT,D>& phi = this->phi();
			const VecDi pos_closest = this->next_closest(pos, dir);
			const FLOAT dist = phi(pos_closest) - dir;
			return dist;
		}

#ifndef _TESTING
	protected:
#endif

	};
}
