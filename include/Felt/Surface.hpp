#ifndef Surface_hpp
#define Surface_hpp

#include <vector>
#include <functional>
#include <limits>
#include <set>
#include <boost/math/special_functions/round.hpp>
#include <eigen3/Eigen/Dense>
#include <omp.h>
#include <iostream>

#include "PartitionedGrid.hpp"

namespace felt
{

/**
 * A n-dimensional sparse-field spatially partitioned level set.
 *
 * @tparam D the number of dimensions of the surface.
 * @tparam L the number of narrow band layers surrounding the zero-level
 * surface.
 */
template <UINT D, UINT L=2>
class Surface
{
public:
	static constexpr INT LAYER_MIN	= (INT)-L;
	static constexpr INT LAYER_MAX	= (INT)L;
	static constexpr UINT NUM_LAYERS = 2*L+1;
	static constexpr FLOAT TINY = 0.00001f;
	static constexpr FLOAT REALLYTINY = std::numeric_limits<FLOAT>::epsilon();

	using Surface_t = Surface<D, L>;
	/**
	 * A delta phi update grid with active (non-zero) grid points tracked.
	 */
	using DeltaPhiGrid = TrackedPartitionedGrid<FLOAT, D, NUM_LAYERS>;
	/**
	 * A level set embedding phi grid, with active grid points (the narrow
	 * band) tracked.
	 */
	using PhiGrid = SharedTrackedPartitionedGrid<FLOAT, D, NUM_LAYERS>;
	/**
	 * A resizable array of D-dimensional grid positions.
	 */
	using PosArray = typename PhiGrid::PosArray;
	/**
	 * D-dimensional unsigned int vector.
	 */
	using VecDu = typename PhiGrid::VecDu;
	/**
	 * D-dimensional integer vector.
	 */
	using VecDi = typename PhiGrid::VecDi;
	/**
	 * D-dimensional float vector.
	 */
	using VecDf = typename PhiGrid::VecDf;
	/**
	 * Grid to track positions that require an update.
	 */
	using AffectedLookupGrid = LookupPartitionedGrid<D, 2*L+1>;


	using Plane = Eigen::Hyperplane<FLOAT, D>;
	using Line = Eigen::ParametrizedLine<FLOAT, D>;


	/**
	 * Storage class for flagging a point in the grid to be moved from one
	 * narrow band layer to another.
	 */
	struct StatusChange
	{
		StatusChange(
			const VecDi& pos_, const INT& from_layer_, const INT& to_layer_
		) : pos(pos_), from_layer(from_layer_), to_layer(to_layer_) {}

		static UINT layer_idx(const INT& layer_id)
		{
			return layer_id + (L + 1);
		}

		static INT layer_id(const INT& layer_idx)
		{
			return layer_idx - (L + 1);
		}

		VecDi pos;
		INT from_layer;
		INT to_layer;
	};

protected:
	struct ChildHit
	{
		ChildHit(const VecDf& pos_intersect_, const VecDi& pos_child_)
		: pos_intersect(pos_intersect_), pos_child(pos_child_)
		{}
		VecDf pos_intersect;
		VecDi pos_child;
	};

	/**
	 * A spatially partitioned array of StatusChange objects.
	 *
	 * There is one extra 'layer' for status change lists, one either side
	 * of the narrow band, for those points that are going out of scope.
	 */
	using StatusChangeGrid = PartitionedArray<StatusChange, D>;

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

	Surface () = default;

	/**
	 * Construct a level set embedding of size dims.
	 *
	 * All points will be marked as outside the surface (i.e. no surface).
	 *
	 * @param dims
	 * @param uborder
	 */
	Surface (
		const VecDu& dims_,
		const VecDu& dims_partition_ = VecDu::Constant(DEFAULT_PARTITION)
	)
	:	m_grid_phi(),
		m_grid_status_change(),
		m_grid_dphi()
	{
		this->init(dims_, dims_partition_);
	}

	/**
	 * Initialise a Surface (e.g. after default-construction).
	 *
	 * @param dims
	 * @param dims_partition
	 */
	void init (
		const VecDu& dims_, const VecDu& dims_partition_ = VecDu::Constant(DEFAULT_PARTITION)
	) {
		this->dims(dims_, dims_partition_);
	}

	/**
	 * Get null position vector for given template typename.
	 *
	 * TODO: c++14 should support variable templates, which is a better solution,
	 * but gcc doesn't yet.
	 *
	 * @return
	 */
	template <typename T>
	constexpr const Eigen::Matrix<T, D, 1> NULL_POS() const
	{
		return Eigen::Matrix<T, D, 1>::Constant(std::numeric_limits<T>::max());
	};

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
	bool inside_band (const ValType& val) const
	{
		return (UINT)std::abs(val) <= L;
	}

	const AffectedLookupGrid& affected ()
	{
		return m_grid_affected;
	}

	/**
	 * Update phi grid point at pos by val.
	 *
	 * Checks if the layer should change and if so adds to the appropriate
	 * status change list.
	 *
	 * Also expands the outer layers as the surface expands/contracts.
	 *
	 * @param pos
	 * @param val
	 * @param layer_id
	 */
	void phi (const VecDi& pos_, const FLOAT& val_, const INT& layer_id_ = 0)
	{
		const INT newlayer_id = this->layer_id(val_);

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		
		if (
			std::abs(this->layer_id(val_)) != std::abs(layer_id_)
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) + 1
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) - 1
		) {
			std::stringstream strs;
			strs << "Phi update value out of bounds at layer "
				<< layer_id_ << " " << felt::format(pos_) << ": " << val_;
			std::string str = strs.str();
			throw std::domain_error(str);
		}
		
		#endif

		m_grid_phi(pos_) = val_;

		this->status_change(pos_, layer_id_, newlayer_id);
	}


	void add_neighs(const VecDi& pos, const INT& side)
	{
		// Get neighbouring points.
		PosArray neighs;
		m_grid_phi.neighs(pos, neighs);

		const INT& to_layer_id = side*INT(L);

		// Get a mutex lock on affected spatial partitions, since neighbouring points might spill
		// over into another thread's partition.
		std::set<std::mutex*> locks = this->lock_children(neighs);

		for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
		{
			const VecDi& pos_neigh = neighs[neighIdx];
			typename PhiGrid::Child& child = m_grid_phi.child(m_grid_phi.pos_child(pos_neigh));
			
			const INT from_layer_id = this->layer_id(pos_neigh);

			// Only add if neighbouring point is not already within the narrow band.
			if (this->inside_band(from_layer_id))
				continue;

			// Get distance of this new point to the zero layer.
			const FLOAT& dist_neigh = this->distance(pos_neigh, side);

			#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

			if (std::abs(this->layer_id(dist_neigh)) != L)
			{
				std::stringstream strs;
				strs << "Expanding outer layer distance out of bounds " << felt::format(pos)
					<< ": " << dist_neigh;
				std::string str = strs.str();
				throw std::domain_error(str);
			}

			#endif

			// Add to status change list to be added to the outer layer.
			m_grid_status_change.add(
				pos_neigh, StatusChange(pos_neigh, from_layer_id, to_layer_id)
			);
			// Set distance in phi grid.
			child(pos_neigh) = dist_neigh;
		} // End for neighbouring points.

		for (std::mutex* lock : locks)
			lock->unlock();
	}


	const StatusChangeGrid& status_change () const
	{
		return m_grid_status_change;
	}


	/**
	 * Add a point to the status change list to eventually be moved
	 * from one layer to another.
	 *
	 * @param pos
	 * @param fromlayer_id
	 * @param tolayer_id
	 */
	bool status_change (const VecDi& pos_, const INT& layer_id_from_, const INT& layer_id_to_)
	{
		if (layer_id_from_ == layer_id_to_)
			return false;

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		{
			m_grid_phi.assert_pos_bounds(pos_, "status_change: ");

			if (
				std::abs(layer_id_from_) != std::abs(layer_id_to_) + 1
				&& std::abs(layer_id_from_) != std::abs(layer_id_to_) - 1
			)
				throw std::domain_error(std::string("Bad status_change"));
		}
		#endif

		m_grid_status_change.add(pos_, StatusChange(pos_, layer_id_from_, layer_id_to_));

		return true;
	}

	/**
	 * Loop through the status change lists moving the referenced points
	 * from one layer to another.
	 */
	void flush_status_change ()
	{
		for (const VecDi& pos_child : m_grid_status_change.branch().list())
			for (const StatusChange& status : m_grid_status_change.child(pos_child))
			{
				const INT& layer_id_to = status.to_layer;
				const INT& layer_id_from = status.from_layer;
				const VecDi& pos_leaf = status.pos;

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (
					std::abs(layer_id_from) != std::abs(layer_id_to) + 1
					&& std::abs(layer_id_from) != std::abs(layer_id_to) - 1
				) {
					std::stringstream strs;
					strs << "flush_status_change layer move invalid "
						<< felt::format(pos_leaf) << ": " << layer_id_from
						<< " -> " << layer_id_to;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				this->layer_move(pos_leaf, layer_id_from, layer_id_to);
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
	void dphi (const VecDi& pos_, FLOAT val_, const INT& layer_id_ = 0)
	{
		if (layer_id_ == 0)
			val_ = clamp(pos_, val_);

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

		if (
			std::abs(this->layer_id(val_)) != std::abs(layer_id_)
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) + 1
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) - 1
		) {
			std::stringstream strs;
			strs << "Outer layer distance update value out of bounds at layer "
				<< layer_id_ << " " << felt::format(pos_) << ": " << val_;
			std::string str = strs.str();
			throw std::domain_error(str);
		}

		#endif
		
		this->dphi().add(pos_, val_, this->layer_idx(layer_id_));
	}

	/**
	 * Clamp delta phi value such that it doesn't breach the grid or cause
	 * instability.
	 *
	 * @param pos
	 * @param val
	 * @return clamped val
	 */
	FLOAT clamp (const VecDi& pos, FLOAT val)
	{
		const VecDi& pos_min = this->pos_min();
		const VecDi& pos_max = this->pos_max();

		if (std::abs(val) > 1.0f)
			val = sgn(val);

		// Cycle each axis.
		for (UINT d = 0; d < D; d++)
			// Check if pos lies at the max bound of this axis.
			if (pos_min(d) == pos(d) || pos_max(d) == pos(d))
			{
				// Get phi at this point.
				const FLOAT& fphi = this->phi(pos);
				// Max value that will not be rounded and thus trigger
				// a layer_move.
				const FLOAT& val_max = -0.5 +
					(std::numeric_limits<FLOAT>::epsilon() * 2);
				// Clamp the value of delta phi.
				val = std::max(val_max - fphi, val);
				break;
			}
		return val;
	}

	/**
	 * Update delta phi grid surrounding a given point with amount determined
	 * by Gaussian distribution about (real-valued) central point.
	 *
	 * Will be normalised so that the total amount distributed over the points
	 * sums to val. Locks mutexes of affected spatial partitions so thread safe.
	 *
	 * @param dist distance from central point to spread over.
	 * @param pos_centre centre of the Gaussian distribution.
	 * @param val amount to spread over the points.
	 * @param stddev standard deviation of Gaussian.
	 */	
	template <UINT Distance>
	FLOAT dphi_gauss (
		const VecDf& pos_centre, const FLOAT& val, const FLOAT& stddev
	) {
		const SharedLookupGrid<D, NUM_LAYERS>& lookup = this->walk_band<Distance>(
			round(pos_centre)
		);
		return this->dphi_gauss(lookup.list(this->layer_idx(0)), pos_centre, val, stddev);

//		PosArray list;
//		const INT w = Distance;
//		const VecDi offset = VecDi::Constant(-w/2);
//		const VecDu dims = VecDu::Constant(w);
//		const VecDi pos_centre_rounded = round(pos_centre);
//		for (UINT i = 0; i < std::pow(w, D); i++)
//		{
//			const VecDi& pos_offset = Grid<UINT, D>::index(i, dims, offset);
//			const VecDi& pos_neigh = pos_offset + pos_centre_rounded;
//			if (this->layer_id(pos_neigh) == 0)
//				list.push_back(pos_neigh);
//		}
//
//		return this->dphi_gauss(list, pos_centre, val, stddev);
	}

	std::set<std::mutex*> lock_children(const PosArray& list_pos_leafs_)
	{
		std::set<std::mutex*> locks;
		for (const VecDi& pos : list_pos_leafs_)
			locks.insert(&m_grid_dphi.child(m_grid_dphi.pos_child(pos)).lookup().mutex());

		for (std::mutex* lock : locks)
			lock->lock();

		return locks;
	}

	/**
	 * Update delta phi grid at given points with amount determined by Gaussian
	 * distribution about (real-valued) central point.
	 *
	 * Will be normalised so that the total amount distributed over the points
	 * sums to val. Locks mutexes of affected spatial partitions so thread safe.
	 *
	 * TODO: better handling of overlapping Gaussians - current just clamps,
	 * so can lose mass.
	 *
	 * @param list the points to spread the amount over.
	 * @param pos_centre centre of the Gaussian distribution.
	 * @param val amount to spread over the points.
	 * @param stddev standard deviation of Gaussian.
	 */
	FLOAT dphi_gauss (
		const typename DeltaPhiGrid::PosArray& list, const VecDf& pos_centre,
		const FLOAT& val, const FLOAT& stddev
	) {
		constexpr FLOAT sqrt2piDinv = 1.0f/sqrt(pow(2*M_PI, D));

		Eigen::VectorXf weights(list.size());
		std::set<std::mutex*> locks = this->lock_children(list);

		// Calculate Gaussian weights.
		for (UINT idx = 0; idx < list.size(); idx++)
		{
			const VecDf& pos = list[idx].template cast<FLOAT>();
			if (this->dphi(list[idx]) != 0)
			{
				weights(idx) = 0;
				continue;
			}

			const FLOAT& dist_sq = (pos - pos_centre).squaredNorm();
			const FLOAT& weight = sqrt2piDinv * exp(-0.5f * dist_sq);
			weights(idx) = weight;
		}

		const FLOAT& weights_sum = weights.sum();
		if (weights_sum > 0)
		{
			// Normalise weights
			weights *= val / weights_sum;

			for (UINT idx = 0; idx < list.size(); idx++)
			{
				const VecDi& pos = list[idx];
				const FLOAT amount = this->clamp(pos, weights(idx) + m_grid_dphi(pos));

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (std::abs(amount) > 1)
				{
					std::stringstream strs;
					strs << "Delta phi value out of bounds at " << felt::format(pos)
						<< ": " << amount;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				if (this->layer_id(pos) != 0)
				{
					std::stringstream strs;
					strs << "Attempting to Gaussian update non-zero layer " << felt::format(pos)
						<< ": " << m_grid_phi(pos) << " += " << amount;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				
				#endif
				
				weights(idx) -= (weights(idx) - amount);

				this->dphi().add(pos, amount, this->layer_idx(0));
			}
		}
		for (std::mutex* lock : locks)
			lock->unlock();

		return val - weights.sum();
	}

	/**
	 * Update delta phi grid surrounding a zero layer point found via a
	 * raycast by amount spread using Gaussian distribution.
	 *
	 * @param pos_origin
	 * @param dir
	 * @param dist
	 * @param amount
	 * @param stddev
	 */
	template <UINT Distance>
	FLOAT dphi_gauss (
		const VecDf& pos_origin, const VecDf& dir,
		const FLOAT& val, const float& stddev
	) {
		const VecDf& pos_hit = this->ray(pos_origin, dir);

		if (pos_hit == NULL_POS<FLOAT>())
		{
			std::cout << "MISS" << std::endl;
			return val;
		}
		return this->dphi_gauss<Distance>(pos_hit, val, stddev);
	}

	/**
	 * Get reference to the active partitions of a given layer of the narrow
	 * band.
	 *
	 * @param layer_id
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
	 * @param layer_id
	 * @return
	 */
	const LeafsContainer<PhiGrid> layer(const UINT& layer_id) const
	{
		return m_grid_phi.leafs(this->layer_idx(layer_id));
	}

	/**
	 * Append position to a layer of the narrow band.
	 *
	 * @param id
	 * @param pos
	 */
	void layer_add (const VecDi& pos, const INT& layer_id)
	{
		// Do nothing if position is outside narrow band.
		if (!this->inside_band(layer_id))
			return;
		m_grid_phi.add(pos, this->layer_idx(layer_id));
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
		const INT& id = this->layer_id(val);
		this->layer_add(pos, id);
	}

	/**
	 * Remove a position from a given layer of the narrow band.
	 *
	 * Does not modify underlying grid value, just the layer list.
	 *
	 * @param pos
	 * @param layer_id
	 */
	void layer_remove(const VecDi& pos, const INT& layer_id)
	{
		// Do nothing if position is outside narrow band.
		if (!this->inside_band(layer_id))
			return;

		m_grid_phi.remove(pos, this->layer_idx(layer_id));
	}

	/**
	 * Remove a point from one layer and move it into another.
	 *
	 * Simply adjusts the lists, does not modify the underlying grid values.
	 *
	 * @param pos
	 * @param fromlayer_id
	 * @param tolayer_id
	 */
	void layer_move(
		const VecDi& pos, const INT& fromlayer_id, const INT& tolayer_id
	) {
		this->layer_remove(pos, fromlayer_id);
		this->layer_add(pos, tolayer_id);
	}

	/**
	 * Get narrow band layer id of location in phi grid.
	 * @param pos
	 * @return
	 */
	template <class PosType>
	INT layer_id(const PosType& pos) const
	{
		return this->layer_id(this->phi()(pos));
	}

	/**
	 * Get narrow band layer id of value.
	 * @param val
	 * @return
	 */
	INT layer_id(const FLOAT& val) const
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
	UINT layer_idx(const INT& id) const
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
			if ((UINT)std::abs(this->layer_id(f_dist)) <= L)
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
		if (this->layer_id(pos) == 0)
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
		const PhiGrid& phi = this->phi();
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
		for (UINT layer_idx = 0; layer_idx < 2*L+1; layer_idx++)
		{
			m_grid_dphi.reset(0, layer_idx);
			m_grid_affected.reset(layer_idx);
		}

		m_grid_status_change.reset();
		
		#if false
		
		for (const VecDi& pos_child : m_grid_dphi.branch())
			for (const VecDi& pos : m_grid_dphi.child(pos_child))
			{
				if (m_grid_dphi(pos) != 0)
				{
					throw std::domain_error(std::string("Delta phi not reset!"));
				}
			}
		#endif
	}


	/**
	 * Apply delta phi to phi along the zero layer.
	 */
	void update_zero_layer ()
	{
		const typename DeltaPhiGrid::BranchGrid&
		branch = this->dphi().branch();

		const UINT& layer_idx = this->layer_idx(0);

		#pragma omp parallel for
		for (
			UINT idx_child = 0; idx_child < branch.list(layer_idx).size();
			idx_child++
		) {
			const VecDi& pos_child = branch.list(layer_idx)[idx_child];
			for (const VecDi& pos : branch(pos_child).list(layer_idx))
			{
				const FLOAT& fphi = this->phi(pos);
				const FLOAT& fdphi = this->dphi(pos);
				const FLOAT& fval = fphi + fdphi;

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (this->layer_id(fphi) != 0)
				{
					std::stringstream strs;
					strs << "Zero layer updated attempted at non-zero layer point "
						<< felt::format(pos) << ": " << fphi << " + " << fdphi << " = " << fval;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				if (std::abs(this->layer_id(fval)) > 1)
				{
					std::stringstream strs;
					strs << "Zero layer phi value out of bounds at " << felt::format(pos)
						<< ": " << fphi << " + " << fdphi << " = " << fval;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

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

		this->update_distance(m_grid_phi);

		this->flush_status_change();
	}

	/**
	 * Update zero layer then update distance transform for affected points in each layer.
	 */
	void update_end_local ()
	{
		// Get points in outer layers that are affected by changes in
		// zero-layer.
		this->calc_affected();

		// Update the zero layer, applying delta to phi.
		this->update_zero_layer();

		this->update_distance(m_grid_affected);

		this->flush_status_change();
	}

	template <typename GridType>
	void update_distance(const GridType& lookup_)
	{
		// Update distance transform for inner layers of the narrow band.
		for (INT layer_id = -1; layer_id >= LAYER_MIN; layer_id--)
			this->update_distance(layer_id, -1, lookup_);

		// Update distance transform for outer layers of the narrow band.
		for (INT layer_id = 1; layer_id <= LAYER_MAX; layer_id++)
			this->update_distance(layer_id, 1, lookup_);
	}

	/**
	 * Update distance transform for all points in given layer.
	 *
	 * @param layer_id
	 * @param side
	 */
	template <typename GridType>
	void update_distance(const INT& layer_id_, const INT& side_, const GridType& lookup_)
	{
		using DeltaPhiChild = typename DeltaPhiGrid::Child;
		using PhiChild = typename PhiGrid::Child;
		const UINT& layer_idx = this->layer_idx(layer_id_);

		#pragma omp parallel for
		for (UINT pos_idx = 0; pos_idx < lookup_.branch().list(layer_idx).size(); pos_idx++)
		{
			const VecDi& pos_child = lookup_.branch().list(layer_idx)[pos_idx];

			m_grid_dphi.add_child(pos_child, layer_idx);

			PhiChild& grid_phi_child = m_grid_phi.branch().get(pos_child);
			DeltaPhiChild& grid_dphi_child = m_grid_dphi.branch().get(pos_child);

			const PosArray& apos_leafs = lookup_.child(pos_child).list(layer_idx);

			// Calculate distance of every point in this layer to the zero
			// layer, and store in delta phi grid.
			// Delta phi grid is used to allow for asynchronous updates, that
			// is, to prevent neighbouring points affecting the distance
			// transform.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance from this position to zero layer.
				const FLOAT dist = this->distance(pos, side_);

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (
					std::abs(this->layer_id(dist)) != std::abs(layer_id_)
					&& std::abs(this->layer_id(dist)) != std::abs(layer_id_) + 1
					&& std::abs(this->layer_id(dist)) != std::abs(layer_id_) - 1
				) {
					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds at layer "
						<< layer_id_ << " " << felt::format(pos) << ": " << dist;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				// Update delta phi grid.
				grid_dphi_child.add(pos, dist, layer_idx);
			}
		}
		
		#pragma omp parallel for
		for (UINT pos_idx = 0; pos_idx < lookup_.branch().list(layer_idx).size(); pos_idx++)
		{
			const VecDi& pos_child = lookup_.branch().list(layer_idx)[pos_idx];

			PhiChild& grid_phi_child = m_grid_phi.branch().get(pos_child);
			DeltaPhiChild& grid_dphi_child = m_grid_dphi.branch().get(pos_child);

			const PosArray& apos_leafs = lookup_.child(pos_child).list(layer_idx);

			// Update distance in phi from delta phi and append any points that
			// move out of their layer to a status change list.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance calculated above.
				const FLOAT& dist = grid_dphi_child(pos);

				grid_phi_child(pos) = dist;

				const INT& newlayer_id = this->layer_id(dist);

				if (this->status_change(pos, layer_id_, newlayer_id))
				{
					// If outside point moving inward, must create new outside points.
					if (std::abs(layer_id_) == L && std::abs(newlayer_id) == L-1)
					{
						// Get which side of the zero-layer this point lies on.
						const INT& side = sgn(newlayer_id);
						this->add_neighs(pos, side);
					}
				}
			} // End for pos_idx.
		}
	}


	/**
	 * Calculate city-block distance from position to zero curve.
	 *
	 * @param pos
	 * @param side
	 * @return
	 */
	FLOAT distance (const VecDi& pos_, const FLOAT& side_) const
	{
		const PhiGrid& phi = this->phi();
		// Get neighbouring point that is next closest to the zero-layer.
		const VecDi pos_closest = this->next_closest(pos_, side_);
		const FLOAT val_closest = phi(pos_closest);
		// This point's distance is then the distance of the closest
		// neighbour +/-1, depending which side of the band we are looking
		// at.
		const FLOAT dist = val_closest + side_;
		return dist;
	}

	struct SortablePos
	{
		VecDi pos;
		UINT hash;
		SortablePos(const VecDi& pos_, const UINT& distance_)
		: pos(pos_), hash(0)
		{
			hash = Grid<UINT, D>::index(
				pos,
				VecDu::Constant(distance_ * 2 + 1),
				pos - VecDi::Constant(distance_)
			);
		}
		bool operator<(const SortablePos& other)
		{
			return this->hash < other.hash;
		}
	};

	/**
	 * Walk the narrow band from given position out to given distance.
	 *
	 * @param pos_
	 * @param distance_
	 * @return SharedLookupGrid with tracking lists for visited points, one list
	 * for each layer.
	 */
	template <UINT Distance>
	SharedLookupGrid<D, NUM_LAYERS>& walk_band (const VecDi& pos_)
	{
		using Lookup = SharedLookupGrid<D, 2*L+1>;
		using PosArray = typename Lookup::PosArray;

		// Box size is: (Distance * 2 directions) + 1 central.
		static Lookup lookup(VecDu::Constant(Distance * 2 + 1));
		lookup.reset_all();
		lookup.offset(pos_ - VecDi::Constant(Distance));

		const INT& layer_id = this->layer_id(pos_);
		if (!this->inside_band(layer_id))
			return lookup;

		lookup.add(pos_, this->layer_idx(layer_id));

		// Arrays to store first and last element in tracking list within
		// each spatial partition of tracking grid.
		std::array<UINT, NUM_LAYERS> aidx_first_neigh;
		std::array<UINT, NUM_LAYERS> aidx_last_neigh;
		aidx_first_neigh.fill(0);

		// Loop round searching outward for zero layer grid nodes.
		for (UINT udist = 1; udist <= Distance; udist++)
		{
			// Copy number of active grid nodes for this partition
			// into relevant index in the list.
			for (UINT i = 0; i < NUM_LAYERS; i++)
				aidx_last_neigh[i] = lookup.list(i).size();

			for (UINT layer_idx = 0; layer_idx < NUM_LAYERS; layer_idx++)
			{
				// Loop over leaf grid nodes, using the cached start and end
				// indices, so that newly added points are skipped.
				for (
					UINT idx_neigh = aidx_first_neigh[layer_idx];
					idx_neigh < aidx_last_neigh[layer_idx];
					idx_neigh++
				) {
					// This leaf grid nodes is the centre to search
					// about.
					const VecDi& pos_centre = lookup.list(layer_idx)[idx_neigh];

					// Use utility method from Grid to get neighbouring
					// grid nodes and call a lambda to add the point
					// to the appropriate tracking list.

					this->phi().neighs(
						pos_centre,
						[this](const VecDi& pos_neigh) {
							// Calculate layer of this neighbouring
							// point from the phi grid.
							const INT& layer_id = this->layer_id(pos_neigh);
							// If the calculated layer lies within the
							// narrow band, then we want to track it.
							if (this->inside_band(layer_id))
							{
								// Add the neighbour point to the
								// tracking grid. Will reject
								// duplicates.
								lookup.add(
									pos_neigh, this->layer_idx(layer_id)
								);
							}
						}
					); // End neighbourhood query.
				} // End for leaf grid node.
			}
			// Set first index to be the previous last index, so we start there
			// on next loop around.
			for (UINT i = 0; i < NUM_LAYERS; i++)
				aidx_first_neigh[i] = aidx_last_neigh[i];
		}

		return lookup;
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
	void calc_affected()
	{
		const UINT& layer_idxZero = this->layer_idx(0);

		// Loop over delta phi modified zero-layer points adding to
		// tracking grid.

		// Loop spatial partitions of dphi for zero-layer.
		for (
			const VecDi& pos_child
			: this->dphi().branch().list(layer_idxZero))
		{
			// Loop leaf grid nodes with spatial partition
			for (
				const VecDi& pos_leaf
				: this->dphi().child(pos_child).list(layer_idxZero)
			)
				// Add zero-layer point to tracking grid.
				m_grid_affected.add(pos_leaf, layer_idxZero);
		}

		// Arrays to store first and last element in tracking list within
		// each spatial partition of tracking grid.
		std::array<std::vector<UINT>, 2*L+1> aidx_first_neigh;
		std::array<std::vector<UINT>, 2*L+1> aidx_last_neigh;

		// Loop round L times, searching outward for affected outer layer
		// grid nodes.
		for (UINT udist = 1; udist <= L; udist++)
		{
			// Reset the first and last element indices for each
			// spatial partition in each layer.

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT& layer_idx = this->layer_idx(layer_id);
				// Get number of spatial partitions for this layer.
				const UINT num_childs = (
					m_grid_affected.branch().list(layer_idx).size()
				);
				// Resize spatial partition index lists for this layer to
				// to include any newly added partitions.
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
					const VecDi& pos_child = m_grid_affected.branch().list(layer_idx)[idx_child];
					// Copy number of active grid nodes for this partition
					// into relevant index in the list.
					aidx_last_neigh[layer_idx][idx_child] = (
						m_grid_affected.child(pos_child).list(layer_idx).size()
					);
				}
			}

			// Loop each layer finding the affected outer layer points
			// for each partition using the start and end points cached
			// above.

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT& layer_idx = this->layer_idx(layer_id);

				// Loop over spatial partitions, ignoring newly added ones
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
						m_grid_affected.branch().list(layer_idx)[idx_child]
					);

					// Loop over leaf grid nodes within this spatial
					// partition, using the cached start and end indices,
					// so that newly added points are skipped.
					for (
						UINT idx_neigh = aidx_first_neigh[layer_idx][idx_child];
						idx_neigh < aidx_last_neigh[layer_idx][idx_child];
						idx_neigh++
					) {
						// Get list of active leaf grid nodes in this
						// spatial partition.
						const PosArray& apos_neigh = (
							m_grid_affected.child(pos_child).list(layer_idx)
						);
						// This leaf grid nodes is the centre to search
						// about.
						const VecDi& pos_centre = apos_neigh[idx_neigh];

						// Use utility method from Grid to get neighbouring
						// grid nodes and call a lambda to add the point
						// to the appropriate tracking list.

						this->phi().neighs(
							pos_centre,
							[this](const VecDi& pos_neigh) {
								// Calculate layer of this neighbouring
								// point from the phi grid.
								const INT& layer_id = this->layer_id(pos_neigh);
								// If the calculated layer lies within the
								// narrow band, then we want to track it.
								if (this->inside_band(layer_id))
								{
									// Add the neighbour point to the
									// tracking grid. Will reject
									// duplicates.
									this->m_grid_affected.add(pos_neigh, this->layer_idx(layer_id));
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
				const UINT& layer_idx = this->layer_idx(layer_id);
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
	 * Perform a a full (parallelised) update of the narrow band.
	 *
	 * Lambda function passed will be given the position to process and
	 * a reference to the phi grid, and is expected to return delta phi to
	 * apply.
	 *
	 * @param fn_ (pos, phi) -> float
	 */
	void update(std::function<FLOAT(const VecDi&, const PhiGrid&)> fn_)
	{
		this->update_start();
		#pragma omp parallel for
		for (UINT part_idx = 0; part_idx < parts().size(); part_idx++)
		{
			const VecDi& pos_part = this->parts()[part_idx];
			for (const VecDi& pos : this->layer(pos_part))
				this->dphi(pos, fn_(pos, m_grid_phi));
		}
		this->update_end();
	}
	
	/**
	 * Get array of offsets to corners of a cube (e.g. 8 for 3D, 4 for 2D).
	 *
	 * @return
	 */
	constexpr std::array<VecDi, 1 << D> corners () const
	{
		std::array<VecDi, 1 << D> acorner;
		for (UINT mask = 0; mask < acorner.size(); mask++)
			for (UINT dim = 0; dim < D; dim++)
				if (mask & (1 << dim))
					acorner[mask](dim) = 0;
				else
					acorner[mask](dim) = -1;
		return acorner;
	}

	/**
	 * Cast a ray to the zero layer.
	 *
	 * @param pos_origin
	 * @param dir
	 * @return NULL_POS if no hit, otherwise interpolated position on zero
	 * curve.
	 */
	const VecDf ray(const VecDf& pos_origin, const VecDf& dir) const
	{

		using ChildHits = std::vector<ChildHit>;

		// If ray is cast from within phi grid, first check child grid containing origin point.
		if (m_grid_phi.inside(pos_origin))
		{
			const VecDf& pos_hit = ray(
				pos_origin, dir,
				m_grid_phi.child(
					m_grid_phi.pos_child(pos_origin.template cast<INT>())
				)
			);
			if (pos_hit != NULL_POS<FLOAT>())
				return pos_hit;
		}

		// Ray to test against.
		Line line(pos_origin, dir);

		// Tracking list for child grids that are hit.
		ChildHits child_hits;

		// Cycle each axis, casting ray to child grid planes marching away from origin.
		for (UINT dim = 0; dim < dir.size(); dim++)
		{
			// Direction +/-1 along this axis.
			FLOAT dir_dim = sgn(dir(dim));
			if (dir_dim == 0)
				continue;

			// Get next child plane along this axis.
			FLOAT pos_plane_dim = round_to_next_child(dim, dir_dim, pos_origin(dim));

			// Construct vector with elements not on this axis at zero.
			VecDf pos_plane = VecDf::Constant(0);
			pos_plane(dim) = pos_plane_dim;

			// If the zero point on this plane is not within the grid, then jump to max/min point
			// on phi grid.
			if (!m_grid_phi.inside(pos_plane))
			{
				FLOAT pos_grid_dim;
				// If casting in -'ve direction, get maximum extent.
				if (dir_dim == -1)
				{
					pos_grid_dim = m_grid_phi.offset()(dim) + m_grid_phi.dims()(dim);
					if (pos_plane_dim < pos_grid_dim)
						continue;
				}
				// Else if casting in +'ve direction, get minimum extent.
				else
				{
					pos_grid_dim = m_grid_phi.offset()(dim);
					if (pos_plane_dim > pos_grid_dim)
						continue;
				}
				// Reset plane position to max/min extent as calculated above.
				pos_plane(dim) = pos_grid_dim;
			}

			// Plane normal is opposite to ray direction.
			VecDf normal = VecDf::Constant(0);
			normal(dim) = -dir_dim;

			// Cast ray to plane and add any child grids hit on the way to tracking list.
			if (!ray_check_add_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
				continue;

			// Round up/down to next child, in case we started at inexact modulo of child grid size
			// (i.e. when phi grid size is not integer multiple of child grid size).
			pos_plane_dim = round_to_next_child(dim, dir_dim, pos_plane(dim));
			// If rounding produced a different plane, then cast to that plane, and potentially add
			// child grid to tracking list.
			if (pos_plane_dim != pos_plane(dim))
			{
				pos_plane(dim) = pos_plane_dim;
				if (!ray_check_add_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
					continue;
			}

			// Keep marching along planes, casting ray to each and adding any candidate child
			// grids to the tracking list.
			const FLOAT& child_size_dim = m_grid_phi.child_dims()(dim);
			while (true)
			{
				pos_plane(dim) += dir_dim * child_size_dim;
				if (!ray_check_add_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
					break;
			}
		}
		// Sort candidate child grids in distance order from front to back.
		std::sort(
			child_hits.begin(), child_hits.end(),
			[&pos_origin](const ChildHit& a, const ChildHit& b) -> bool {
				const VecDf& dist_a = (a.pos_intersect - pos_origin);
				const VecDf& dist_b = (b.pos_intersect - pos_origin);
				return dist_a.squaredNorm() < dist_b.squaredNorm();
			}
		);
		// Remove any duplicate child grids from the sorted list (i.e. where ray intersects
		// precisely at the interesction of two or more planes).
		child_hits.erase(std::unique(
			child_hits.begin(), child_hits.end(),
			[&pos_origin](const ChildHit& a, const ChildHit& b) -> bool {
				return a.pos_child == b.pos_child;
			}
		), child_hits.end());

		// For each candidate child, cast ray through until the zero-curve is hit.
		for (const ChildHit& child_hit : child_hits)
		{
			const VecDf& pos_hit = this->ray(
				child_hit.pos_intersect, dir,
				this->phi().child(child_hit.pos_child)
			);

			if (pos_hit != NULL_POS<FLOAT>())
				return pos_hit;
		}

		return NULL_POS<FLOAT>();
	}

#ifndef _TESTING
protected:
#endif

	/**
	 * Cast a ray to the zero layer within a given child grid.
	 *
	 * @param pos_origin
	 * @param dir
	 * @return NULL_POS if no hit, otherwise interpolated position on zero
	 * curve.
	 */
	const VecDf ray(
		VecDf pos_sample, const VecDf& dir, const typename PhiGrid::Child& child
	) const {
		using Line = Eigen::ParametrizedLine<FLOAT, D>;

		const Line line_leaf(pos_sample, dir);
		FLOAT t_leaf = 0;

		while (child.inside(pos_sample))
		{
			const INT& layer_id = this->layer_id(pos_sample);

			if (abs(layer_id) == 0)
			{
				VecDf normal = this->phi().grad(pos_sample);
				normal.normalize();

				if (normal.dot(dir) < 0)
				{
					static const UINT MAX_CONVERGE_STEPS = 100;
					UINT num_converge_steps = 0;
					FLOAT dist = 0;
					for (; num_converge_steps < MAX_CONVERGE_STEPS; num_converge_steps++)
					{
						dist = this->phi().interp(pos_sample);

						pos_sample -= normal*dist;

						if (std::abs(dist) <= TINY || normal.dot(dir) >= 0)
							break;

						normal = this->phi().grad(pos_sample);
						normal.normalize();
					}

					#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
					if (num_converge_steps == MAX_CONVERGE_STEPS)
					{
						std::cerr << "WARNING: raycast failed to get to distance < "
							<< TINY << ".\n"
							<< "dir(" << felt::format(dir) << ");"
							<< " normal(" << felt::format(normal) << ");"
							<< " sample(" << felt::format(pos_sample) << ")"
							<< " at dist " << dist;
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
	 * Cast ray to plane, get child at that point, and add to list if it contains zero-curve.
	 *
	 * @param child_hits
	 * @param line
	 * @param plane
	 * @return
	 */
	bool ray_check_add_child(
		std::vector<ChildHit>& child_hits, const Line& line, const Plane& plane
	) const {
		const VecDf& pos_intersect = (
			line.intersectionPoint(plane) + line.direction() * FLOAT(TINY)
		);

		if (!this->m_grid_phi.inside(pos_intersect))
			return false;

		const VecDi& pos_floor = floor(pos_intersect);
		const VecDi& pos_child = this->m_grid_phi.pos_child(pos_floor);

		if (this->layer(pos_child, 0).size())
		{
			child_hits.push_back(ChildHit(pos_intersect, pos_child));
		}
		return true;
	}

	/**
	 * Along a given dimension at given position, round up or down to border of next child grid.
	 *
	 * @param dim
	 * @param dir
	 * @param pos
	 * @return
	 */
	FLOAT round_to_next_child(const UINT& dim, const FLOAT& dir, const FLOAT& pos) const
	{
		// Real-valued child pos translated to [0, 2*childsize) space.
		FLOAT pos_plane_dim = (
			FLOAT(pos - m_grid_phi.offset()(dim)) / m_grid_phi.child_dims()(dim)
		);
		// Round to next child en route in [0, 2*childsize) space.
		pos_plane_dim = (dir == -1) ?
			std::floor(pos_plane_dim) : std::ceil(pos_plane_dim);
		// Scale back to phi grid in [0, 2*fullsize) space.
		pos_plane_dim *= m_grid_phi.child_dims()(dim);
		// Translate back to phi grid in [-fullsize, fullsize) space.
		pos_plane_dim += m_grid_phi.offset()(dim);

		return pos_plane_dim;
	}
};
}
#endif
