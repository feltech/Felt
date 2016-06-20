#ifndef Surface_hpp
#define Surface_hpp

#include <vector>
#include <functional>
#include <limits>
#include <set>
#include <iostream>
#include <boost/math/special_functions/round.hpp>
#include <boost/math/constants/constants.hpp>
#include <eigen3/Eigen/Dense>
#include <omp.h>

#include "Util.hpp"
#include "LookupPartitionedGrid.hpp"
#include "TrackedPartitionedGrid.hpp"
#include "PartitionedArray.hpp"

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
	/// Furthest layer from the zero-layer on the inside of the volume.
	static constexpr INT LAYER_MIN	= (INT)-L;
	/// Furthest layer from the zero-layer on the outside of the volume.
	static constexpr INT LAYER_MAX	= (INT)L;
	/// Total number of layers.
	static constexpr UINT NUM_LAYERS = 2*L+1;
	/// A tiny number used for error margin when raycasting.
	static constexpr FLOAT TINY = 0.00001f;
	/**
	 * A delta isogrid update grid with active (non-zero) grid points tracked.
	 */
	using DeltaIsoGrid = LazySingleTrackedPartitionedGrid<FLOAT, D, NUM_LAYERS>;
	/**
	 * A level set embedding isogrid grid, with active grid points (the narrow
	 * band) tracked.
	 */
	using IsoGrid = LazySingleTrackedPartitionedGrid<FLOAT, D, NUM_LAYERS>;
	/**
	 * A resizable array of D-dimensional grid positions.
	 */
	using PosArray = typename IsoGrid::PosArray;
	/**
	 * D-dimensional unsigned int vector.
	 */
	using VecDu = typename IsoGrid::VecDu;
	/**
	 * D-dimensional integer vector.
	 */
	using VecDi = typename IsoGrid::VecDi;
	/**
	 * D-dimensional float vector.
	 */
	using VecDf = typename IsoGrid::VecDf;
	/**
	 * Grid to track positions that require an update.
	 */
	using AffectedMultiLookupGrid = LazySingleLookupPartitionedGrid<D, 2*L+1>;

	/// D-dimensional hyperplane type (using Eigen library), for raycasting.
	using Plane = Eigen::Hyperplane<FLOAT, D>;
	/// D-dimensional parameterised line, for raycasting.
	using Line = Eigen::ParametrizedLine<FLOAT, D>;

	/**
	 * Storage class for flagging a point in the grid to be moved from one
	 * narrow band layer to another.
	 */
	struct StatusChange
	{
		StatusChange(
			const VecDi& pos_, const INT from_layer_, const INT to_layer_
		) : pos(pos_), from_layer(from_layer_), to_layer(to_layer_) {}

		static UINT layer_idx(const INT layer_id_)
		{
			return layer_id_ + (L + 1);
		}

		static INT layer_id(const INT layer_idx_)
		{
			return layer_idx_ - (L + 1);
		}

		VecDi pos;
		INT from_layer;
		INT to_layer;
	};

protected:
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
	 * Additional space is required at the border of the isogrid embedding for
	 * the outer layers, so this is the effective minimum position a point
	 * on the surface can occupy.
	 */
	VecDi m_pos_min;
	/**
	 * The maximum usable position in the grid for the surface (zero-layer).
	 *
	 * Additional space is required at the border of the isogrid embedding for
	 * the outer layers, so this is the effective maximum position a point
	 * on the surface can occupy.
	 */
	VecDi m_pos_max;

	/**
	 * Grid for preventing duplicates when doing neighbourhood queries.
	 */
	AffectedMultiLookupGrid 	m_grid_affected;

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


public:

	/**
	 * Explicitly defined default constructor.
	 */
	Surface () = default;

	/**
	 * Construct a level set embedding of size size.
	 *
	 * All points will be marked as outside the surface (i.e. no surface).
	 *
	 * @param size_ size of the isogrid.
	 * @param size_partition_ size of each spatial partition of the isogrid.
	 */
	Surface (const VecDu& size_, const VecDu& size_partition_ = VecDu::Constant(8)) :
		m_grid_isogrid(),
		m_grid_status_change(),
		m_grid_delta(),
		m_grid_affected()
	{
		this->init(size_, size_partition_);
	}

	/**
	 * Initialise a Surface (e.g. after default-construction).
	 *
	 * @param size_ size of the isogrid.
	 * @param size_partition_ size of each spatial partition of the isogrid.
	 */
	void init (const VecDu& size_, const VecDu& size_partition_)
	{
		this->size(size_, size_partition_);
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
	constexpr const felt::VecDT<T, D> NULL_POS() const
	{
		return felt::VecDT<T, D>::Constant(std::numeric_limits<T>::max());
	};

	/**
	 * Shorthand for accessing the isogrid (level set) grid at a given position.
	 *
	 * @param pos
	 * @return
	 */
	FLOAT& operator() (const VecDi& pos_)
	{
		return this->isogrid()(pos_);
	}

	/**
	 * Shorthand for accessing the isogrid (level set) grid (const version).
	 *
	 * @param pos
	 * @return
	 */
	const FLOAT operator() (const VecDi& pos) const
	{
		return this->isogrid()(pos);
	}

	/**
	 * Initialise level set embedding with given dimensions.
	 *
	 * Initialises the various lookup grids and (indirectly) the level set
	 * sparse field layers. Also calculates and stores the spatial limits
	 * of the grid accounting for the narrow band space required.
	 *
	 * @param usize
	 */
	void size (const VecDu& usize_, const VecDu& size_partition)
	{
		const VecDi isize = usize_.template cast<INT>();
		const VecDi offset = -1 * isize / 2;

		// Configure isogrid embedding, initialising to all outside values.
		m_grid_isogrid.init(usize_, offset, LAYER_MAX+1, size_partition);
		// Configure delta isogrid embedding, initialising to zero delta.
		m_grid_delta.init(usize_, offset, 0, size_partition);
		// Configure status change partitioned lists.
		m_grid_status_change.init(usize_, offset, size_partition);
		// Configure de-dupe grid for neighbourhood queries.
		m_grid_affected.init(usize_, offset, size_partition);

		// Store min and max usable positions in isogrid embedding.
		this->pos_min(
			VecDi::Constant(L + 1) + m_grid_isogrid.offset()
		);
		this->pos_max(
			(isize - VecDi::Constant(L + 1))
			+ m_grid_isogrid.offset() - VecDi::Constant(1)
		);
	}

	/**
	 * Get minimum usable position in isogrid grid.
	 *
	 * @return minimum position a zero layer point could be found.
	 */
	const VecDi& pos_min () const
	{
		return m_pos_min;
	}

	/**
	 * Get maximum usable position in isogrid grid.
	 *
	 * @return maximum position a zero layer point could be found.
	 */
	const VecDi& pos_max () const
	{
		return m_pos_max;
	}


	/**
	 * Set minimum usable position in isogrid grid.
	 *
	 * @param pos
	 */
	void pos_min (const VecDi& pos)
	{
		m_pos_min = pos;
	}

	/**
	 * Set maximum usable position in isogrid grid.
	 *
	 * @param pos
	 */
	void pos_max (const VecDi& pos)
	{
		m_pos_max = pos;
	}


	/**
	 * Get reference to isogrid grid.
	 *
	 * @return
	 */
	IsoGrid& isogrid ()
	{
		return m_grid_isogrid;
	}

	/**
	 * Get reference to isogrid grid.
	 *
	 * @return
	 */
	const IsoGrid& isogrid () const
	{
		return m_grid_isogrid;
	}

	/**
	 * Shorthand for accessing isogrid grid at given position (const
	 * version).
	 *
	 * @param pos
	 * @return
	 */
	const FLOAT isogrid (const VecDi& pos) const
	{
		return m_grid_isogrid(pos);
	}

	/**
	 * Shorthand for accessing isogrid grid at given position.
	 *
	 * @param pos
	 * @return
	 */
	FLOAT& isogrid (const VecDi& pos)
	{
		return m_grid_isogrid(pos);
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

	const AffectedMultiLookupGrid& affected ()
	{
		return m_grid_affected;
	}

	/**
	 * Update isogrid grid point at pos by val.
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
	void isogrid (const VecDi& pos_, const FLOAT val_, const INT layer_id_ = 0)
	{
		const INT newlayer_id = this->layer_id(val_);

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		
		if (
			std::abs(this->layer_id(val_)) != std::abs(layer_id_)
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) + 1
			&& std::abs(this->layer_id(val_)) != std::abs(layer_id_) - 1
		) {
			std::stringstream strs;
			strs << "Iso update value out of bounds at layer "
				<< layer_id_ << " " << felt::format(pos_) << ": " << val_;
			std::string str = strs.str();
			throw std::domain_error(str);
		}
		
		#endif

		m_grid_isogrid(pos_) = val_;

		this->status_change(pos_, layer_id_, newlayer_id);
	}


	void add_neighs(const VecDi& pos, const INT side)
	{
		// Get neighbouring points.
		using Child = typename IsoGrid::Child;
		PosArray neighs;
		m_grid_isogrid.neighs(pos, neighs);

		const INT to_layer_id = side*INT(L);

		// Get a mutex lock on affected spatial partitions, since neighbouring points might spill
		// over into another thread's partition.
		for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
		{
			const VecDi& pos_neigh = neighs[neighIdx];
			const VecDi pos_child = m_grid_isogrid.pos_child(pos_neigh);

			m_grid_isogrid.add_child(pos_child);
			Child& child = m_grid_isogrid.children().get(pos_child);
			
			const INT from_layer_id = this->layer_id(pos_neigh);

			// Only add if neighbouring point is not already within the narrow band.
			if (this->inside_band(from_layer_id))
				continue;

			// Get distance of this new point to the zero layer.
			const FLOAT dist_neigh = this->distance(pos_neigh, side);

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
			m_grid_status_change.add_safe(
				pos_neigh, StatusChange(pos_neigh, from_layer_id, to_layer_id)
			);
			// Set distance in isogrid grid.
			child(pos_neigh) = dist_neigh;
		} // End for neighbouring points.
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
	bool status_change (const VecDi& pos_, const INT layer_id_from_, const INT layer_id_to_)
	{
		if (layer_id_from_ == layer_id_to_)
			return false;

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		{
			m_grid_isogrid.assert_pos_bounds(pos_, "status_change: ");

			if (
				std::abs(layer_id_from_) != std::abs(layer_id_to_) + 1
				&& std::abs(layer_id_from_) != std::abs(layer_id_to_) - 1
			) {

				std::stringstream strs;
				strs << "Bad status_change: " << felt::format(pos_) << " from layer " <<
					layer_id_from_ << " to layer " << layer_id_to_;
				std::string str = strs.str();
				throw std::domain_error(str);
			}
		}
		#endif

		m_grid_status_change.add_safe(pos_, StatusChange(pos_, layer_id_from_, layer_id_to_));

		return true;
	}

	/**
	 * Loop through the status change lists moving the referenced points
	 * from one layer to another.
	 */
	void flush_status_change ()
	{
		for (const VecDi& pos_child : m_grid_status_change.children().list())
			for (const StatusChange& status : m_grid_status_change.children().get(pos_child))
			{
				const INT layer_id_to = status.to_layer;
				const INT layer_id_from = status.from_layer;
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
	 * Get reference to delta isogrid grid.
	 *
	 * @return
	 */
	DeltaIsoGrid& delta ()
	{
		return m_grid_delta;
	}

	/**
	 * Get reference to delta isogrid grid (const version).
	 *
	 * @return
	 */
	const DeltaIsoGrid& delta () const
	{
		return m_grid_delta;
	}

	/**
	 * Shorthand for access to the delta isogrid grid at given position.
	 *
	 * @param pos
	 * @return
	 */
	FLOAT& delta (const VecDi& pos)
	{
		return m_grid_delta(pos);
	}

	/**
	 * Shorthand for access to the delta isogrid grid at given position (const
	 * version)
	 *
	 * @param pos
	 * @return
	 */
	const FLOAT delta (const VecDi& pos) const
	{
		return m_grid_delta(pos);
	}

	/**
	 * Update delta isogrid grid and append point to change list for given thread.
	 *
	 * @snippet test_Surface.cpp Delta isogrid clamping
	 *
	 * @snippet test_Surface.cpp Simple delta isogrid update
	 *
	 * @param pos
	 * @param val
	 */
	void delta (const VecDi& pos_, FLOAT val_, const INT layer_id_ = 0)
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
		
		this->delta().add(pos_, val_, this->layer_idx(layer_id_));
	}

	/**
	 * Clamp delta isogrid value such that it doesn't breach the grid or cause
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
				// Get isogrid at this point.
				const FLOAT fisogrid = this->isogrid(pos);
				// Max value that will not be rounded and thus trigger
				// a layer_move.
				const FLOAT val_max = -0.5 +
					(std::numeric_limits<FLOAT>::epsilon() * 2);
				// Clamp the value of delta isogrid.
				val = std::max(val_max - fisogrid, val);
				break;
			}
		return val;
	}

	/**
	 * Update delta isogrid grid surrounding a given point with amount determined
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
	FLOAT delta_gauss (
		const VecDf& pos_centre, const FLOAT val, const FLOAT stddev
	) {
		const SingleLookupGrid<D, NUM_LAYERS>& lookup = this->walk_band<Distance>(
			felt::round(pos_centre)
		);
		return this->delta_gauss(lookup.list(this->layer_idx(0)), pos_centre, val, stddev);

//		PosArray list;
//		const INT w = Distance;
//		const VecDi offset = VecDi::Constant(-w/2);
//		const VecDu size = VecDu::Constant(w);
//		const VecDi pos_centre_rounded = round(pos_centre);
//		for (UINT i = 0; i < std::pow(w, D); i++)
//		{
//			const VecDi& pos_offset = Grid<UINT, D>::index(i, size, offset);
//			const VecDi& pos_neigh = pos_offset + pos_centre_rounded;
//			if (this->layer_id(pos_neigh) == 0)
//				list.push_back(pos_neigh);
//		}
//
//		return this->delta_gauss(list, pos_centre, val, stddev);
	}

//	std::set<std::mutex*> lock_children(const PosArray& list_pos_leafs_)
//	{
//		using MutexSet = std::set<std::mutex*>;
//		using MutexIter = boost::indirect_iterator<MutexSet::iterator>;
//
//		MutexSet mutexes;
//		for (const VecDi& pos : list_pos_leafs_)
//			mutexes.insert(&m_grid_delta.children().get(m_grid_delta.pos_child(pos)).lookup().mutex());
//
//		MutexIter first(mutexes.begin()), last(mutexes.end());
//		boost::lock(first, last);
//
//		return mutexes;
//	}

	/**
	 * Update delta isogrid grid at given points with amount determined by Gaussian
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
	FLOAT delta_gauss (
		const typename DeltaIsoGrid::PosArray& list, const VecDf& pos_centre,
		const FLOAT val, const FLOAT stddev
	) {
		constexpr FLOAT sqrt2piDinv = 1.0f/sqrt(pow(2*boost::math::constants::pi<FLOAT>(), D));

		Eigen::VectorXf weights(list.size());

		// Calculate Gaussian weights.
		for (UINT idx = 0; idx < list.size(); idx++)
		{
			const VecDf& pos = list[idx].template cast<FLOAT>();
			if (this->delta(list[idx]) != 0)
			{
				weights(idx) = 0;
				continue;
			}

			const FLOAT dist_sq = (pos - pos_centre).squaredNorm();
			const FLOAT weight = sqrt2piDinv * exp(-0.5f * dist_sq);
			weights(idx) = weight;
		}

		std::mutex* pmutex_current = NULL;
		const FLOAT weights_sum = weights.sum();

		if (weights_sum > 0)
		{
			// Normalise weights
			weights *= val / weights_sum;

			VecDi pos_child_current = VecDi::Constant(std::numeric_limits<INT>::max());


			for (UINT idx = 0; idx < list.size(); idx++)
			{
				const VecDi& pos = list[idx];
				const VecDi& pos_child = m_grid_delta.pos_child(pos);

				if (pos_child_current != pos_child)
				{
					if (pmutex_current)
					{
						pmutex_current->unlock();
					}
					pos_child_current = pos_child;
					pmutex_current = &m_grid_delta.children().get(pos_child_current).lookup().mutex();
					pmutex_current->lock();
				}

				const FLOAT amount = this->clamp(pos, weights(idx) + m_grid_delta(pos));

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (std::abs(amount) > 1)
				{
					std::stringstream strs;
					strs << "Delta isogrid value out of bounds at " << felt::format(pos)
						<< ": " << amount;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				if (this->layer_id(pos) != 0)
				{
					std::stringstream strs;
					strs << "Attempting to Gaussian update non-zero layer " << felt::format(pos)
						<< ": " << m_grid_isogrid(pos) << " += " << amount;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				
				#endif
				
				weights(idx) -= (weights(idx) - amount);

				this->delta().add(pos, amount, this->layer_idx(0));
			} // End for list of leafs.
		} // End if weights > 0

		if (pmutex_current)
		{
			pmutex_current->unlock();
		}

		return val - weights.sum();
	}

	/**
	 * Update delta isogrid grid surrounding a zero layer point found via a
	 * raycast by amount spread using Gaussian distribution.
	 *
	 * @param pos_origin
	 * @param dir
	 * @param dist
	 * @param amount
	 * @param stddev
	 */
	template <UINT Distance>
	FLOAT delta_gauss (
		const VecDf& pos_origin, const VecDf& dir,
		const FLOAT val, const float stddev
	) {
		const VecDf& pos_hit = this->ray(pos_origin, dir);

		if (pos_hit == NULL_POS<FLOAT>())
		{
//			std::cout << "MISS" << std::endl;
			return val;
		}
		return this->delta_gauss<Distance>(pos_hit, val, stddev);
	}

	/**
	 * Get reference to the active partitions of a given layer of the narrow
	 * band.
	 *
	 * @param layer_id
	 * @return
	 */
	PosArray& parts (const INT id = 0)
	{
		return m_grid_isogrid.children().list(id+L);
	}

	/**
	 * Get reference to a single layer of the narrow band at a given
	 * spatial partition.
	 *
	 * @param pos_child
	 * @param id
	 * @return
	 */
	PosArray& layer (const VecDi& pos_child, const INT id = 0)
	{
		return m_grid_isogrid.children().get(pos_child).list(id+L);
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
		const VecDi& pos_child, const INT id = 0
	) const
	{
		return m_grid_isogrid.children().get(pos_child).list(id+L);
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
	const LeafsContainer<IsoGrid> layer(const UINT layer_id) const
	{
		return m_grid_isogrid.leafs(this->layer_idx(layer_id));
	}

	/**
	 * Append position to a layer of the narrow band.
	 *
	 * @param id
	 * @param pos
	 */
	void layer_add (const VecDi& pos, const INT layer_id)
	{
		// Do nothing if position is outside narrow band.
		if (!this->inside_band(layer_id))
			return;
		m_grid_isogrid.add(pos, this->layer_idx(layer_id));
	}

	/**
	 * Append position to a layer of the narrow band.
	 *
	 * Skip the lookup in the isogrid grid by assuming passed val is the isogrid
	 * grid value of this point.
	 *
	 * @param pos
	 * @param val
	 */
	void layer_add (const FLOAT val, const VecDi& pos)
	{
		const INT layer_id = this->layer_id(val);
		// Do nothing if position is outside narrow band.
		if (!this->inside_band(layer_id))
			return;
		m_grid_isogrid.add(pos, val, this->layer_idx(layer_id));
	}

	/**
	 * Remove a position from a given layer of the narrow band.
	 *
	 * Does not modify underlying grid value, just the layer list.
	 *
	 * @param pos
	 * @param layer_id
	 */
	void layer_remove(const VecDi& pos, const INT layer_id)
	{
		// Do nothing if position is outside narrow band.
		if (!this->inside_band(layer_id))
			return;

		m_grid_isogrid.remove(pos, this->layer_idx(layer_id), layer_id + sgn(layer_id));
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
		const VecDi& pos, const INT fromlayer_id, const INT tolayer_id
	) {
		this->layer_remove(pos, fromlayer_id);
		this->layer_add(pos, tolayer_id);
	}

	/**
	 * Get narrow band layer id of location in isogrid grid.
	 * @param pos
	 * @return
	 */
	template <class PosType>
	INT layer_id(const PosType& pos) const
	{
		return this->layer_id(this->isogrid()(pos));
	}

	/**
	 * Get narrow band layer id of value.
	 * @param val
	 * @return
	 */
	INT layer_id(const FLOAT val) const
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
	UINT layer_idx(const INT id) const
	{
		return id+L;
	}

	/**
	 * Create a single singularity seed point in the isogrid grid.
	 * NOTE: does not handle overwriting of points currently already on the
	 * surface/in the volume.
	 *
	 * @param pos_centre
	 */
	void seed (const VecDi& pos_centre)
	{
		IsoGrid& isogrid = this->isogrid();
		const UINT dims = isogrid.size().size();

		// Width of seed.
		const VecDi vec_width = VecDi::Constant(L);

		// Min and max positions affected by placing seed point.
		const VecDi pos_min = pos_centre - vec_width;
		const VecDi pos_max = pos_centre + vec_width;

		// Get vector size of window formed by pos_min and pos_max.
		const VecDu window_size = (
			pos_max - pos_min + VecDi::Constant(1)
		).template cast<UINT>(); //+1 for zero coord.

		// Calculate number of grid points to be cycled through within window.
		UINT total_size = 1;
		for (INT axis = 0; axis < dims; axis++)
			total_size *= window_size(axis);

		// Cycle through each point in window.
		for (UINT idx_pos = 0; idx_pos <= total_size; idx_pos++)
		{
			// Calculate vector position from integer index, using Felt::Grid utility function,
			// index().
			VecDi pos = IsoGrid::index(idx_pos, window_size);
			// Translate position into isogrid grid space.
			pos += pos_min;
			// Calculate vector distance from this position to seed centre.
			const VecDi vec_dist = pos - pos_centre;
			// Sum of absolute distance along each axis == city-block distance.
			FLOAT f_dist = (FLOAT)vec_dist.template lpNorm<1>();
			// Check distance indicates that this point is within the narrow band.
			if ((UINT)std::abs(this->layer_id(f_dist)) <= L)
			{
				// Append point to a narrow band layer (if applicable).
				this->layer_add(f_dist, pos);
			}
		}
	}

	/**
	 * Get neighbouring position in isogrid grid that is closest to
	 * zero-curve.
	 *
	 * @param pos
	 * @param side
	 * @return
	 */
	VecDi next_closest (const VecDi& pos, const FLOAT side) const
	{
		// Trivially return if this is already a zero-layer point.
		if (this->layer_id(pos) == 0)
			return pos;

		const IsoGrid& isogrid = this->isogrid();

		// Get all neighbours of this point.
		PosArray neighs;
		isogrid.neighs(pos, neighs);

		VecDi pos_nearest = VecDi(pos);
		FLOAT val_nearest = isogrid(pos)*side;
		// Cycle neighbours finding one that is closest to zero-layer.
		for (UINT neighIdx = 0; neighIdx < neighs.size(); neighIdx++)
		{
			const VecDi pos_neigh = neighs[neighIdx];
			const FLOAT val_neigh = isogrid(pos_neigh);
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
//			const VecDf vec_dir = isogrid.grad(pos) * dir;
//			const UINT axis_best = ublas::index_norm_inf(vec_dir);
//			VecDi pos_nearest = VecDi(pos);
//			const FLOAT val_best = vec_dir(axis_best);
//			if (fabs(val_best) > 0)
//				pos_nearest(axis_best) += sgn(val_best);

		return pos_nearest;
	}

	/**
	 * Get neighbouring position in isogrid grid that is closest to the
	 * zero-curve.
	 *
	 * @param pos
	 * @return
	 */
	VecDi next_closest (const VecDi& pos) const
	{
		const IsoGrid& isogrid = this->isogrid();
		const FLOAT val_centre = isogrid(pos);
		// Direction multiplier for gradient toward zero-curve.
		const FLOAT side = sgn(val_centre);

		return this->next_closest(pos, side);
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
		m_grid_status_change.reset();
		
		#if false
		
		for (const VecDi& pos_child : m_grid_delta.children())
			for (const VecDi& pos : m_grid_delta.children().get(pos_child))
			{
				if (m_grid_delta(pos) != 0)
				{
					throw std::domain_error(std::string("Delta isogrid not reset!"));
				}
			}
		#endif
	}


	/**
	 * Apply delta isogrid to isogrid along the zero layer.
	 */
	void update_zero_layer ()
	{
		const typename DeltaIsoGrid::ChildrenGrid&
		children = this->delta().children();

		const UINT layer_idx = this->layer_idx(0);

		#pragma omp parallel for
		for (
			UINT idx_child = 0; idx_child < children.list(layer_idx).size();
			idx_child++
		) {
			const VecDi& pos_child = children.list(layer_idx)[idx_child];

			m_grid_isogrid.add_child(pos_child);

			for (const VecDi& pos : children(pos_child).list(layer_idx))
			{
				const FLOAT fisogrid = this->isogrid(pos);
				const FLOAT fdelta = this->delta(pos);
				const FLOAT fval = fisogrid + fdelta;

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (this->layer_id(fisogrid) != 0)
				{
					std::stringstream strs;
					strs << "Zero layer updated attempted at non-zero layer point " <<
						felt::format(pos) << ": " << fisogrid << " + " << fdelta << " = " << fval;
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				if (std::abs(this->layer_id(fval)) > 1)
				{
					std::stringstream strs;
					strs << "Zero layer isogrid value out of bounds at " << felt::format(pos)
						<< ": " << fisogrid << " + " << fdelta << " = " << fval;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif

				this->isogrid(pos, fval);
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

		this->update_distance(m_grid_isogrid);

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

		// Update the zero layer, applying delta to isogrid.
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
	void update_distance(const INT layer_id_, const INT side_, const GridType& lookup_)
	{
		using DeltaIsoChild = typename DeltaIsoGrid::Child;
		using IsoChild = typename IsoGrid::Child;
		const UINT layer_idx = this->layer_idx(layer_id_);

		#pragma omp parallel for
		for (UINT pos_idx = 0; pos_idx < lookup_.children().list(layer_idx).size(); pos_idx++)
		{
			const VecDi& pos_child = lookup_.children().list(layer_idx)[pos_idx];

			m_grid_delta.add_child(pos_child, layer_idx);

			IsoChild& grid_isogrid_child = m_grid_isogrid.children().get(pos_child);
			DeltaIsoChild& grid_delta_child = m_grid_delta.children().get(pos_child);

			const PosArray& apos_leafs = lookup_.children().get(pos_child).list(layer_idx);

			// Calculate distance of every point in this layer to the zero
			// layer, and store in delta isogrid grid.
			// Delta isogrid grid is used to allow for asynchronous updates, that
			// is, to prevent neighbouring points affecting the distance
			// transform.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance from this position to zero layer.
				const FLOAT dist = this->distance(pos, side_);

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)

				if (
					this->layer_id(dist) != layer_id_ &&
					this->layer_id(dist) != layer_id_ + 1 &&
					this->layer_id(dist) != layer_id_ - 1
				) {
					std::stringstream strs;
					strs << "Outer layer distance update value out of bounds at layer "
						<< layer_id_ << " " << felt::format(pos) << ": " << dist;
					std::string str = strs.str();
					throw std::domain_error(str);
				}

				#endif


				// Update delta isogrid grid.
				grid_delta_child.add(pos, dist, layer_idx);
			}
		}
		
		#pragma omp parallel for
		for (UINT pos_idx = 0; pos_idx < lookup_.children().list(layer_idx).size(); pos_idx++)
		{
			const VecDi& pos_child = lookup_.children().list(layer_idx)[pos_idx];

			IsoChild& grid_isogrid_child = m_grid_isogrid.children().get(pos_child);
			DeltaIsoChild& grid_delta_child = m_grid_delta.children().get(pos_child);

			const PosArray& apos_leafs = lookup_.children().get(pos_child).list(layer_idx);

			// Update distance in isogrid from delta isogrid and append any points that
			// move out of their layer to a status change list.
			for (const VecDi& pos : apos_leafs)
			{
				// Distance calculated above.
				const FLOAT dist = grid_delta_child(pos);

				grid_isogrid_child(pos) = dist;

				const INT newlayer_id = this->layer_id(dist);

				if (this->status_change(pos, layer_id_, newlayer_id))
				{
					// If outside point moving inward, must create new outside points.
					if (std::abs(layer_id_) == L && std::abs(newlayer_id) == L-1)
					{
						// Get which side of the zero-layer this point lies on.
						const INT side = sgn(newlayer_id);
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
	FLOAT distance (const VecDi& pos_, const FLOAT side_) const
	{
		const IsoGrid& isogrid = this->isogrid();
		// Get neighbouring point that is next closest to the zero-layer.
		const VecDi pos_closest = this->next_closest(pos_, side_);
		const FLOAT val_closest = isogrid(pos_closest);
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
		SortablePos(const VecDi& pos_, const UINT distance_)
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
	 * @return SingleLookupGrid with tracking lists for visited points, one list
	 * for each layer.
	 */
	template <UINT Distance>
	SingleLookupGrid<D, NUM_LAYERS>& walk_band (const VecDi& pos_)
	{
		using MultiLookup = SingleLookupGrid<D, 2*L+1>;
		using PosArray = typename MultiLookup::PosArray;

		// Box size is: (Distance * 2 directions) + 1 central.
		static MultiLookup lookup(VecDu::Constant(Distance * 2 + 1));
		lookup.reset_all();
		lookup.offset(pos_ - VecDi::Constant(Distance));

		const INT layer_id = this->layer_id(pos_);
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

					this->isogrid().neighs(
						pos_centre,
						[this](const VecDi& pos_neigh) {
							// Calculate layer of this neighbouring
							// point from the isogrid grid.
							const INT layer_id = this->layer_id(pos_neigh);
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
	 * Find all outer layer points who's distance transform is affected by modified zero-layer
	 * points.
	 *
	 * @snippet test_Surface.cpp Calculate affected outer layers for localised narrow band updates
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
		const UINT layer_idxZero = this->layer_idx(0);

		// Loop over delta isogrid modified zero-layer points adding to
		// tracking grid.

		// Loop spatial partitions of delta for zero-layer.
		for (
			const VecDi& pos_child
			: this->delta().children().list(layer_idxZero))
		{
			// Loop leaf grid nodes with spatial partition
			for (
				const VecDi& pos_leaf
				: this->delta().children().get(pos_child).list(layer_idxZero)
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
				const UINT layer_idx = this->layer_idx(layer_id);
				// Get number of spatial partitions for this layer.
				const UINT num_childs = (
					m_grid_affected.children().list(layer_idx).size()
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
					const VecDi& pos_child = m_grid_affected.children().list(layer_idx)[idx_child];
					// Copy number of active grid nodes for this partition
					// into relevant index in the list.
					aidx_last_neigh[layer_idx][idx_child] = (
						m_grid_affected.children().get(pos_child).list(layer_idx).size()
					);
				}
			}

			// Loop each layer finding the affected outer layer points
			// for each partition using the start and end points cached
			// above.

			for (INT layer_id = LAYER_MIN; layer_id <= LAYER_MAX; layer_id++)
			{
				const UINT layer_idx = this->layer_idx(layer_id);

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
						m_grid_affected.children().list(layer_idx)[idx_child]
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
							m_grid_affected.children().get(pos_child).list(layer_idx)
						);
						// This leaf grid nodes is the centre to search
						// about.
						const VecDi& pos_centre = apos_neigh[idx_neigh];

						// Use utility method from Grid to get neighbouring
						// grid nodes and call a lambda to add the point
						// to the appropriate tracking list.

						this->isogrid().neighs(
							pos_centre,
							[this](const VecDi& pos_neigh) {
								// Calculate layer of this neighbouring
								// point from the isogrid grid.
								const INT layer_id = this->layer_id(pos_neigh);
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
	 * Perform a a full (parallelised) update of the narrow band.
	 *
	 * Lambda function passed will be given the position to process and
	 * a reference to the isogrid grid, and is expected to return delta isogrid to
	 * apply.
	 *
	 * @param fn_ (pos, isogrid) -> float
	 */
	void update(std::function<FLOAT(const VecDi&, const IsoGrid&)> fn_)
	{
		this->update_start();
		#pragma omp parallel for
		for (UINT part_idx = 0; part_idx < parts().size(); part_idx++)
		{
			const VecDi& pos_part = this->parts()[part_idx];
			for (const VecDi& pos : this->layer(pos_part))
				this->delta(pos, fn_(pos, m_grid_isogrid));
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

		// If ray is cast from within isogrid grid, first check child grid containing origin point.
		if (m_grid_isogrid.inside(pos_origin))
		{
			const VecDf& pos_hit = ray(
				pos_origin, dir,
				m_grid_isogrid.children().get(
					m_grid_isogrid.pos_child(pos_origin.template cast<INT>())
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

			// Cast ray to plane and add any child grids hit on the way to tracking list.
			if (!ray_check_add_child(child_hits, line, Plane(normal, pos_plane(dim) * dir_dim)))
				continue;

			// Round up/down to next child, in case we started at inexact modulo of child grid size
			// (i.e. when isogrid grid size is not integer multiple of child grid size).
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
			const FLOAT child_size_dim = m_grid_isogrid.child_size()(dim);
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
			[](const ChildHit& a, const ChildHit& b) -> bool {
				return a.pos_child == b.pos_child;
			}
		), child_hits.end());

		// For each candidate child, cast ray through until the zero-curve is hit.
		for (const ChildHit& child_hit : child_hits)
		{
			const VecDf& pos_hit = this->ray(
				child_hit.pos_intersect, dir,
				this->isogrid().children().get(child_hit.pos_child)
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
		VecDf pos_sample, const VecDf& dir, const typename IsoGrid::Child& child
	) const {
		using Line = Eigen::ParametrizedLine<FLOAT, D>;

		const Line line_leaf(pos_sample, dir);
		FLOAT t_leaf = 0;

		while (child.inside(pos_sample))
		{
			const INT layer_id = this->layer_id(pos_sample);

			if (abs(layer_id) == 0)
			{
				VecDf normal = this->isogrid().grad(pos_sample);

				#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
				if (std::isnan(normal.normalized()[0]))
				{
					std::stringstream strs;
					strs << "ERROR: raycast isosurface gradient normal is NaN " <<
						"when normalising " << felt::format(normal);
					std::string str = strs.str();
					throw std::domain_error(str);
				}
				#endif

				normal.normalize();

				if (normal.dot(dir) < 0)
				{
					static const UINT MAX_CONVERGE_STEPS = 100;
					UINT num_converge_steps = 0;
					FLOAT dist = 0;
					for (; num_converge_steps < MAX_CONVERGE_STEPS; num_converge_steps++)
					{
						dist = this->isogrid().interp(pos_sample);

						pos_sample -= normal*dist;

						if (std::abs(dist) <= TINY || normal.dot(dir) >= 0)
							break;

						normal = this->isogrid().grad(pos_sample);
						normal.normalize();
					}

					#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
					if (num_converge_steps == MAX_CONVERGE_STEPS)
					{
//						std::cerr << "WARNING: raycast failed to get to distance < "
//							<< TINY << ".\n"
//							<< "dir(" << felt::format(dir) << ");"
//							<< " normal(" << felt::format(normal) << ");"
//							<< " sample(" << felt::format(pos_sample) << ")"
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

		if (!this->m_grid_isogrid.inside(pos_intersect))
			return false;

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
	 * Along a given dimension at given position, round up or down to border of next child grid.
	 *
	 * @param dim
	 * @param dir
	 * @param pos
	 * @return
	 */
	FLOAT round_to_next_child(const UINT dim, const FLOAT dir, const FLOAT pos) const
	{
		// Real-valued child pos translated to [0, 2*childsize) space.
		FLOAT pos_plane_dim = (
			FLOAT(pos - m_grid_isogrid.offset()(dim)) / m_grid_isogrid.child_size()(dim)
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
};
}
#endif
