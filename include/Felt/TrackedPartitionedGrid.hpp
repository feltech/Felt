#ifndef INCLUDE_FELT_SINGLETRACKEDPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_SINGLETRACKEDPARTITIONEDGRID_HPP_

#include "TrackedGrid.hpp"
#include "TrackingPartitionedGridBase.hpp"


namespace felt
{

/**
 * SingleTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
class TrackedPartitionedGrid
	: public TrackingPartitionedGridBase < TrackedPartitionedGrid<T, D, N> >
{
public:
	/// This class.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// Base class
	using Base = TrackingPartitionedGridBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::ChildrenGrid;
	using LeafType = typename GridTraits<ThisType>::LeafType;

	using Base::TrackingPartitionedGridBase;
	using Base::add;
	using Base::reset;

	/**
	 * Set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add(const VecDi& pos_, const LeafType& val_, const UINT list_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		this->add_child(pos_child, list_idx_);
		return this->m_grid_children.get(pos_child).add(pos_, val_, list_idx_);
	}

	/**
	 * Thread-safely set value in grid at given position and add position to lookup grid.
	 *
	 * Will set value regardless whether lookup grid already set for given
	 * position + tracking list.
	 *
	 * Descend to relevant child grid to add to their tracking structure.
	 *
	 * @param pos_ position in grid.
	 * @param val_ value to set.
	 * @param list_idx_ tracking list id.
	 * @return true if grid node set in child lookup grid and position added
	 * to tracking list, false if child grid node was already set so
	 * position already in a list.
	 */
	bool add_safe(const VecDi& pos_, const LeafType& val_, const UINT list_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		this->add_child(pos_child, list_idx_);
		Child& child = this->m_grid_children.get(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		return child.add(pos_, val_, list_idx_);
	}

	/**
	 * Set every active grid node (i.e. those referenced by lookup grid)
	 * to given value and reset the lookup grid.
	 *
	 * MultiLookup grid will then be full of NULL indices for given tracking list
	 * and the relevant tracking list(s) will be empty.
	 *
	 * @param val_ value to set in main grid.
	 * @param list_idx_ tracking list id to cycle over and clear.
	 */
	void reset(const LeafType& val_, const UINT list_idx_ = 0)
	{
		for (const VecDi& pos_child : this->m_grid_children.list(list_idx_))
			this->m_grid_children(pos_child).reset(val_, list_idx_);

		Base::Base::reset(list_idx_);
	}

	/**
	 * Remove a leaf position from relevant child tracking structure and
	 * remove child from tracking list if child's list is now empty.
	 *
	 * @param pos_ leaf position to remove.
	 * @param list_idx_ tracking list id.
	 */
	void remove(const VecDi& pos_, const UINT list_idx_, const LeafType& background_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->m_grid_children.get(pos_child);
		child.remove(pos_, list_idx_);
		if (child.list(list_idx_).size() == 0)
			remove_child(pos_child, list_idx_, background_);
	}

	/**
	 * Move a tracked point from one tracking list to another.
	 *
	 * @param pos_ leaf position to remove.
	 * @param list_idx_ tracking list id.
	 */
	void move(
		const VecDi& pos_, const UINT list_idx_from_, const UINT list_idx_to_
	) {
		const VecDi& pos_child = this->pos_child(pos_);
		Child& child = this->m_grid_children.get(pos_child);
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		{
			if (!this->is_child_active(pos_child))
			{

				std::stringstream strs;
				strs << "Attempting to move lists within an inactive child: " <<
					felt::format(pos_) << " from list " << list_idx_from_ << " to list " <<
					list_idx_to_ << " in partition " << felt::format(pos_child);
				std::string str = strs.str();
				throw std::domain_error(str);
			}

		}
		#endif
		child.remove(pos_, list_idx_from_);
		child.add(pos_, list_idx_to_);

		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		this->m_grid_children.add(pos_child, list_idx_to_);
		if (child.list(list_idx_from_).size() == 0)
			this->m_grid_children.remove(pos_child, list_idx_from_);
	}

	/**
	 * Remove a spatial partition from children grid's tracking subgrid.
	 *
	 * Uses mutex for thread safety. Deactivates the child grid if not being tracked by any list.
	 *
	 * @param pos_ position of spatial partition to stop tracking.
	 * @param list_idx_ index of tracking list used to track position.
	 */
	void remove_child(
		const VecDi& pos_child_, const UINT list_idx_, const LeafType& background_
	) {
		std::lock_guard<std::mutex> lock(this->m_mutex_update_branch);
		if (!this->m_grid_children.is_active(pos_child_, list_idx_))
			return;
		this->m_grid_children.remove(pos_child_, list_idx_);
		if (!this->is_child_active(pos_child_))
		{
			Child& child = this->m_grid_children.get(pos_child_);
			child.background() = background_;
			this->m_grid_children.get(pos_child_).deactivate();
		}
	}

	/**
	 * Check if spatial partition is tracked.
	 *
	 * @param pos_child_ spartial partition location.
	 * @return true if spatial partition has active tracking lists, false otherwise.
	 */
	bool is_child_active(const VecDi& pos_child_) const
	{
		using NULL_IDX_TYPE = typename ChildrenGrid::Lookup::Traits::NULL_IDX_TYPE;

		const NULL_IDX_TYPE idxs = this->children().lookup().get(pos_child_);

		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		{
			if (idxs == ChildrenGrid::Lookup::Traits::NULL_IDX_DATA)
				for (UINT list_idx = 0; list_idx < Base::NUM_LISTS; list_idx++)
				{
					if (this->m_grid_children.get(pos_child_).list(list_idx).size())
					{
						std::stringstream strs;
						strs << "Children lookup is NULL_IDX, but child still contains active"
							" lists: ";
						for (UINT list_idx = 0; list_idx < Base::NUM_LISTS; list_idx++)
						{
							strs << "[" << list_idx << "] = " <<
								this->m_grid_children.get(pos_child_).list(list_idx).size() << "; ";
						}
						std::string str = strs.str();
						throw std::domain_error(str);
					}
				}

		}
		#endif

		return idxs != ChildrenGrid::Lookup::Traits::NULL_IDX_DATA;
	}

private:
	/**
	 * Ensure base class implementation cannot be called, since we must set a background value.
	 */
	void remove(const VecDi& pos_, const UINT list_idx_ = 0)
	{}
 };


/**
 * Traits for LazySingleTrackedPartitionedGrid.
 *
 * @tparam T type of value stored in leaf grid nodes.
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<TrackedPartitionedGrid<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = TrackedPartitionedGrid<T, D, N>;
	/// Child grid class, in this case LazySingleTrackedGrid.
	using ChildType = LazyTrackedGrid<T, D, N>;
	/// The type of lookup grid to use for tracking active grid nodes.
	using LookupType = typename ChildType::Traits::LookupType;
	/// Base grid class to partition, in this case TrackedGridBase.
	using MixinType = TrackedGridBase<ThisType>;
	/// Number of tracking lists, from template parameter.
	static const UINT NumLists = N;
	/// Parent grid is eagerly constructed.
	static const Laziness IsLazy = Laziness::EAGER;
};

} // End namespace felt.
#endif /* INCLUDE_FELT_SINGLETRACKEDPARTITIONEDGRID_HPP_ */
