
#ifndef INCLUDE_FELT_LOOKUPPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_LOOKUPPARTITIONEDGRID_HPP_

#include "TrackingPartitionedGridBase.hpp"

namespace felt
{

/**
 * Spatially partitioned wrapper for MultiLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N=1>
class MultiLookupPartitionedGrid
	: public TrackingPartitionedGridBase<MultiLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = MultiLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::LeafType;
	using Base::TrackingPartitionedGridBase;

	/**
	 * Construct a spatially partitioned MultiLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	MultiLookupPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_
	) : Base(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_)
	{}
};


/**
 * Spatially partitioned wrapper for SingleLookupGrid.
 *
 * @snippet test_PartitionedGrid.cpp LazySingleLookupPartitionedGrid initialisation
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class SingleLookupPartitionedGrid
	: public TrackingPartitionedGridBase<SingleLookupPartitionedGrid<D, N> >
{
public:
	using ThisType = SingleLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;

	using Base::TrackingPartitionedGridBase;

	/**
	 * Construct a spatially partitioned SingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	SingleLookupPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_
	) :
		Base(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_)
	{}
};


/**
 * Spatially partitioned wrapper for LazySingleLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class LazySingleLookupPartitionedGrid
	: public TrackingPartitionedGridBase<LazySingleLookupPartitionedGrid<D, N>, true>
{
public:
	using ThisType = LazySingleLookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType, true>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::Child;

	LazySingleLookupPartitionedGrid() : Base()
	{}

	/**
	 * Construct a spatially partitioned LazySingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	LazySingleLookupPartitionedGrid (
		const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_
	) : Base()
	{
		this->init(size_, offset_, partition_size_);
	}

	/**
	 * Initialise a spatially partitioned LazySingleLookupGrid.
	 *
	 * Initialises grid data to NULL_IDX.
	 *
	 * @param size_ size of the grid.
	 * @param offset_ spatial offset of the grid.
	 * @param partition_size_ size of a spatial partition.
	 */
	void init (const VecDu& size_, const VecDi& offset_, const VecDu& partition_size_)
	{
		Base::init(size_, offset_, Base::Child::Traits::NULL_IDX_DATA, partition_size_);
	}
};


/**
 * Traits for MultiLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<MultiLookupPartitionedGrid<D, N> > : DefaultMultiLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = MultiLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case MultiLookupGrid.
	using ChildType = MultiLookupGrid<D, N>;
	/// Base grid class to partition, in this case LookupGridBase.
	using MixinType = LookupGridBase<ThisType>;
};


/**
 * Traits for SingleLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<SingleLookupPartitionedGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = SingleLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case SingleLookupGrid.
	using ChildType = SingleLookupGrid<D, N>;
	/// Grid class whose interface to copy via CRTP mixin.
	using MixinType = SingleLookupGridBase<ThisType>;
};


/**
 * Traits for LazySingleLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LazySingleLookupPartitionedGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = LazySingleLookupPartitionedGrid<D, N>;
	/// Child grid class, in this case LazySingleLookupGrid.
	using ChildType = LazySingleLookupGrid<D, N>;
	/// Grid class whose interface to copy via CRTP mixin.
	using MixinType = SingleLookupGridBase<ThisType>;
};

} // End namespace felt
#endif
