
#ifndef INCLUDE_FELT_SINGLELOOKUPPARTITIONEDGRID_HPP_
#define INCLUDE_FELT_SINGLELOOKUPPARTITIONEDGRID_HPP_

#include "TrackingPartitionedGridBase.hpp"


namespace felt
{

/**
 * Spatially partitioned wrapper for LazySingleLookupGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
class LookupPartitionedGrid
	: public TrackingPartitionedGridBase< LookupPartitionedGrid<D, N> >
{
public:
	using ThisType = LookupPartitionedGrid<D, N>;
	using Base = TrackingPartitionedGridBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::Child;

	LookupPartitionedGrid() : Base()
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
	LookupPartitionedGrid (
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
 * Traits for LazySingleLookupPartitionedGrid.
 *
 * @tparam D dimension of the grid.
 * @tparam N number of tracking lists to use.
 */
template <UINT D, UINT N>
struct GridTraits<LookupPartitionedGrid<D, N> > : DefaultSingleLookupGridTraits<D, N>
{
	/// The class inheriting from TrackingPartitionedGrid.
	using ThisType = LookupPartitionedGrid<D, N>;
	/// Child grid class, in this case LazySingleLookupGrid.
	using ChildType = LazyLookupGrid<D, N>;
	/// Grid class whose interface to copy via CRTP mixin.
	using MixinType = EagerSingleLookupGridBase<ThisType>;
	/// Parent grid is eagerly constructed.
	static const Laziness IsLazy = Laziness::EAGER;
};

} // End namespace felt
#endif
