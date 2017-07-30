#ifndef INCLUDE_FELT_IMPL_PARTITIONED_HPP_
#define INCLUDE_FELT_IMPL_PARTITIONED_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Tracked.hpp>
#include <Felt/Impl/Mixin/PartitionedMixin.hpp>
#include <Felt/Impl/Mixin/NumericMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Partitioned
{

template <Dim D, TupleIdx N>
class Lookup :
	FELT_MIXINS(
		(Lookup<D, N>),
		(Grid::Size)(Partitioned::Children)(Partitioned::Lookup)
	)
private:
	using VecDi = Felt::VecDi<D>;

	using ThisType = Lookup<D, N>;
	using TraitsType = Traits<ThisType>;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<ThisType>;
	using LookupImpl = Impl::Mixin::Partitioned::Lookup<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

public:
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using ChildType = typename TraitsType::ChildType;

	Lookup(const VecDi& size_, const VecDi& offset_, const VecDi& child_size_) :
		SizeImpl{size_, offset_},
		ChildrenImpl{size_, offset_, child_size_, ChildType()}
	{}

	using LookupImpl::track;
	using ChildrenImpl::children;
	using ChildrenImpl::leafs;
	using ChildrenImpl::reset;
};


namespace Tracked
{

template <typename T, Dim D, TupleIdx N>
class Simple :
	FELT_MIXINS(
		(Simple<T, D, N>),
		(Grid::Size)(Partitioned::Children)(Partitioned::Tracked)
	)
private:
	using VecDi = Felt::VecDi<D>;

	using ThisType = Simple<T, D, N>;
	using TraitsType = Traits<ThisType>;
	using LeafType = typename TraitsType::LeafType;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using TrackedImpl = Impl::Mixin::Partitioned::Tracked<ThisType>;
public:
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using ChildType = typename TraitsType::ChildType;

	Simple(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const LeafType background_
	) :
		SizeImpl{size_, offset_},
		ChildrenImpl{size_, offset_, child_size_, ChildType(background_)}
	{}

	using ChildrenImpl::children;
	using ChildrenImpl::reset;
	using ChildrenImpl::track_children;
	using TrackedImpl::track;
};


template <typename T, Dim D, TupleIdx N>
class Numeric :
	FELT_MIXINS(
		(Numeric<T, D, N>),
		(Grid::Size)(Numeric::Spatial)(Partitioned::Accessor)(Partitioned::Children)
		(Partitioned::Snapshot)(Partitioned::Tracked)(Partitioned::Untrack)
	)
private:
	using VecDi = Felt::VecDi<D>;

	using ThisType = Numeric<T, D, N>;
	using TraitsType = Traits<ThisType>;
	using LeafType = typename TraitsType::LeafType;

	using AccessorImpl = Impl::Mixin::Partitioned::Accessor<ThisType>;
	using ChildrenImpl = Impl::Mixin::Partitioned::Children<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using SpatialImpl = Impl::Mixin::Numeric::Spatial<ThisType>;
	using SnapshotImpl = Impl::Mixin::Partitioned::Snapshot<ThisType>;
	using TrackedImpl = Impl::Mixin::Partitioned::Tracked<ThisType>;
	using UntrackImpl = Impl::Mixin::Partitioned::Untrack<ThisType>;
public:
	using ChildType = typename TraitsType::ChildType;
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using SnapshotPtr = typename SnapshotImpl::SnapshotGridPtr;
public:
	Numeric(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const LeafType background_
	) :
		SizeImpl{size_, offset_},
		ChildrenImpl{size_, offset_, child_size_, ChildType(background_)}
	{}

	using AccessorImpl::get;
	using AccessorImpl::set;
	using ChildrenImpl::children;
	using ChildrenImpl::pos_child;
	using ChildrenImpl::reset;
	using ChildrenImpl::track_children;
	using ChildrenImpl::leafs;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::size;
	using SpatialImpl::curv;
	using SpatialImpl::divergence;
	using SpatialImpl::dx;
	using SpatialImpl::get;
	using SpatialImpl::grad;
	using SpatialImpl::gradB;
	using SpatialImpl::gradC;
	using SpatialImpl::gradE;
	using SpatialImpl::gradF;
	using SpatialImpl::interp;
	using SpatialImpl::neighs;
	using SnapshotImpl::operator=;
	using SnapshotImpl::snapshot;
	using TrackedImpl::track;
	using UntrackImpl::untrack;
	using UntrackImpl::retrack;
};

} // Tracked.
} // Partitioned.


/**
 * Traits for Partitioned::Lookup.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Partitioned::Lookup<D, N> > : public DefaultLookupTraits<D, N>
{
	using ChildType = Impl::Lookup::LazySingle<D, N>;
};


/**
 * Traits for Partitioned::Tracked::simple.
 *
 * @tparam T type to store in leaf nodes.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Partitioned::Tracked::Simple<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	using ChildType = Impl::Tracked::LazySingleByValue<T, D, N>;
};

/**
 * Traits for Partitioned::Tracked::Numeric.
 *
 * @tparam T type to store in leaf nodes.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Partitioned::Tracked::Numeric<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	using ChildType = Impl::Tracked::LazySingleByValue<T, D, N>;
};

} // Impl.
} // Felt.



#endif /* INCLUDE_FELT_IMPL_PARTITIONED_HPP_ */
