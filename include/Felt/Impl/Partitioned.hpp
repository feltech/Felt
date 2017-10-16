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
		(Grid::Size)(Partitioned::Children)(Partitioned::Leafs)(Partitioned::Lookup)
		(Partitioned::Reset::MultiList)
	)
private:
	using VecDi = Felt::VecDi<D>;

	using This = Lookup<D, N>;
	using Traits = Impl::Traits<This>;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<This>;
	using LeafsImpl = Impl::Mixin::Partitioned::Leafs<This>;
	using LookupImpl = Impl::Mixin::Partitioned::Lookup<This>;
	using ResetImpl = Impl::Mixin::Partitioned::Reset::MultiList<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

public:
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using Child = typename Traits::Child;

	Lookup(const VecDi& size_, const VecDi& offset_, const VecDi& child_size_) :
		SizeImpl{size_, offset_},
		ChildrenImpl{size_, offset_, child_size_, Child()}
	{}

	using LookupImpl::track;
	using ChildrenImpl::children;
	using LeafsImpl::leafs;
	using ResetImpl::reset;
};


namespace Tracked
{

template <typename T, Dim D, TupleIdx N>
class Simple :
	FELT_MIXINS(
		(Simple<T, D, N>),
		(Grid::Size)(Partitioned::Children)(Partitioned::Leafs)(Partitioned::Reset::MultiList)
		(Partitioned::Tracked)
	)
private:
	using VecDi = Felt::VecDi<D>;

	using This = Simple<T, D, N>;
	using Traits = Impl::Traits<This>;
	using Leaf = typename Traits::Leaf;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<This>;
	using LeafsImpl = Impl::Mixin::Partitioned::Leafs<This>;
	using ResetImpl = Impl::Mixin::Partitioned::Reset::MultiList<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using TrackedImpl = Impl::Mixin::Partitioned::Tracked<This>;
public:
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using Child = typename Traits::Child;

	Simple(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const Leaf background_
	) :
		SizeImpl{size_, offset_},
		ChildrenImpl{size_, offset_, child_size_, Child(background_)}
	{}

	using ChildrenImpl::children;
	using ResetImpl::reset;
	using ResetImpl::track_children;
	using TrackedImpl::track;
};


template <typename T, Dim D, TupleIdx N>
class Numeric :
	FELT_MIXINS(
		(Numeric<T, D, N>),
		(Grid::Size)(Numeric::Spatial)(Partitioned::Access)(Partitioned::Children)
		(Partitioned::Leafs)(Partitioned::Reset::MultiList)(Partitioned::Snapshot)
		(Partitioned::Tracked)(Partitioned::Untrack)
	)
	friend class cereal::access;
private:
	using VecDi = Felt::VecDi<D>;

	using This = Numeric<T, D, N>;
	using Traits = Impl::Traits<This>;
	using Leaf = typename Traits::Leaf;

	using AccessImpl = Impl::Mixin::Partitioned::Access<This>;
	using ChildrenImpl = Impl::Mixin::Partitioned::Children<This>;
	using LeafsImpl = Impl::Mixin::Partitioned::Leafs<This>;
	using ResetImpl = Impl::Mixin::Partitioned::Reset::MultiList<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using SpatialImpl = Impl::Mixin::Numeric::Spatial<This>;
	using SnapshotImpl = Impl::Mixin::Partitioned::Snapshot<This>;
	using TrackedImpl = Impl::Mixin::Partitioned::Tracked<This>;
	using UntrackImpl = Impl::Mixin::Partitioned::Untrack<This>;
public:
	using Child = typename Traits::Child;
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
	using SnapshotPtr = typename SnapshotImpl::SnapshotGridPtr;

	static constexpr TupleIdx num_lists = N;

	Numeric(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const Leaf background_
	) :
		SizeImpl{size_, offset_},
		AccessImpl{background_},
		ChildrenImpl{size_, offset_, child_size_, Child{background_}}
	{}

	Numeric(Numeric&& other_) :
		SizeImpl{std::move(other_)},
		SpatialImpl{std::move(other_)},
		AccessImpl{std::move(other_)},
		ChildrenImpl{std::move(other_)},
		LeafsImpl{std::move(other_)},
		ResetImpl{std::move(other_)},
		SnapshotImpl{std::move(other_)},
		TrackedImpl{std::move(other_)},
		UntrackImpl{std::move(other_)}
	{}

	/**
	 * Serialisation hook for cereal library.
	 *
	 * @param ar
	 */
	template<class Archive>
	void serialize(Archive & ar)
	{
		AccessImpl::serialize(ar);
		ChildrenImpl::serialize(ar);
		SizeImpl::serialize(ar);
		SpatialImpl::serialize(ar);
	}

	using AccessImpl::get;
	using AccessImpl::set;
	using ChildrenImpl::child_size;
	using ChildrenImpl::children;
	using LeafsImpl::pos_child;
	using LeafsImpl::pos_idx_child;
	using LeafsImpl::leafs;
	using ResetImpl::reset;
	using ResetImpl::track_children;
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
	using SnapshotImpl::load;
	using SnapshotImpl::operator=;
	using SnapshotImpl::save;
	using SnapshotImpl::snapshot;
	using TrackedImpl::track;
	using UntrackImpl::untrack;
	using UntrackImpl::retrack;

private:
	Numeric() {};
};


} // Tracked.
} // Partitioned.
} // Impl.
} // Felt.


namespace Felt
{
namespace Impl
{
/**
 * Traits for Partitioned::Lookup.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Partitioned::Lookup<D, N> > : public DefaultLookupTraits<D, N>
{
	using Child = Impl::Lookup::LazyMultiListSingleIdx<D, N>;
	using Children = Impl::Tracked::MultiListMultiIdxByRef<Child, D, N>;
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
	using Child = Impl::Tracked::LazyMultiListSingleIdxByValue<T, D, N>;
	using Children = Impl::Tracked::MultiListMultiIdxByRef<Child, D, N>;
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
	using Child = Impl::Tracked::LazyMultiListSingleIdxByValue<T, D, N>;
	using Children = Impl::Tracked::MultiListMultiIdxByRef<Child, D, N>;
};

} // Impl.
} // Felt.



#endif /* INCLUDE_FELT_IMPL_PARTITIONED_HPP_ */
