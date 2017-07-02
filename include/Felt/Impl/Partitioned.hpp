#ifndef INCLUDE_FELT_IMPL_PARTITIONED_HPP_
#define INCLUDE_FELT_IMPL_PARTITIONED_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Mixin/PartitionedMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Partitioned
{

template <UINT D, UINT N>
class Lookup :
	FELT_MIXINS(
		(Lookup<D, N>),
		(Grid::Size)(Partitioned::Children)(Partitioned::Lookup)
	)
private:
	using ThisType = Lookup<D, N>;
	using TraitsType = Traits<ThisType>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<ThisType>;
	using LookupImpl = Impl::Mixin::Partitioned::Lookup<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using VecDi = Felt::VecDi<D>;

public:
	using ChildrenGrid = typename ChildrenImpl::ChildrenGrid;
public:
	Lookup(const VecDi& size_, const VecDi& offset_, const VecDi& child_size_) :
		ChildrenImpl{size_, offset_, child_size_, ChildType()}, SizeImpl{size_, offset_}
	{}

	using LookupImpl::track;
	using ChildrenImpl::children;
	using ChildrenImpl::reset;
};

} // Partitioned.


/**
 * Traits for Lookup::Partitioned.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Partitioned::Lookup<D, N> > : public DefaultLookupTraits<D, N>
{
	using ChildType = Impl::Lookup::LazySingle<D, N>;
};

} // Impl.
} // Felt.



#endif /* INCLUDE_FELT_IMPL_PARTITIONED_HPP_ */
