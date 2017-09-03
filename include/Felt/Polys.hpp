#ifndef INCLUDE_FELT_POLYS_HPP_
#define INCLUDE_FELT_POLYS_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/PolyMixin.hpp>
#include <Felt/Impl/Poly.hpp>

namespace Felt
{
template <class TSurface>
class Polys :
	FELT_MIXINS(
		(Polys<TSurface>),
		(Poly::Children)(Poly::Update)
	)
private:
	using ThisType = Polys<TSurface>;
	using TraitsType = Impl::Traits<ThisType>;
	using ChildType = typename TraitsType::ChildType;

	using ChildrenImpl = Impl::Mixin::Poly::Children<ThisType>;
	using UpdateImpl = Impl::Mixin::Poly::Update<ThisType>;
public:
	Polys(const TSurface& surface_) :
		ChildrenImpl{surface_.isogrid()}, UpdateImpl{surface_}
	{}

	using ChildrenImpl::children;
	using UpdateImpl::changes;
	using UpdateImpl::invalidate;
	using UpdateImpl::notify;
	using UpdateImpl::march;
};

namespace Impl
{
/**
 * Traits for Polys.
 *
 * @tparam SurfaceType surface type to polygonise.
 */
template <class TSurface>
struct Traits< Polys<TSurface> >
{
	/// Type of surface to polygonise.
	using SurfaceType = TSurface;
	/// Isogrid type that the surface wraps.
	using IsoGridType = typename SurfaceType::IsoGrid;
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<IsoGridType>::t_dims;
	/// Child poly type to polygonise a single spatial partition.
	using ChildType = Impl::Poly::Single<IsoGridType>;
	/// Children grid type to store and track active child polys.
	using ChildrenType = Impl::Tracked::SingleListSingleIdxByRef<ChildType, t_dims>;
};
} // Impl.

} // Felt.

#endif /* INCLUDE_FELT_POLYS_HPP_ */
