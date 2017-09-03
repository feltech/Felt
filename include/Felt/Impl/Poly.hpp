#ifndef INCLUDE_FELT_IMPL_POLY_HPP_
#define INCLUDE_FELT_IMPL_POLY_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/PolyMixin.hpp>
#include <Felt/Impl/Mixin/TrackedMixin.hpp>
#include <Felt/Impl/Surface.hpp>

namespace Felt
{
namespace Impl
{
namespace Poly
{

template <class TIsoGrid>
class Single :
	FELT_MIXINS(
		(Single<TIsoGrid>),
		(Grid::Access::ByRef)(Grid::Data)(Poly::Activate)(Poly::Reset)(Poly::Resize)(Poly::March)
		(Tracked::LookupInterface),
		(Grid::Activate)(Grid::Index)(Grid::Resize)(Grid::Size)(Tracked::Resize)
		(Tracked::Activate)(Tracked::SingleList::Reset)
	)
private:
	using ThisType = Single<TIsoGrid>;
	using TraitsType = Impl::Traits<ThisType>;

	using ActivateImpl = Impl::Mixin::Poly::Activate<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using MarchImpl = Impl::Mixin::Poly::March<ThisType>;
	using ResetterImpl = Impl::Mixin::Poly::Reset<ThisType>;
	using SizeImpl = Impl::Mixin::Poly::Resize<ThisType>;

	using LookupType = typename TraitsType::LookupType;
	using IsoGridType = typename TraitsType::IsoGridType;

public:
	using Simplex = typename TraitsType::Simplex;

	Single(const IsoGridType& isogrid_) :
		MarchImpl{isogrid_}, LookupInterfaceImpl{LookupType{}}
	{}

	using ActivateImpl::activate;
	using ActivateImpl::deactivate;
	using ActivateImpl::is_active;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
	using MarchImpl::bind;
	using MarchImpl::march;
	using MarchImpl::spxs;
	using MarchImpl::vtxs;
	using ResetterImpl::reset;
};


template <class TSurface>
class Grid :
	FELT_MIXINS(
		(Grid<TSurface>),
		(Poly::Children)(Poly::Update)
	)
private:
	using ThisType = Grid<TSurface>;
	using TraitsType = Traits<ThisType>;
	using ChildType = typename TraitsType::ChildType;

	using ChildrenImpl = Impl::Mixin::Poly::Children<ThisType>;
	using UpdateImpl = Impl::Mixin::Poly::Update<ThisType>;
public:
	Grid(const TSurface& surface_) :
		ChildrenImpl{surface_.isogrid()}, UpdateImpl{surface_}
	{}

	using ChildrenImpl::children;
	using UpdateImpl::changes;
	using UpdateImpl::invalidate;
	using UpdateImpl::notify;
	using UpdateImpl::march;
};


} // Poly.


/**
 * Traits for Poly::Single.
 *
 * @tparam IsoGrid isogrid type to polygonise.
 */
template <class TIsoGrid>
struct Traits< Poly::Single<TIsoGrid> > : Mixin::Poly::Traits<Traits<TIsoGrid>::t_dims, void>
{
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<TIsoGrid>::t_dims;
	/// A vertex index for each positively directed edge stored at each grid node.
	using LeafType = Felt::VecDu<t_dims>;
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::LazySingleListSingleIdx<t_dims>;
	/// TIsoGrid type that will be polygonised.
	using IsoGridType = TIsoGrid;
};


/**
  * Traits for Poly::Grid.
 *
 * @tparam SurfaceType surface type to polygonise.
 */
template <class TSurface>
struct Traits< Poly::Grid<TSurface> >
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

#endif /* INCLUDE_FELT_IMPL_POLY_HPP_ */
