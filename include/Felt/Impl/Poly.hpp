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

template <class IsoGrid>
class Single :
	FELT_MIXINS(
		(Single<IsoGrid>),
		(Grid::Access::ByRef)(Grid::Data)(Poly::Activate)(Poly::Reset)(Poly::Resize)(Poly::March)
		(Tracked::LookupInterface),
		(Grid::Activate)(Grid::Index)(Grid::Resize)(Grid::Size)(Tracked::Resize)
		(Tracked::Activate)(Tracked::Single::Reset)
	)
private:
	using ThisType = Single<IsoGrid>;
	using TraitsType = Impl::Traits<ThisType>;

	using ActivateImpl = Impl::Mixin::Poly::Activate<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using MarchImpl = Impl::Mixin::Poly::March<ThisType>;
	using ResetterImpl = Impl::Mixin::Poly::Reset<ThisType>;
	using SizeImpl = Impl::Mixin::Poly::Resize<ThisType>;

	using LookupType = typename TraitsType::LookupType;
	using IsoGridType = typename TraitsType::IsoGridType;

public:
	Single(const IsoGridType& isogrid_) :
		MarchImpl{isogrid_}, LookupInterfaceImpl{LookupType{}}
	{}

	using ActivateImpl::activate;
	using ActivateImpl::deactivate;
	using SizeImpl::resize;
	using MarchImpl::bind;
	using MarchImpl::march;
	using MarchImpl::spxs;
	using MarchImpl::vtxs;
	using ResetterImpl::reset;
};


template <class SurfaceType>
class Grid :
	FELT_MIXINS(
		(Grid<SurfaceType>),
		(Partitioned::Children)
	)
private:
	using ThisType = Grid<SurfaceType>;
	using TraitsType = Traits<ThisType>;
	using ChildType = typename TraitsType::ChildType;

	using ChildrenImpl = Impl::Mixin::Partitioned::Children<ThisType>;
public:
	Grid(const SurfaceType& surface_) :
		ChildrenImpl{
			surface_.isogrid().size(), surface_.isogrid().offset(),
			surface_.isogrid().child_size(), ChildType(surface_.isogrid())
		}
	{}

	using ChildrenImpl::children;
};


} // Poly.


/**
 * Traits for Poly::Single.
 *
 * @tparam IsoGrid isogrid type to polygonise.
 */
template <class IsoGrid>
struct Traits< Poly::Single<IsoGrid> > : Mixin::Poly::Traits<Traits<IsoGrid>::t_dims, void>
{
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<IsoGrid>::t_dims;
	/// A vertex index for each positively directed edge stored at each grid node.
	using LeafType = Felt::VecDu<t_dims>;
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::LazySimple<t_dims>;
	/// Isogrid type that will be polygonised.
	using IsoGridType = IsoGrid;
};


/**
  * Traits for Poly::Grid.
 *
 * @tparam SurfaceType surface type to polygonise.
 */
template <class SurfaceType>
struct Traits< Poly::Grid<SurfaceType> >
{
	using IsoGridType = typename SurfaceType::IsoGrid;
	/// Dimension of grid.
	static constexpr Dim t_dims = Traits<IsoGridType>::t_dims;
	/// Child poly type to polygonise a single spatial partition.
	using ChildType = Impl::Poly::Single<IsoGridType>;
	using ChildrenType = Impl::Tracked::SingleByRef<ChildType, t_dims>;
};

} // Impl.
} // Felt.

#endif /* INCLUDE_FELT_IMPL_POLY_HPP_ */
