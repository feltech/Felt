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
	using MarchImpl::march;
	using MarchImpl::spxs;
	using MarchImpl::vtxs;
	using ResetterImpl::reset;
};

} // Poly.


/**
 * Traits for Poly.
 *
 * @tparam D number of dimensions of the isogrid.
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


} // Impl.
} // Felt.

#endif /* INCLUDE_FELT_IMPL_POLY_HPP_ */
