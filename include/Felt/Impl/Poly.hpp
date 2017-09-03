#ifndef INCLUDE_FELT_IMPL_POLY_HPP_
#define INCLUDE_FELT_IMPL_POLY_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/PolyMixin.hpp>
#include <Felt/Impl/Mixin/TrackedMixin.hpp>
#include <Felt/Surface.hpp>

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
	/// IsoGrid type that will be polygonised.
	using IsoGridType = TIsoGrid;
};

} // Impl.

} // Felt.

#endif /* INCLUDE_FELT_IMPL_POLY_HPP_ */
