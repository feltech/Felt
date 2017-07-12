#ifndef FELT_PUBLIC_TRACKED_HPP_
#define FELT_PUBLIC_TRACKED_HPP_

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/TrackedMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Tracked
{

template <typename T, UINT D, UINT N>
class LazySingleByValue :
	FELT_MIXINS(
		(LazySingleByValue<T, D, N>),
		(Grid::Accessor::LazyByValue)(Grid::Data)(Tracked::Activator)(Tracked::ByValue)
		(Tracked::LookupInterface)(Tracked::Resetter)(Tracked::Resize),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using ThisType = LazySingleByValue<T, D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Tracked::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using ResetterImpl = Impl::Mixin::Tracked::Resetter<ThisType>;
	using SizeImpl = Impl::Mixin::Tracked::Resize<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::ByValue<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;
	using LookupType = typename ThisTraits::LookupType;

public:
	LazySingleByValue(const LeafType background_) :
		ActivatorImpl{background_}, LookupInterfaceImpl{LookupType{}}
	{}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using AccessorImpl::set;
	using ActivatorImpl::activate;
	using ActivatorImpl::background;
	using ActivatorImpl::deactivate;
	using ActivatorImpl::is_active;
	using DataImpl::data;
	using LookupInterfaceImpl::lookup;
	using ResetterImpl::reset;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
	using TrackedImpl::track;

};


template <typename T, UINT D, UINT N>
class MultiByRef :
	FELT_MIXINS(
		(MultiByRef<T, D, N>),
		(Grid::Accessor::ByRef)(Grid::Activator)(Grid::Data)(Grid::Size)(Tracked::ByRef)
		(Tracked::LookupInterface),
		(Grid::Index)
	)
private:
	using ThisType = MultiByRef<T, D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::ByRef<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;
	using LookupType = typename ThisTraits::LookupType;

public:
	using LookupInterfaceImpl::PosArray;
	static constexpr ListIdx NumLists = ThisTraits::NumLists;

public:
	MultiByRef(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		ActivatorImpl{background_}, LookupInterfaceImpl{LookupType{size_, offset_}},
		SizeImpl{size_, offset_}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using ActivatorImpl::activate;
	using DataImpl::data;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::offset;
	using SizeImpl::size;
	using LookupInterfaceImpl::lookup;
	using TrackedImpl::track;
};

} // Tracked.


/**
 * Common base traits for tracked grids with multiple tracking lists but single index per grid node.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, UINT D, UINT N>
struct DefaultTrackedTraits
{
	/// Single index stored in each grid node.
	using LeafType = T;
	/// Dimension of grid.
	static constexpr UINT Dims = D;
	/// Number of lists tracking grid nodes.
	static constexpr UINT NumLists = N;
};


/**
 * Traits for Tracked::LazySingle.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, UINT D, UINT N>
struct Traits< Tracked::LazySingleByValue<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	using LookupType = Lookup::LazySingle<D, N>;
};

/**
 * Traits for Tracked::Multi.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, UINT D, UINT N>
struct Traits< Tracked::MultiByRef<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	using LookupType = Lookup::Multi<D, N>;
};

} // Impl
} // Felt



#endif /* FELT_PUBLIC_TRACKED_HPP_ */
