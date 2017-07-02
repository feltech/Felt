#ifndef FELT_PUBLIC_LOOKUP_HPP_
#define FELT_PUBLIC_LOOKUP_HPP_

#include <type_traits>
#include <Felt/Impl/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/LookupMixin.hpp>

namespace Felt
{
namespace Impl
{
namespace Lookup
{

template <UINT D>
class Simple :
	FELT_MIXINS(
		(Simple<D>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Grid::Size)
		(Lookup::Simple),
		(Grid::Index)
	)
private:
	using ThisType = Simple<D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Simple<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Simple(const VecDi& size_, const VecDi& offset_) :
		SizeImpl{size_, offset_}, ActivatorImpl{NULL_IDX}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class Single :
	FELT_MIXINS(
		(Single<D,N>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Grid::Size)
		(Lookup::Single),
		(Grid::Index)
	)
private:
	using ThisType = Single<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Single(const VecDi& size, const VecDi& offset) :
		SizeImpl{size, offset}, ActivatorImpl{NULL_IDX}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class LazySingle :
	FELT_MIXINS(
		(LazySingle<D,N>),
		(Grid::Accessor::LazyByValue)(Grid::Accessor::Ref)(Grid::Data)(Grid::Resize)
		(Lookup::Activator)
		(Lookup::Single),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)(Grid::Size)
	)
private:
	using ThisType = LazySingle<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Resize<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	using LookupImpl::PosArray;

public:
	LazySingle() :
		ActivatorImpl{NULL_IDX}
	{}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using ActivatorImpl::activate;
	using ActivatorImpl::deactivate;
	using ActivatorImpl::is_active;
	using DataImpl::data;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
};


template <UINT D, UINT N>
class Multi :
	FELT_MIXINS(
		(Multi<D, N>),
		(Grid::Accessor::ByRef)(Grid::Data)(Grid::Size)(Lookup::Activator)(Lookup::Multi),
		(Grid::Activator)(Grid::Index)
	)
private:
	using ThisType = Multi<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Multi<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	using LookupImpl::PosArray;

public:
	Multi(const VecDi& size, const VecDi& offset) :
		SizeImpl{size, offset}, ActivatorImpl{LookupImpl::NULL_IDX_TUPLE}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using ActivatorImpl::activate;
	using ActivatorImpl::deactivate;
	using DataImpl::data;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
	using LookupImpl::NULL_IDX_TUPLE;
};

} // Lookup

/**
 * Traits for Lookup::Simple.
 *
 * @tparam D number of dimensions of the grid.
 */
template <UINT D>
struct Traits< Lookup::Simple<D> >
{
	using LeafType = UINT;
	static constexpr UINT Dims = D;
};

/**
 * Common base traits for lookup grids with multiple tracking lists but single index per grid node.
 */
template <UINT D, UINT N>
struct DefaultLookupTraits
{
	/// Single index stored in each grid node.
	using LeafType = UINT;
	/// Dimension of grid.
	static constexpr UINT Dims = D;
	/// Number of lists tracking grid nodes.
	static constexpr UINT NumLists = N;
};

/**
 * Traits for Lookup::Single.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::Single<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::LazySingle.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::LazySingle<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::Multi.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::Multi<D, N> > : public DefaultLookupTraits<D, N>
{
	/// Multiple indices stored at each grid node, one per tracking list.
	using LeafType = VecDu<D>;
};

} // Impl
} // Felt

#endif
