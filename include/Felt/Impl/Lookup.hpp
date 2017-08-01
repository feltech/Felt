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

template <Dim D>
class Simple :
	FELT_MIXINS(
		(Simple<D>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Grid::Size)
		(Lookup::Single::Single),
		(Grid::Index)
	)
private:
	using ThisType = Simple<D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single::Single<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Simple(const VecDi& size_, const VecDi& offset_) :
		ActivatorImpl{null_idx}, SizeImpl{size_, offset_}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
};


template <Dim D>
class LazySimple :
	FELT_MIXINS(
		(LazySimple<D>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Resize)(Lookup::Single::Activator)
		(Grid::Data)(Lookup::Single::Single),
		(Grid::Activator)(Grid::Size)(Grid::Index)
	)
private:
	using ThisType = LazySimple<D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Single::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single::Single<ThisType>;
	using ResizeImpl = Impl::Mixin::Grid::Resize<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	LazySimple() :
		ActivatorImpl{null_idx}
	{}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using ActivatorImpl::activate;
	using ActivatorImpl::deactivate;
	using DataImpl::data;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
	using ResizeImpl::resize;
	using SizeImpl::offset;
	using SizeImpl::size;
};


template <Dim D, TupleIdx N>
class Single :
	FELT_MIXINS(
		(Single<D,N>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Grid::Size)
		(Lookup::Multi::Single),
		(Grid::Index)
	)
private:
	using ThisType = Single<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Multi::Single<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Single(const VecDi& size, const VecDi& offset) :
		ActivatorImpl{null_idx}, SizeImpl{size, offset}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
};


template <Dim D, TupleIdx N>
class LazySingle :
	FELT_MIXINS(
		(LazySingle<D,N>),
		(Grid::Accessor::LazyByValue)(Grid::Accessor::Ref)(Grid::Data)(Grid::Resize)
		(Lookup::Multi::Activator)
		(Lookup::Multi::Single),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)(Grid::Size)
	)
private:
	using ThisType = LazySingle<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Multi::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Multi::Single<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Resize<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	LazySingle() :
		ActivatorImpl{null_idx}
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
	using LookupImpl::untrack;
	using LookupImpl::reset;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
};


template <Dim D, TupleIdx N>
class Multi :
	FELT_MIXINS(
		(Multi<D, N>),
		(Grid::Accessor::ByRef)(Grid::Data)(Grid::Size)(Lookup::Multi::Activator)
		(Lookup::Multi::Multi),
		(Grid::Activator)(Grid::Index)
	)
private:
	using ThisType = Multi<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Multi::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Multi::Multi<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Multi(const VecDi& size, const VecDi& offset) :
		SizeImpl{size, offset}, ActivatorImpl{LookupImpl::s_null_idxs}
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
	using LookupImpl::untrack;
	using LookupImpl::reset;
	using LookupImpl::s_null_idxs;
};

} // Lookup

/**
 * Traits for Lookup::Simple.
 *
 * @tparam D number of dimensions of the grid.
 */
template <Dim D>
struct Traits< Lookup::Simple<D> >
{
	using LeafType = ListIdx;
	static constexpr UINT t_dims = D;
};

/**
 * Traits for Lookup::LazySimple.
 *
 * @tparam D number of dimensions of the grid.
 */
template <Dim D>
struct Traits< Lookup::LazySimple<D> >
{
	using LeafType = ListIdx;
	static constexpr UINT t_dims = D;
};

/**
 * Common base traits for lookup grids with multiple tracking lists but single index per grid node.
 */
template <Dim D, TupleIdx N>
struct DefaultLookupTraits
{
	/// Single index stored in each grid node.
	using LeafType = ListIdx;
	/// Dimension of grid.
	static constexpr Dim t_dims = D;
	/// Number of lists tracking grid nodes.
	static constexpr TupleIdx t_num_lists = N;
};

/**
 * Traits for Lookup::Single.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Lookup::Single<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::LazySingle.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Lookup::LazySingle<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::Multi.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Lookup::Multi<D, N> > : public DefaultLookupTraits<D, N>
{
	/// Multiple indices stored at each grid node, one per tracking list.
	using LeafType = Tuple<ListIdx, N>;
};


} // Impl
} // Felt

#endif
