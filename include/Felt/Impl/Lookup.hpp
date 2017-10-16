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
class SingleListSingleIdx :
	FELT_MIXINS(
		(SingleListSingleIdx<D>),
		(Grid::Access::ByValue)(Grid::Access::Ref)(Grid::Activate)(Grid::Data)(Grid::Size)
		(Lookup::SingleList::SingleIdx),
		(Grid::Index)
	)
private:
	using This = SingleListSingleIdx<D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupImpl = Impl::Mixin::Lookup::SingleList::SingleIdx<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	SingleListSingleIdx(const VecDi& size_, const VecDi& offset_) :
		ActivateImpl{null_idx}, SizeImpl{size_, offset_}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
};


template <Dim D>
class LazySingleListSingleIdx :
	FELT_MIXINS(
		(LazySingleListSingleIdx<D>),
		(Grid::Access::ByValue)(Grid::Access::Ref)(Grid::Resize)(Lookup::SingleList::Activate)
		(Grid::Data)(Lookup::SingleList::SingleIdx),
		(Grid::Activate)(Grid::Size)(Grid::Index)
	)
private:
	using This = LazySingleListSingleIdx<D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Lookup::SingleList::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupImpl = Impl::Mixin::Lookup::SingleList::SingleIdx<This>;
	using ResizeImpl = Impl::Mixin::Grid::Resize<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	LazySingleListSingleIdx() :
		ActivateImpl{null_idx}
	{}

	using AccessImpl::get;
	using AccessImpl::index;
	using ActivateImpl::activate;
	using ActivateImpl::deactivate;
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
class MultiListSingleIdx :
	FELT_MIXINS(
		(MultiListSingleIdx<D,N>),
		(Grid::Access::ByValue)(Grid::Access::Ref)(Grid::Activate)(Grid::Data)(Grid::Size)
		(Lookup::MultiList::SingleIdx),
		(Grid::Index)
	)
private:
	using This = MultiListSingleIdx<D, N>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupImpl = Impl::Mixin::Lookup::MultiList::SingleIdx<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	MultiListSingleIdx(const VecDi& size, const VecDi& offset) :
		ActivateImpl{null_idx}, SizeImpl{size, offset}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
};


template <Dim D, TupleIdx N>
class LazyMultiListSingleIdx :
	FELT_MIXINS(
		(LazyMultiListSingleIdx<D,N>),
		(Grid::Access::LazyByValue)(Grid::Access::Ref)(Grid::Data)(Grid::Resize)
		(Lookup::MultiList::Activate)
		(Lookup::MultiList::SingleIdx),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Index)(Grid::Size)
	)
private:
	using This = LazyMultiListSingleIdx<D, N>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::LazyByValue<This>;
	using ActivateImpl = Impl::Mixin::Lookup::MultiList::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupImpl = Impl::Mixin::Lookup::MultiList::SingleIdx<This>;
	using SizeImpl = Impl::Mixin::Grid::Resize<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	static constexpr TupleIdx num_lists = Traits::t_num_lists;

	LazyMultiListSingleIdx() :
		ActivateImpl{null_idx}
	{}

	/**
	 * Serialisation hook for cereal library.
	 *
	 * @param ar
	 */
	template<class Archive>
	void serialize(Archive & ar)
	{
		ActivateImpl::serialize(ar);
		DataImpl::serialize(ar);
		LookupImpl::serialize(ar);
		SizeImpl::serialize(ar);
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using ActivateImpl::activate;
	using ActivateImpl::deactivate;
	using ActivateImpl::is_active;
	using DataImpl::assert_pos_idx_bounds;
	using DataImpl::data;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
};


template <Dim D, TupleIdx N>
class MultiListMultiIdx :
	FELT_MIXINS(
		(MultiListMultiIdx<D, N>),
		(Grid::Access::ByRef)(Grid::Data)(Grid::Size)(Lookup::MultiList::Activate)
		(Lookup::MultiList::MultiIdx),
		(Grid::Activate)(Grid::Index)
	)
private:
	using This = MultiListMultiIdx<D, N>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByRef<This>;
	using ActivateImpl = Impl::Mixin::Lookup::MultiList::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupImpl = Impl::Mixin::Lookup::MultiList::MultiIdx<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	static constexpr TupleIdx num_lists = Traits::t_num_lists;

	MultiListMultiIdx() = default;

	MultiListMultiIdx(const VecDi& size, const VecDi& offset) :
		SizeImpl{size, offset}, ActivateImpl{LookupImpl::s_null_idxs}
	{
		this->activate();
	}

	template<class Archive>
	void serialize(Archive & ar)
	{
		ActivateImpl::serialize(ar);
		DataImpl::serialize(ar);
		LookupImpl::serialize(ar);
		SizeImpl::serialize(ar);
	}


	using AccessImpl::get;
	using AccessImpl::index;
	using ActivateImpl::activate;
	using ActivateImpl::deactivate;
	using DataImpl::data;
	using LookupImpl::track;
	using LookupImpl::is_tracked;
	using LookupImpl::list;
	using LookupImpl::untrack;
	using LookupImpl::reset;
	using LookupImpl::s_null_idxs;
};

} // Lookup
} // Impl
} // Felt


namespace Felt
{
namespace Impl
{
/**
 * Traits for Lookup::Simple.
 *
 * @tparam D number of dimensions of the grid.
 */
template <Dim D>
struct Traits< Lookup::SingleListSingleIdx<D> >
{
	using Leaf = ListIdx;
	static constexpr UINT t_dims = D;
};

/**
 * Traits for Lookup::LazySimple.
 *
 * @tparam D number of dimensions of the grid.
 */
template <Dim D>
struct Traits< Lookup::LazySingleListSingleIdx<D> >
{
	using Leaf = ListIdx;
	static constexpr UINT t_dims = D;
};

/**
 * Common base traits for lookup grids with multiple tracking lists but single index per grid node.
 */
template <Dim D, TupleIdx N>
struct DefaultLookupTraits
{
	/// Single index stored in each grid node.
	using Leaf = ListIdx;
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
struct Traits< Lookup::MultiListSingleIdx<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::LazySingle.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Lookup::LazyMultiListSingleIdx<D, N> > : public DefaultLookupTraits<D, N>
{};

/**
 * Traits for Lookup::Multi.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <Dim D, TupleIdx N>
struct Traits< Lookup::MultiListMultiIdx<D, N> > : public DefaultLookupTraits<D, N>
{
	/// Multiple indices stored at each grid node, one per tracking list.
	using Leaf = Tuple<ListIdx, N>;
};


} // Impl
} // Felt

#endif
