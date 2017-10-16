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

template <typename T, Dim D, TupleIdx N>
class LazyMultiListSingleIdxByValue :
	FELT_MIXINS(
		(LazyMultiListSingleIdxByValue<T, D, N>),
		(Grid::Access::LazyByValue)(Grid::Data)(Tracked::Activate)(Tracked::MultiList::ByValue)
		(Tracked::MultiList::LookupInterface)(Tracked::MultiList::Reset)(Tracked::Resize),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using This = LazyMultiListSingleIdxByValue<T, D, N>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::LazyByValue<This>;
	using ActivateImpl = Impl::Mixin::Tracked::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::MultiList::LookupInterface<This>;
	using ResetterImpl = Impl::Mixin::Tracked::MultiList::Reset<This>;
	using SizeImpl = Impl::Mixin::Tracked::Resize<This>;
	using TrackedImpl = Impl::Mixin::Tracked::MultiList::ByValue<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	using Lookup = typename Traits::Lookup;

	LazyMultiListSingleIdxByValue() = default;

	LazyMultiListSingleIdxByValue(const Leaf background_) :
		ActivateImpl{background_}, LookupInterfaceImpl{Lookup{}}
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
		LookupInterfaceImpl::serialize(ar);
		SizeImpl::serialize(ar);
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using AccessImpl::set;
	using ActivateImpl::activate;
	using ActivateImpl::background;
	using ActivateImpl::deactivate;
	using ActivateImpl::is_active;
	using DataImpl::assert_pos_idx_bounds;
	using DataImpl::data;
	using LookupInterfaceImpl::list;
	using LookupInterfaceImpl::lookup;
	using ResetterImpl::reset;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
	using TrackedImpl::track;
};


template <typename T, Dim D>
class SingleListSingleIdxByRef :
	FELT_MIXINS(
		(SingleListSingleIdxByRef<T, D>),
		(Grid::Access::ByRef)(Grid::Activate)(Grid::Data)(Grid::Size)(Tracked::SingleList::ByRef)
		(Tracked::LookupInterface),
		(Grid::Index)
	)
private:
	using This = SingleListSingleIdxByRef<T, D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByRef<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using TrackedImpl = Impl::Mixin::Tracked::SingleList::ByRef<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;
	using Lookup = typename Traits::Lookup;

public:
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	SingleListSingleIdxByRef(const VecDi& size_, const VecDi& offset_, const Leaf background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{Lookup{size_, offset_}}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using ActivateImpl::activate;
	using DataImpl::assert_pos_idx_bounds;
	using DataImpl::data;
	using SizeImpl::offset;
	using SizeImpl::size;
	using LookupInterfaceImpl::lookup;
	using TrackedImpl::track;
};


template <typename T, Dim D>
class SingleListSingleIdxByValue :
	FELT_MIXINS(
		(SingleListSingleIdxByValue<T, D>),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Data)(Grid::Size)
		(Tracked::SingleList::ByValue)
		(Tracked::LookupInterface),
		(Grid::Index)
	)
private:
	using This = SingleListSingleIdxByValue<T, D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using TrackedImpl = Impl::Mixin::Tracked::SingleList::ByValue<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;
	using Lookup = typename Traits::Lookup;

public:
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	SingleListSingleIdxByValue(const VecDi& size_, const VecDi& offset_, const Leaf background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{Lookup{size_, offset_}}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using AccessImpl::set;
	using ActivateImpl::activate;
	using DataImpl::assert_pos_idx_bounds;
	using DataImpl::data;
	using SizeImpl::offset;
	using SizeImpl::size;
	using LookupInterfaceImpl::lookup;
	using TrackedImpl::track;
};


template <typename T, Dim D, TupleIdx N>
class MultiListMultiIdxByRef :
	FELT_MIXINS(
		(MultiListMultiIdxByRef<T, D, N>),
		(Grid::Access::ByRef)(Grid::Activate)(Grid::Data)(Grid::Size)(Tracked::MultiList::ByRef)
		(Tracked::MultiList::LookupInterface),
		(Grid::Index)
	)
private:
	using This = MultiListMultiIdxByRef<T, D, N>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByRef<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::MultiList::LookupInterface<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using TrackedImpl = Impl::Mixin::Tracked::MultiList::ByRef<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;
	using Lookup = typename Traits::Lookup;

public:
	static constexpr TupleIdx t_num_lists = Traits::t_num_lists;

	MultiListMultiIdxByRef() {};

	MultiListMultiIdxByRef(const VecDi& size_, const VecDi& offset_, const Leaf background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{Lookup{size_, offset_}}
	{
		this->activate();
	}

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
		LookupInterfaceImpl::serialize(ar);
		SizeImpl::serialize(ar);
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using ActivateImpl::activate;
	using DataImpl::assert_pos_idx_bounds;
	using DataImpl::data;
	using SizeImpl::offset;
	using SizeImpl::size;
	using LookupInterfaceImpl::lookup;
	using TrackedImpl::track;
};


} // Tracked.
} // Impl
} // Felt


namespace Felt
{
namespace Impl
{

/**
 * Traits for Tracked::SingleListSingleIdxByRef.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Tracked::SingleListSingleIdxByRef<T, D> >
{
	/// Single index stored in each grid node.
	using Leaf = T;
	/// Dimension of grid.
	static constexpr Dim t_dims = D;
	/// Type of lookup grid for tracking active positions.
	using Lookup = Lookup::SingleListSingleIdx<D>;
};

/**
 * Traits for Tracked::SingleListSingleIdxByValue.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Tracked::SingleListSingleIdxByValue<T, D> > :
	Traits< Tracked::SingleListSingleIdxByRef<T, D> >
{};

/**
 * Common base traits for tracked grids with multiple tracking lists but single index per grid node.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct DefaultTrackedTraits
{
	/// Single index stored in each grid node.
	using Leaf = T;
	/// Dimension of grid.
	static constexpr Dim t_dims = D;
	/// Number of lists tracking grid nodes.
	static constexpr TupleIdx t_num_lists = N;
};

/**
 * Traits for Tracked::LazySingle.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Tracked::LazyMultiListSingleIdxByValue<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	/// Type of lookup grid for tracking active positions.
	using Lookup = Lookup::LazyMultiListSingleIdx<D, N>;
};

/**
 * Traits for Tracked::Multi.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Tracked::MultiListMultiIdxByRef<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	/// Type of lookup grid for tracking active positions.
	using Lookup = Lookup::MultiListMultiIdx<D, N>;
};

} // Impl
} // Felt

#endif /* FELT_PUBLIC_TRACKED_HPP_ */
