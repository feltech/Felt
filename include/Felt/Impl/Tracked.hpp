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

template <typename T, Dim D>
class LazySimpleByRef :
	FELT_MIXINS(
		(LazySimpleByRef<T, D>),
		(Grid::Accessor::ByRef)(Grid::Data)(Tracked::Activator)(Tracked::Single::ByRef)
		(Tracked::LookupInterface)(Tracked::Single::Resetter)(Tracked::Resize),
		(Grid::Activator)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using ThisType = LazySimpleByRef<T, D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using ActivatorImpl = Impl::Mixin::Tracked::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using ResetterImpl = Impl::Mixin::Tracked::Single::Resetter<ThisType>;
	using SizeImpl = Impl::Mixin::Tracked::Resize<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Single::ByRef<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;
	using LookupType = typename ThisTraits::LookupType;

public:
	LazySimpleByRef(const LeafType background_) :
		ActivatorImpl{background_}, LookupInterfaceImpl{LookupType{}}
	{}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using ActivatorImpl::activate;
	using ActivatorImpl::background;
	using ActivatorImpl::deactivate;
	using ActivatorImpl::is_active;
	using DataImpl::data;
	using LookupInterfaceImpl::lookup;
	using ResetterImpl::reset;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
	using TrackedImpl::track;
};


template <typename T, Dim D, TupleIdx N>
class LazySingleByValue :
	FELT_MIXINS(
		(LazySingleByValue<T, D, N>),
		(Grid::Accessor::LazyByValue)(Grid::Data)(Tracked::Activator)(Tracked::Multi::ByValue)
		(Tracked::Multi::LookupInterface)(Tracked::Multi::Resetter)(Tracked::Resize),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using ThisType = LazySingleByValue<T, D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Tracked::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::Multi::LookupInterface<ThisType>;
	using ResetterImpl = Impl::Mixin::Tracked::Multi::Resetter<ThisType>;
	using SizeImpl = Impl::Mixin::Tracked::Resize<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Multi::ByValue<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
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
	using LookupInterfaceImpl::list;
	using LookupInterfaceImpl::lookup;
	using ResetterImpl::reset;
	using SizeImpl::assert_pos_idx_bounds;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::resize;
	using SizeImpl::size;
	using TrackedImpl::track;
};


template <typename T, Dim D, TupleIdx N>
class MultiByRef :
	FELT_MIXINS(
		(MultiByRef<T, D, N>),
		(Grid::Accessor::ByRef)(Grid::Activator)(Grid::Data)(Grid::Size)(Tracked::Multi::ByRef)
		(Tracked::Multi::LookupInterface),
		(Grid::Index)
	)
private:
	using ThisType = MultiByRef<T, D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::Multi::LookupInterface<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Multi::ByRef<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::t_dims>;
	using LeafType = typename ThisTraits::LeafType;
	using LookupType = typename ThisTraits::LookupType;

public:
	static constexpr TupleIdx t_num_lists = ThisTraits::t_num_lists;

	MultiByRef(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		ActivatorImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{LookupType{size_, offset_}}
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
template <typename T, Dim D, TupleIdx N>
struct DefaultTrackedTraits
{
	/// Single index stored in each grid node.
	using LeafType = T;
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
template <typename T, Dim D>
struct Traits< Tracked::LazySimpleByRef<T, D> >
{
	/// Single index stored in each grid node.
	using LeafType = T;
	/// Dimension of grid.
	static constexpr Dim t_dims = D;
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::LazySimple<D>;
};

/**
 * Traits for Tracked::LazySingle.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Tracked::LazySingleByValue<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::LazySingle<D, N>;
};

/**
 * Traits for Tracked::Multi.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, Dim D, TupleIdx N>
struct Traits< Tracked::MultiByRef<T, D, N> > : public DefaultTrackedTraits<T, D, N>
{
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::Multi<D, N>;
};

} // Impl
} // Felt



#endif /* FELT_PUBLIC_TRACKED_HPP_ */
