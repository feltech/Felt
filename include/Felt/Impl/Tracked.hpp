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
class LazySingleByValue :
	FELT_MIXINS(
		(LazySingleByValue<T, D, N>),
		(Grid::Access::LazyByValue)(Grid::Data)(Tracked::Activate)(Tracked::Multi::ByValue)
		(Tracked::Multi::LookupInterface)(Tracked::Multi::Reset)(Tracked::Resize),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Index)(Grid::Resize)(Grid::Size)
	)
private:
	using ThisType = LazySingleByValue<T, D, N>;
	using TraitsType = Impl::Traits<ThisType>;

	using AccessImpl = Impl::Mixin::Grid::Access::LazyByValue<ThisType>;
	using ActivateImpl = Impl::Mixin::Tracked::Activate<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::Multi::LookupInterface<ThisType>;
	using ResetterImpl = Impl::Mixin::Tracked::Multi::Reset<ThisType>;
	using SizeImpl = Impl::Mixin::Tracked::Resize<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Multi::ByValue<ThisType>;

	using VecDi = Felt::VecDi<TraitsType::t_dims>;
	using LeafType = typename TraitsType::LeafType;

public:
	using LookupType = typename TraitsType::LookupType;

	LazySingleByValue(const LeafType background_) :
		ActivateImpl{background_}, LookupInterfaceImpl{LookupType{}}
	{}

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
class SingleByRef :
	FELT_MIXINS(
		(SingleByRef<T, D>),
		(Grid::Access::ByRef)(Grid::Activate)(Grid::Data)(Grid::Size)(Tracked::Single::ByRef)
		(Tracked::LookupInterface),
		(Grid::Index)
	)
private:
	using ThisType = SingleByRef<T, D>;
	using TraitsType = Impl::Traits<ThisType>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByRef<ThisType>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::LookupInterface<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Single::ByRef<ThisType>;

	using VecDi = Felt::VecDi<TraitsType::t_dims>;
	using LeafType = typename TraitsType::LeafType;
	using LookupType = typename TraitsType::LookupType;

public:
	static constexpr TupleIdx t_num_lists = TraitsType::t_num_lists;

	SingleByRef(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{LookupType{size_, offset_}}
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


template <typename T, Dim D, TupleIdx N>
class MultiByRef :
	FELT_MIXINS(
		(MultiByRef<T, D, N>),
		(Grid::Access::ByRef)(Grid::Activate)(Grid::Data)(Grid::Size)(Tracked::Multi::ByRef)
		(Tracked::Multi::LookupInterface),
		(Grid::Index)
	)
private:
	using ThisType = MultiByRef<T, D, N>;
	using TraitsType = Impl::Traits<ThisType>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByRef<ThisType>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupInterfaceImpl = Impl::Mixin::Tracked::Multi::LookupInterface<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using TrackedImpl = Impl::Mixin::Tracked::Multi::ByRef<ThisType>;

	using VecDi = Felt::VecDi<TraitsType::t_dims>;
	using LeafType = typename TraitsType::LeafType;
	using LookupType = typename TraitsType::LookupType;

public:
	static constexpr TupleIdx t_num_lists = TraitsType::t_num_lists;

	MultiByRef(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_},
		LookupInterfaceImpl{LookupType{size_, offset_}}
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


} // Tracked.

/**
 * Traits for Tracked::SingleByRef.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Tracked::SingleByRef<T, D> >
{
	/// Single index stored in each grid node.
	using LeafType = T;
	/// Dimension of grid.
	static constexpr Dim t_dims = D;
	/// Type of lookup grid for tracking active positions.
	using LookupType = Lookup::Simple<D>;
};

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
