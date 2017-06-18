#ifndef FELT_PUBLIC_TRACKED_HPP_
#define FELT_PUBLIC_TRACKED_HPP_

#include <Felt/Impl/Lookup.hpp>
#include <Felt/Impl/Mixin/Grid.hpp>
#include <Felt/Impl/Mixin/Tracked.hpp>

namespace Felt
{
namespace Tracked
{

template <typename T, UINT D, UINT N>
class LazySingle :
	FELT_MIXINS(
		(LazySingle<T, D, N>),
		(Grid::Accessor::LazyByValue)(Grid::Data)(Tracked::Activator)(Tracked::LazySingleByValue),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)
	)
private:
	using ThisType = LazySingle<T, D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Tracked::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using TrackerImpl = Impl::Mixin::Tracked::LazySingleByValue<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	LazySingle(const VecDi& size, const VecDi& offset, const LeafType background) :
		DataImpl{size, offset, background},
		TrackerImpl{size, offset}
	{}

	using AccessorImpl::get;
	using AccessorImpl::set;
	using ActivatorImpl::activate;
	using ActivatorImpl::deactivate;
	using DataImpl::data;
	using TrackerImpl::add;
	using TrackerImpl::is_active;
	using TrackerImpl::list;
	using TrackerImpl::lookup;
	using TrackerImpl::remove;
	using TrackerImpl::reset;
};

} // Tracked.

namespace Impl
{

/**
 * Traits for Grid.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <typename T, UINT D, UINT N>
struct Traits< Felt::Tracked::LazySingle<T, D, N> >
{
	using LeafType = T;
	static constexpr UINT Dims = D;
	static constexpr UINT NumLists = N;
};

} // Impl
} // Felt



#endif /* FELT_PUBLIC_TRACKED_HPP_ */
