#ifndef FELT_PUBLIC_GRID_HPP_
#define FELT_PUBLIC_GRID_HPP_

#include <type_traits>
#include <Felt/Impl/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/Grid.hpp>
#include <Felt/Impl/Mixin/Lookup.hpp>

namespace Felt
{
namespace Impl
{
namespace Grid
{
template <typename T, UINT D>
class Simple :
	FELT_MIXINS(
		(Simple<T, D>),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Data)(Grid::Size),
		(Grid::Index)
	)
private:
	using ThisType = Simple<T, D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Simple(const VecDi& size, const VecDi& offset, const LeafType background) :
		SizeImpl{size, offset}, ActivatorImpl{background}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using AccessorImpl::set;
	using DataImpl::data;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::size;
};
} // Grid
/**
 * Traits for Grid.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, UINT D>
struct Traits< Grid::Simple<T, D> >
{
	using LeafType = T;
	static constexpr UINT Dims = D;
};

} // Impl
} // Felt

#endif
