#ifndef FELT_PUBLIC_GRID_HPP_
#define FELT_PUBLIC_GRID_HPP_

#include <type_traits>
#include <Felt/Impl/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/GridMixin.hpp>
#include <Felt/Impl/Mixin/NumericMixin.hpp>

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
	Simple(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		SizeImpl{size_, offset_}, ActivatorImpl{background_}
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


template <typename T, UINT D>
class Snapshot :
	FELT_MIXINS(
		(Snapshot<T, D>),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Data)(Grid::Size)(Numeric::Snapshot),
		(Grid::Index)
	)
private:
	using ThisType = Snapshot<T, D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using ActivatorImpl = Impl::Mixin::Grid::Activator<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using SnapthotImpl = Impl::Mixin::Numeric::Snapshot<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	using SnapthotImpl::VArrayData;

	Snapshot(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
		SizeImpl{size_, offset_}, ActivatorImpl{background_}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using AccessorImpl::set;
	using DataImpl::data;
	using SizeImpl::offset;
	using SizeImpl::size;
	using SnapthotImpl::vdata;
};
} // Grid


/**
 * Traits for Grid::Simple.
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

/**
 * Traits for Grid::Snapshot.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, UINT D>
struct Traits< Grid::Snapshot<T, D> >
{
	using LeafType = T;
	static constexpr UINT Dims = D;
};

} // Impl
} // Felt

#endif
