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

/**
 * A simple grid storing arbitrary types with by-value getter/setter.
 */
template <typename T, Dim D>
class Simple :
	FELT_MIXINS(
		(Simple<T, D>),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Data)(Grid::Size),
		(Grid::Index)
	)
private:
	using This = Simple<T, D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	Simple(const VecDi& size_, const VecDi& offset_, const Leaf background_) :
		ActivateImpl{background_}, SizeImpl{size_, offset_}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using AccessImpl::set;
	using DataImpl::data;
	using SizeImpl::inside;
	using SizeImpl::offset;
	using SizeImpl::size;
};


/**
 * A simple grid storing arbitrary types with by-value getter/setter.
 *
 * Includes Eigen::Array access to underlying data.
 */
template <typename T, Dim D>
class Snapshot :
	FELT_MIXINS(
		(Snapshot<T, D>),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Data)(Grid::Size)(Numeric::Snapshot),
		(Grid::Index)
	)
private:
	using This = Snapshot<T, D>;
	using Traits = Impl::Traits<This>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<This>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<This>;
	using DataImpl = Impl::Mixin::Grid::Data<This>;
	using SizeImpl = Impl::Mixin::Grid::Size<This>;
	using SnapshotImpl = Impl::Mixin::Numeric::Snapshot<This>;

	using VecDi = Felt::VecDi<Traits::t_dims>;
	using Leaf = typename Traits::Leaf;

public:
	using SnapshotImpl::VArrayData;

	Snapshot(const VecDi& size_, const VecDi& offset_, const Leaf background_) :
		 ActivateImpl{background_}, SizeImpl{size_, offset_}
	{
		this->activate();
	}

	using AccessImpl::get;
	using AccessImpl::index;
	using AccessImpl::set;
	using DataImpl::data;
	using SizeImpl::offset;
	using SizeImpl::size;
	using SnapshotImpl::array;
};
} // Grid
} // Impl
} // Felt


namespace Felt
{
namespace Impl
{
/**
 * Traits for Grid::Simple.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Grid::Simple<T, D> >
{
	using Leaf = T;
	static constexpr Dim t_dims = D;
};

/**
 * Traits for Grid::Snapshot.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Grid::Snapshot<T, D> >
{
	using Leaf = T;
	static constexpr Dim t_dims = D;
};

} // Impl
} // Felt

#endif
