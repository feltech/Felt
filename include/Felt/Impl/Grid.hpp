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
template <typename T, Dim D>
class Simple :
	FELT_MIXINS(
		(Simple<T, D>),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Data)(Grid::Size),
		(Grid::Index)
	)
private:
	using ThisType = Simple<T, D>;
	using TraitsType = Impl::Traits<ThisType>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<ThisType>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;

	using VecDi = Felt::VecDi<TraitsType::t_dims>;
	using LeafType = typename TraitsType::LeafType;

public:
	Simple(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
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


template <typename T, Dim D>
class Snapshot :
	FELT_MIXINS(
		(Snapshot<T, D>),
		(Grid::Access::ByValue)(Grid::Activate)(Grid::Data)(Grid::Size)(Numeric::Snapshot),
		(Grid::Index)
	)
private:
	using ThisType = Snapshot<T, D>;
	using TraitsType = Impl::Traits<ThisType>;

	using AccessImpl = Impl::Mixin::Grid::Access::ByValue<ThisType>;
	using ActivateImpl = Impl::Mixin::Grid::Activate<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using SizeImpl = Impl::Mixin::Grid::Size<ThisType>;
	using SnapshotImpl = Impl::Mixin::Numeric::Snapshot<ThisType>;

	using VecDi = Felt::VecDi<TraitsType::t_dims>;
	using LeafType = typename TraitsType::LeafType;

public:
	using SnapshotImpl::VArrayData;

	Snapshot(const VecDi& size_, const VecDi& offset_, const LeafType background_) :
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


/**
 * Traits for Grid::Simple.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, Dim D>
struct Traits< Grid::Simple<T, D> >
{
	using LeafType = T;
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
	using LeafType = T;
	static constexpr Dim t_dims = D;
};

} // Impl
} // Felt

#endif
