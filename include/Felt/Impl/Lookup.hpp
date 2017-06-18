#ifndef FELT_PUBLIC_LOOKUP_HPP_
#define FELT_PUBLIC_LOOKUP_HPP_

#include <type_traits>
#include <Felt/public/Util.hpp>
#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Mixin/Grid.hpp>
#include <Felt/Impl/Mixin/Lookup.hpp>

namespace Felt
{
namespace Impl
{
namespace Lookup
{

template <UINT D>
class Simple :
	FELT_MIXINS(
		(Simple<D>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Lookup::Simple),
		(Grid::Index)
	)
private:
	using ThisType = Simple<D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Simple<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Simple(const VecDi& size, const VecDi& offset) :
		DataImpl{size, offset, NULL_IDX}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using LookupImpl::add;
	using LookupImpl::is_active;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class Single :
	FELT_MIXINS(
		(Single<D,N>),
		(Grid::Accessor::ByValue)(Grid::Accessor::Ref)(Grid::Activator)(Grid::Data)(Lookup::Single),
		(Grid::Index)
	)
private:
	using ThisType = Single<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByValue<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Single(const VecDi& size, const VecDi& offset) :
		DataImpl{size, offset, NULL_IDX}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using LookupImpl::add;
	using LookupImpl::is_active;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class LazySingle :
	FELT_MIXINS(
		(LazySingle<D,N>),
		(Grid::Accessor::LazyByValue)(Grid::Accessor::Ref)(Grid::Data)(Lookup::Activator)
		(Lookup::Single),
		(Grid::Accessor::ByValue)(Grid::Activator)(Grid::Index)
	)
private:
	using ThisType = LazySingle<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::LazyByValue<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using ActivatorImpl = Impl::Mixin::Lookup::Activator<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Single<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	using LookupImpl::PosArray;

public:
	LazySingle(const VecDi& size, const VecDi& offset) :
		DataImpl{size, offset, NULL_IDX}
	{}

	using AccessorImpl::get;
	using ActivatorImpl::activate;
	using ActivatorImpl::deactivate;
	using DataImpl::data;
	using LookupImpl::add;
	using LookupImpl::is_active;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class Multi :
	FELT_MIXINS(
		(Multi<D, N>),
		(Grid::Accessor::ByRef)(Grid::Activator)(Grid::Data)(Lookup::Multi),
		(Grid::Index)
	)
private:
	using ThisType = Multi<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Mixin::Grid::Accessor::ByRef<ThisType>;
	using DataImpl = Impl::Mixin::Grid::Data<ThisType>;
	using LookupImpl = Impl::Mixin::Lookup::Multi<ThisType>;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Multi(const VecDi& size, const VecDi& offset) :
		DataImpl{size, offset, LookupImpl::NULL_IDX_TUPLE}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using LookupImpl::add;
	using LookupImpl::is_active;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};

} // Lookup

/**
 * Traits for SimpleLookupGrid.
 *
 * @tparam D number of dimensions of the grid.
 */
template <UINT D>
struct Traits< Lookup::Simple<D> >
{
	using LeafType = UINT;
	static constexpr UINT Dims = D;
};

/**
 * Traits for SingleLookupGrid.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::Single<D, N> >
{
	using LeafType = UINT;
	static constexpr UINT Dims = D;
	static constexpr UINT NumLists = N;
};

/**
 * Traits for LazySingleLookupGrid.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::LazySingle<D, N> >
{
	using LeafType = UINT;
	static constexpr UINT Dims = D;
	static constexpr UINT NumLists = N;
};

/**
 * Traits for MultiLookupGrid.
 *
 * @tparam D number of dimensions of the grid.
 * @tparam N number of tracking lists.
 */
template <UINT D, UINT N>
struct Traits< Lookup::Multi<D, N> >
{
	using LeafType = VecDu<D>;
	static constexpr UINT Dims = D;
	static constexpr UINT NumLists = N;
};


} // Impl
} // Felt

#endif
