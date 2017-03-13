#include <type_traits>
#include <Felt/public/Util.hpp>
#include <Felt/Impl/Base.hpp>
#include <Felt/Impl/Grid.hpp>
#include <Felt/Impl/Lookup.hpp>

namespace Felt
{

template <typename T, UINT D>
class Grid :
		private Impl::Grid::Accessor::ByValue< Grid<T, D> >,
		private Impl::Grid::Data< Grid<T, D> >
{
private:
	using ThisType = Grid<T, D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Grid::Accessor::ByValue<ThisType>;
	using DataImpl = Impl::Grid::Data<ThisType>;

	friend AccessorImpl;
	friend DataImpl;
	friend typename AccessorImpl::Base;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	Grid(const VecDi& size, const VecDi& offset, const LeafType background) :
		DataImpl{size, offset, background}
	{
		this->activate();
	}

	using AccessorImpl::get;
	using AccessorImpl::index;
	using AccessorImpl::set;
	using DataImpl::data;
	using DataImpl::inside;
	using DataImpl::offset;
	using DataImpl::size;
};


template <UINT D>
class SimpleLookupGrid :
		private Impl::Grid::Accessor::ByValue< SimpleLookupGrid<D> >,
		private Impl::Grid::Accessor::Ref< SimpleLookupGrid<D> >,
		private Impl::Grid::Data< SimpleLookupGrid<D> >,
		private Impl::Lookup::Simple< SimpleLookupGrid<D> >
{
private:
	using ThisType = SimpleLookupGrid<D>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Grid::Accessor::ByValue<ThisType>;
	using AccessorRefImpl = Impl::Grid::Accessor::Ref<ThisType>;
	using DataImpl = Impl::Grid::Data<ThisType>;
	using LookupImpl = Impl::Lookup::Simple<ThisType>;

	friend AccessorImpl;
	friend AccessorRefImpl;
	friend DataImpl;
	friend LookupImpl;
	friend typename AccessorImpl::Base;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	SimpleLookupGrid(const VecDi& size, const VecDi& offset) :
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
class SingleLookupGrid :
		private Impl::Grid::Accessor::ByValue< SingleLookupGrid<D, N> >,
		private Impl::Grid::Accessor::Ref< SingleLookupGrid<D, N> >,
		private Impl::Grid::Data< SingleLookupGrid<D, N> >,
		private Impl::Lookup::Single< SingleLookupGrid<D, N> >
{
private:
	using ThisType = SingleLookupGrid<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Grid::Accessor::ByValue<ThisType>;
	using AccessorRefImpl = Impl::Grid::Accessor::Ref<ThisType>;
	using DataImpl = Impl::Grid::Data<ThisType>;
	using LookupImpl = Impl::Lookup::Single<ThisType>;

	friend AccessorImpl;
	friend AccessorRefImpl;
	friend DataImpl;
	friend LookupImpl;
	friend typename AccessorImpl::Base;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	SingleLookupGrid(const VecDi& size, const VecDi& offset) :
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
class LazySingleLookupGrid :
		private Impl::Grid::Accessor::LazyByValue< LazySingleLookupGrid<D, N> >,
		private Impl::Grid::Accessor::Ref< LazySingleLookupGrid<D, N> >,
		private Impl::Grid::Data< LazySingleLookupGrid<D, N> >,
		private Impl::Lookup::Deactivator< LazySingleLookupGrid<D, N> >,
		private Impl::Lookup::Single< LazySingleLookupGrid<D, N> >
{
private:
	using ThisType = LazySingleLookupGrid<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Grid::Accessor::LazyByValue<ThisType>;
	using AccessorRefImpl = Impl::Grid::Accessor::Ref<ThisType>;
	using DataImpl = Impl::Grid::Data<ThisType>;
	using DeactivatorImpl = Impl::Lookup::Deactivator<ThisType>;
	using LookupImpl = Impl::Lookup::Single<ThisType>;

	friend AccessorImpl;
	friend AccessorRefImpl;
	friend DataImpl;
	friend DeactivatorImpl;
	friend LookupImpl;
	friend typename AccessorImpl::Base;
	friend typename AccessorImpl::Base::Base;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	LazySingleLookupGrid(const VecDi& size, const VecDi& offset) :
		DataImpl{size, offset, NULL_IDX}
	{}

	using AccessorImpl::get;
	using DataImpl::activate;
	using DataImpl::data;
	using DeactivatorImpl::deactivate;
	using LookupImpl::add;
	using LookupImpl::is_active;
	using LookupImpl::list;
	using LookupImpl::remove;
	using LookupImpl::reset;
};


template <UINT D, UINT N>
class MultiLookupGrid :
		private Impl::Grid::Accessor::ByRef< MultiLookupGrid<D, N> >,
		private Impl::Grid::Data< MultiLookupGrid<D, N> >,
		private Impl::Lookup::Multi< MultiLookupGrid<D, N> >
{
private:
	using ThisType = MultiLookupGrid<D, N>;
	using ThisTraits = Impl::Traits<ThisType>;

	using AccessorImpl = Impl::Grid::Accessor::ByRef<ThisType>;
	using DataImpl = Impl::Grid::Data<ThisType>;
	using LookupImpl = Impl::Lookup::Multi<ThisType>;

	friend AccessorImpl;
	friend DataImpl;
	friend LookupImpl;
	friend typename AccessorImpl::Base;

	using VecDi = Felt::VecDi<ThisTraits::Dims>;
	using LeafType = typename ThisTraits::LeafType;

public:
	MultiLookupGrid(const VecDi& size, const VecDi& offset) :
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



namespace Impl
{
/**
 * Traits for Grid.
 *
 * @tparam T type of data to store in the grid.
 * @tparam D number of dimensions of the grid.
 */
template <typename T, UINT D>
struct Traits< Felt::Grid<T, D> >
{
	using LeafType = T;
	static constexpr UINT Dims = D;
};

/**
 * Traits for SimpleLookupGrid.
 *
 * @tparam D number of dimensions of the grid.
 */
template <UINT D>
struct Traits< Felt::SimpleLookupGrid<D> >
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
struct Traits< Felt::SingleLookupGrid<D, N> >
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
struct Traits< Felt::LazySingleLookupGrid<D, N> >
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
struct Traits< Felt::MultiLookupGrid<D, N> >
{
	using LeafType = VecDu<D>;
	static constexpr UINT Dims = D;
	static constexpr UINT NumLists = N;
};
} // Impl

} // Felt
