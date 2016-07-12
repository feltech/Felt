#ifndef Grid_hpp
#define Grid_hpp

#include "EagerGridBase.hpp"
#include "LazyGridBase.hpp"

namespace felt
{

enum Laziness
{
	EAGER, LAZY
};

/**
 * Placeholder GridBase to be specialised based on IsLazy template parameter.
 *
 * @tparam Derived the CRTP derived class.
 * @tparam Lazy EAGER or LAZY depending whether the grid is lazy activated/deactivated, or
 * always active.
 */
template <class Derived, Laziness IsLazy=Laziness::EAGER>
class GridBase
{};


/**
 * Specialisation of GridBase that inherits from StaticGridBase.
 *
 * Used to provide templated choice of base class.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class GridBase<Derived, Laziness::EAGER> : public EagerGridBase<Derived>
{
public:
	using EagerGridBase<Derived>::EagerGridBase;
};


/**
 * Specialisation of GridBase that inherits from LazyGridBase.
 *
 * Used to provide templated choice of base class.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class GridBase<Derived, Laziness::LAZY> : public LazyGridBase<Derived>
{
public:
	using LazyGridBase<Derived>::LazyGridBase;
};


/**
 * A standard D-dimensional grid for storing values of type T.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class Grid : public EagerGridBase< Grid<T, D> >
{
public:
	using ThisType = Grid<T, D>;
	using Base = EagerGridBase< Grid<T, D> >;
	using Base::EagerGridBase;
};


/**
 * A lazy-loaded D-dimensional grid for storing values of type T.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class LazyGrid : public LazyGridBase< LazyGrid<T, D> >
{
public:
	using Base = LazyGridBase< LazyGrid<T, D> >;
	using Base::LazyGridBase;
};


/**
 * Default traits for grids to derive from.clipse clip
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct DefaultGridTraits
{
	/// Dimensions of the grid, from template parameter.
	static const UINT Dims = D;
	/// The data type to store at leaf grid nodes, from template parameter.
	using LeafType = T;
};


/**
 * Traits for Grid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridTraits<Grid<T, D> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = Grid<T, D>;
};


/**
 * Traits for LazyGrid.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the number of dimensions of the grid.
 */
template <typename T, UINT D>
struct GridTraits<LazyGrid<T, D> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = LazyGrid<T, D>;
};

/** @} */ // End group Grid.

}// End namespace felt.
#endif


