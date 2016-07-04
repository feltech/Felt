#ifndef Grid_hpp
#define Grid_hpp

#include "StaticGridBase.hpp"
#include "LazyGridBase.hpp"

namespace felt
{

/**
 * Placeholder GridBase to be specialised based on IsLazy template parameter.
 *
 * @tparam Derived the CRTP derived class.
 * @tparam IsLazy true if base class is a LazyGridBase, false if base class is a StaticGridBase.
 */
template <class Derived, bool IsLazy=false>
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
class GridBase<Derived, false> : public StaticGridBase<Derived>
{
public:
	using StaticGridBase<Derived>::StaticGridBase;
};


/**
 * Specialisation of GridBase that inherits from LazyGridBase.
 *
 * Used to provide templated choice of base class.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class GridBase<Derived, true> : public LazyGridBase<Derived>
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
class Grid : public GridBase<Grid<T, D>, false>
{
public:
	using ThisType = Grid<T, D>;
	using Base = GridBase<ThisType, false>;
	using Base::GridBase;
};


/**
 * A lazy-loaded D-dimensional grid for storing values of type T.
 *
 * @tparam T the type of data to store in the grid.
 * @tparam D the dimensions of the grid.
 */
template <typename T, UINT D>
class LazyGrid : public GridBase<LazyGrid<T, D>, true>
{
public:
	using Base = GridBase<LazyGrid<T, D>, true>;
	using Base::GridBase;
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


