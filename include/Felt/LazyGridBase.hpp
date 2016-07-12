#ifndef INCLUDE_FELT_LAZYGRIDBASE_HPP_
#define INCLUDE_FELT_LAZYGRIDBASE_HPP_

#include "EagerGridBase.hpp"

namespace felt
{

/**
 * A lazy-loaded D-dimensional grid for storing values of type T.
 *
 * @tparam Derived the CRTP derived class.
 */
template <class Derived>
class LazyGridBase : public EagerGridBase<Derived>
{
public:
	using ThisType = LazyGridBase<Derived>;
	using DerivedType = typename GridTraits<Derived>::ThisType;
	using Base = EagerGridBase<Derived>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::VecDf;
	using typename Base::LeafType;
	using Base::EagerGridBase;
public:

	/**
	 * Default constructor
	 */
	LazyGridBase () : Base()
	{}

	/**
	 * Get whether this grid is active or not.
	 *
	 * An inactive grid stores no data and always returns the background value.
	 *
	 * @return boolean of the current active state of this grid.
	 */
	const bool is_active() const
	{
		return this->m_data.size() > 0;
	}

	using Base::size;

	/**
	 * Set the dimensions of the grid, but do not alter the data array.
	 *
	 * @param size_ new size of the grid.
	 */
	void size (const VecDu& size_)
	{
		this->m_size = size_;
	}

	/**
	 * Create the internal data array and fill with background value.
	 *
	 * @snippet test_Grid.cpp LazyGrid activation
	 */
	void activate()
	{
		Base::activate();
	}

	/**
	 * Destroy the internal data array.
	 *
	 * @snippet test_Grid.cpp LazyGrid deactivation
	 */
	void deactivate()
	{
		this->m_data.clear();
		this->m_data.shrink_to_fit();
	}

	/**
	 * Get value at position in grid, or return background value if grid is inactive.
	 *
	 * @snippet test_Grid.cpp LazyGrid activation
	 * @snippet test_Grid.cpp LazyGrid deactivation
	 *
	 * @param pos_ position in grid to query.
	 * @return reference to value represented at given position.
	 */
	LeafType& get (const VecDi& pos_)
	{
		if (is_active())
			return this->get_internal(pos_);
		return this->m_background;
	}

	/**
	 * @copydoc LazyGrid::get(const VecDi&)
	 * @brief Const version.
	 */
	const LeafType& get (const VecDi& pos_) const
	{
		if (is_active())
			return this->get_internal(pos_);
		return this->m_background;
	}
};


/**
 * Traits for LazyGridBase.
 *
 * Just forward the traits defined for LazyGridBase subclasses.
 */
template <class Derived>
struct GridTraits<LazyGridBase<Derived> > : GridTraits<Derived>
{};

}

#endif /* INCLUDE_FELT_LAZYGRIDBASE_HPP_ */
