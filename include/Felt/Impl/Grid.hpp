#ifndef FELT_IMPL_GRID_HPP
#define FELT_IMPL_GRID_HPP

#include <vector>
#include <Felt/Impl/Base.hpp>
#include <Felt/public/Util.hpp>

namespace Felt
{
namespace Impl
{
namespace Grid
{


template <class Derived>
class Data : private Base<Derived>
{
protected:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;

	/**
	 * D-dimensional signed integer vector.
	 */
	using VecDi = Felt::VecDi<Dims>;
	/**
	 * Dynamic 1D vector (a resizeable array of data) for storage of grid data.
	 */
	using ArrayData = std::vector<LeafType>;

private:
	/**
	 * The dimensions (size) of the grid.
	 */
	VecDi	m_size;
	/**
	 * The translational offset of the grid's zero coordinate.
	 */
	VecDi	m_offset;

protected:
	/// The actual grid data store.
	ArrayData	m_data;
	/// Default/initial value of grid nodes.
	LeafType	m_background;

protected:

	Data(const VecDi& size, const VecDi& offset, const LeafType background) :
		m_size{size}, m_offset{offset}, m_background{background}
	{}

	ArrayData& data ()
	{
		return m_data;
	}
	const ArrayData& data () const
	{
		return m_data;
	}
	const VecDi& size () const
	{
		return m_size;
	}
	const VecDi& offset () const
	{
		return m_offset;
	}

	/**
	 * Construct the internal data array, initialising nodes to the background value.
	 */
	void activate()
	{
		INT arr_size = m_size(0);
		for (INT i = 1; i < m_size.size(); i++)
			arr_size *= m_size(i);
		m_data.resize(arr_size, m_background);
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
	 * Test if a position is inside the grid bounds.
	 *
	 * @tparam PosType the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename PosType>
	bool inside (const Felt::VecDT<PosType, Dims>& pos_) const
	{
		return inside(pos_, cself->offset(), cself->offset() + cself->size());
	}

	/**
	 * Test if a position is inside given bounds.
	 *
	 * @tparam PosType the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @param pos_min_ minimum allowed position.
	 * @param pos_max_ one more than the maximum allowed position.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename ElemType>
	static bool inside (
		const Felt::VecDT<ElemType, Dims>& pos_, const VecDi& pos_min_, const VecDi& pos_max_
	) {
		for (INT i = 0; i < pos_.size(); i++)
		{
			if (pos_(i) >= static_cast<ElemType>(pos_max_(i)))
				return false;
			if (pos_(i) < static_cast<ElemType>(pos_min_(i)))
				return false;
		}
		return true;
	}

	/**
	 * Check if given position is within the grid and raise a domain_error if not.
	 *
	 * @param pos_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_bounds (const VecDi& pos_, std::string title_) const
	{
		if (!cself->inside(pos_))
		{
			const VecDi& pos_min = cself->offset();
			const VecDi& pos_max = (
				cself->size() + pos_min - VecDi::Constant(1)
			);
			std::stringstream err;
			err << title_ << format(pos_.transpose())
				<< " is outside grid "
				<< format(pos_min) << "-" << format(pos_max)
				<< std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}
	}
};


template <class Derived>
class Index : private Base<Derived>
{
protected:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/**
	 * D-dimensional signed integer vector.
	 */
	using VecDi = Felt::VecDi <Dims>;

protected:

	/**
	 * Get index in data array of position vector.
	 *
	 * The grid is packed in a 1D array, so this method is required to
	 * get the index in that array of the D-dimensional position.
	 *
	 * @snippet test_Grid.cpp Position index
	 *
	 * @param pos_ position in grid to query.
	 * @return index in internal data array of this grid position.
	 */
	INT index (const VecDi& pos_) const
	{
		return index(pos_, cself->size(), cself->offset());
	}

	/**
	 * Get position of index.
	 *
	 * Given an index in the 1D grid data array, calculate the position vector that it pertains to.
	 *
	 * @snippet test_Grid.cpp Position index
	 *
	 * @param idx_ index in internal data array to query.
	 * @return the position in the grid represented in the data array at given index.
	 */
	VecDi index (const UINT idx_) const
	{
		return index(idx_, cself->size(), cself->offset());
	}


	/**
	 * Get index in data array of position vector.
	 *
	 * The grid is packed in a 1D array, so this method is required to
	 * get the index in that array of the D-dimensional position.
	 *
	 * @snippet test_Grid.cpp Position index
	 *
	 * @param pos_ position in grid.
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @return index in data array of pos in grid of given size and offset.
	 */
	static INT index (
		const VecDi& pos_, const VecDi& size_, const VecDi& offset_ = VecDi::Constant(0)
	) {
		INT idx = 0;
		for (INT i = 0; i < size_.size(); i++)
		{
			INT u_pos = pos_(i) - offset_(i);

			for (INT j = i+1; j < size_.size(); j++)
				u_pos *= size_(j);

			idx += u_pos;
		}
		return idx;
	}

	/**
	 * Get position of index.
	 *
	 * Given an index and the dimensions and offset of a grid, calculate
	 * the position vector that the index pertains to in a representative
	 * 1D array.
	 *
	 * @snippet test_Grid.cpp Position index
	 *
	 * @param idx_ index in to query.
	 * @param size_ size of grid.
	 * @param offset_ spatial offset of grid.
	 * @return position that the given index would represent in a grid of given size and offset.
	 */
	static VecDi index (UINT idx_, const VecDi& size_, const VecDi& offset_ = VecDi::Zero())
	{
/*
Eg. 2D: row major order (3x4=12): (x,y)[idx] =>
	(0,0)[0], (0,1)[1], (0,2)[2],  (0,3)[3]
	(1,0)[4], (1,1)[5], (1,2)[6],  (1,3)[7]
	(2,0)[8], (2,1)[9], (2,2)[10], (2,3)[11]

E.g. 3D:
z = idx % Dz
y = (idx/Dz) % Dy
x = (idx/Dz)/Dy % Dx
*/
		VecDi pos;

		for (INT axis = size_.size()-1; axis >= 0; axis--)
		{
			pos(axis) = idx_ % size_(axis) + offset_(axis);
			idx_ /= size_(axis);
		}

		return pos;
	}
};



namespace Accessor
{

template <class Derived>
class Ref : private Base<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/**
	 * Get a reference to the underlying data in the grid.
	 *
	 * @param pos_ position to query.
	 * @return reference to grid value.
	 */
	LeafType& ref(const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = cself->index(pos_);
		return nself->data()[idx];
	}

	/**
	 * @copydoc LeafType& ref(const VecDi&)
	 */
	const LeafType& ref(const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = cself->index(pos_);
		return cself->data()[idx];
	}
};


template <class Derived>
class ByValue : protected Index<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// Base class, for subclass friending.
	using Base = Index<Derived>;

protected:

	/**
	 * Get a copy of the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = cself->index(pos_);
		return cself->data()[idx];
	}


	/**
	 * Set the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	void set (const VecDi& pos_, LeafType val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = cself->index(pos_);
		nself->data()[idx] = val_;
	}
};


template <class Derived>
class ByRef : protected Index<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// Base class, for subclass friending.
	using Base = Index<Derived>;

protected:

	/**
	 * Get a reference to the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType& get (const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = cself->index(pos_);
		return nself->data()[idx];
	}

	/**
	 * Get a const reference to the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	const LeafType& get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = cself->index(pos_);
		return cself->data()[idx];
	}
};


template <class Derived>
class LazyByValue : protected ByValue<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// Base class, for subclass friending.
	using Base = ByValue<Derived>;

protected:

	/**
	 * Get a copy of the value stored in the grid or background value if inactive.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		cself->assert_pos_bounds(pos_, "get: ");
		#endif
		if (cself->data().size())
		{
			const UINT idx = cself->index(pos_);
			return cself->data()[idx];
		}
		else
		{
			return cself->m_background;
		}
	}
};

} // Accessor


} // Grid
} // Impl
} // Felt

#endif // FELT_IMPL_GRID_HPP
