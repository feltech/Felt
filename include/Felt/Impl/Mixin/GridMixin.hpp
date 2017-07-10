#ifndef FELT_IMPL_GRID_HPP
#define FELT_IMPL_GRID_HPP

#include <vector>

#include <Felt/Impl/Common.hpp>
#include <Felt/Impl/Util.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Grid
{

template <class Derived>
class Activator
{
private:
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;

	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// Default/initial value of grid nodes.
	LeafType	m_background;

protected:
	Activator(const LeafType background_) :
		m_background{background_}
	{}

	/**
	 * Get whether this grid has been activated (data allocated) or not.
	 *
	 * @return true if data allocated, false if not.
	 */
	bool is_active() const
	{
		return bool(pself->data().size());
	}

	/**
	 * Construct the internal data array, initialising nodes to the background value.
	 */
	void activate()
	{
		INT arr_size = pself->m_size(0);
		for (INT i = 1; i < pself->m_size.size(); i++)
			arr_size *= pself->m_size(i);
		pself->m_data.resize(arr_size, m_background);
	}
};


template <class Derived>
class Data
{
private:
	/// Type of data to store in grid nodes.
	using LeafType = typename Traits<Derived>::LeafType;
	/// Dynamic 1D vector (a resizeable array of data) for storage of grid data.
	using ArrayData = std::vector<LeafType>;

protected:
	ArrayData	m_data;

protected:

	ArrayData& data ()
	{
		return m_data;
	}

	const ArrayData& data () const
	{
		return m_data;
	}
};


template <class Derived>
class Index
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
	PosIdx index (const VecDi& pos_) const
	{
		return index(pos_, pself->size(), pself->offset());
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
	VecDi index (const PosIdx idx_) const
	{
		return index(idx_, pself->size(), pself->offset());
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
	static PosIdx index (
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
	static VecDi index (PosIdx idx_, const VecDi& size_, const VecDi& offset_ = VecDi::Zero())
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


template <class Derived>
class Size
{
private:
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	/// The translational offset of the grid's zero coordinate.
	VecDi	m_offset;

	///The dimensions (size) of the grid.
	VecDi	m_size;

protected:

	Size(const VecDi& size, const VecDi& offset) :
		m_size{size}, m_offset{offset}
	{}

	const VecDi& size () const
	{
		return m_size;
	}

	const VecDi& offset () const
	{
		return m_offset;
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
		return inside(pos_, pself->offset(), pself->offset() + pself->size());
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
	 * Check if given position's index is within the data array and raise a domain_error if not.
	 *
	 * @param pos_idx_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_idx_bounds (const PosIdx pos_idx_, std::string title_) const
	{
		assert_pos_idx_bounds(pself->index(pos_idx_), title_);
	}

	/**
	 * Check if given position's index is within the data array and raise a domain_error if not.
	 *
	 * @param pos_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_idx_bounds (const VecDi& pos_, std::string title_) const
	{
		const PosIdx pos_idx = pself->index(pos_);

		if (pos_idx > pself->data().size())
		{
			const VecDi& pos_min = pself->offset();
			const VecDi& pos_max = (
				pself->size() + pos_min - VecDi::Constant(1)
			);
			std::stringstream err;
			err << title_ << format(pos_.transpose())
				<< " data index " << pos_idx << " is less than data size " <<
				pself->data().size() << " for grid " <<
				format(pos_min) << "-" << format(pos_max) << std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}

		assert_pos_bounds(pos_, title_);
	}

	/**
	 * Check if given position is within the grid and raise a domain_error if not.
	 *
	 * @param pos_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_bounds (const PosIdx pos_idx_, std::string title_) const
	{
		assert_pos_bounds(pself->index(pos_idx_), title_);
	}

	/**
	 * Check if given position is within the grid and raise a domain_error if not.
	 *
	 * @param pos_ position in grid to query.
	 * @param title_ message to include in generated exception.
	 */
	void assert_pos_bounds (const VecDi& pos_, std::string title_) const
	{
		if (!pself->inside(pos_))
		{
			const VecDi& pos_min = pself->offset();
			const VecDi& pos_max = (
				pself->size() + pos_min - VecDi::Constant(1)
			);
			std::stringstream err;
			err << title_ << format(pos_.transpose()) << " is outside grid " <<
				format(pos_min) << "-" << format(pos_max) << std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}
	}
};


template <class Derived>
class Resize : protected Size<Derived>
{
private:
	using Base = Size<Derived>;
	/// Dimension of the grid.
	static const UINT Dims = Traits<Derived>::Dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<Dims>;

protected:
	Resize() : Base::Size(VecDi(), VecDi())
	{}

	void resize(const VecDi& size_, const VecDi& offset_)
	{
		this->m_size = size_;
		this->m_offset = offset_;
	}
};


namespace Accessor
{

template <class Derived>
class Ref
{
private:
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
	 * @param pos_idx_ data index of position to query.
	 * @return reference to grid value.
	 */
	LeafType& ref(const UINT pos_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pself->index(pos_idx_), "ref: ");
		#endif
		return pself->data()[pos_idx_];
	}

	/**
	 * @copydoc LeafType& ref(const UINT)
	 */
	const LeafType& ref(const UINT pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pself->index(pos_idx_), "ref: ");
		#endif
		return pself->data()[pos_idx_];
	}

	/**
	 * Get a reference to the underlying data in the grid.
	 *
	 * @param pos_ position to query.
	 * @return reference to grid value.
	 */
	LeafType& ref(const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = pself->index(pos_);
		return ref(idx);
	}

	/**
	 * @copydoc LeafType& ref(const VecDi&)
	 */
	const LeafType& ref(const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = pself->index(pos_);
		return ref(idx);
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

	/**
	 * Get a copy of the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = this->index(pos_);
		return get(idx);
	}

	/**
	 * Get a copy of the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return internally stored value at given grid position
	 */
	LeafType get (const UINT pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(this->index(pos_idx_), "get: ");
		#endif
		return pself->data()[pos_idx_];
	}

	/**
	 * Set the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @param val_ value to copy into grid at pos_.
	 */
	void set (const VecDi& pos_, LeafType val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "set: ");
		#endif
		const UINT idx = this->index(pos_);
		set(idx, val_);
	}

	/**
	 * Set the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @param val_ value to copy into grid at pos_.
	 */
	void set (const UINT pos_idx_, LeafType val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_idx_, "set: ");
		#endif
		pself->data()[pos_idx_] = val_;
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

	/**
	 * Get a reference to the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType& get (const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = this->index(pos_);
		return get(idx);
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
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const UINT idx = this->index(pos_);
		return get(idx);
	}

	/**
	 * Get a reference to the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return internally stored value at given grid position
	 */
	LeafType& get (const UINT pos_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(this->index(pos_idx_), "get: ");
		#endif
		return pself->data()[pos_idx_];
	}

	/**
	 * Get a const reference to the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return internally stored value at given grid position
	 */
	const LeafType& get (const UINT pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(this->index(pos_idx_), "get: ");
		#endif
		return pself->data()[pos_idx_];
	}

};


template <class Derived>
class LazyByValue : protected ByValue<Derived>
{
private:
	/// CRTP derived class.
	using DerivedType = Derived;
	/// Base class
	using BaseType = ByValue<Derived>;
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
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		if (pself->data().size())
		{
			return BaseType::get(pos_);
		}
		else
		{
			return pself->m_background;
		}
	}

	/**
	 * Get a copy of the value stored in the grid or background value if inactive.
	 *
	 * @param pos_idx_ index of position in grid to query.
	 * @return internally stored value at given grid position
	 */
	LeafType get (const UINT pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pself->index(pos_idx_), "get: ");
		#endif
		if (pself->data().size())
		{
			return BaseType::get(pos_idx_);
		}
		else
		{
			return pself->m_background;
		}
	}
};

} // Accessor


} // Grid
} // Mixin
} // Impl
} // Felt

#endif // FELT_IMPL_GRID_HPP