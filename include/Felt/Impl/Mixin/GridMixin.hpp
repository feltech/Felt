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

template <class TDerived>
class Activate
{
private:
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits<TDerived>::Leaf;

	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/// Default/initial value of grid nodes.
	Leaf	m_background;

protected:
	Activate(const Leaf background_) :
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
	 * Get the background value used to initially fill the grid.
	 *
	 * @return backround value.
	 */
	Leaf background() const
	{
		return m_background;
	}

	/**
	 * Construct the internal data array, initialising nodes to the background value.
	 */
	void activate()
	{
		INT arr_size = pself->m_size(0);
		for (INT i = 1; i < pself->m_size.size(); i++)
			arr_size *= pself->m_size(i);
		pself->m_data.resize(PosIdx(arr_size), m_background);
	}

	/**
	 * Throw exception if grid is inactive
	 *
	 * @param title_ text to prefix exception message with.
	 */
	void assert_is_active(std::string title_) const
	{
		if (!is_active())
		{
			const VecDi& pos_min = pself->offset();
			const VecDi& pos_max = (
				pself->size() + pos_min - VecDi::Constant(1)
			);
			std::stringstream err;
			err << title_ << "inactive grid " <<
				format(pos_min) << "-" << format(pos_max) << std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}
	}
};


template <class TDerived>
class Data
{
private:
	using Traits = Impl::Traits<TDerived>;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits::Leaf;
	static constexpr Dim t_dims = Traits::t_dims;
	using VecDi = Felt::VecDi<t_dims>;
	using DataArray = Felt::DataArray<Leaf>;
protected:
	DataArray	m_data;

protected:

	DataArray& data ()
	{
		return m_data;
	}

	const DataArray& data () const
	{
		return m_data;
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
				<< " data index " << pos_idx << " is greater than data size " <<
				pself->data().size() << " for grid " <<
				format(pos_min) << "-" << format(pos_max) << std::endl;
			std::string err_str = err.str();
			throw std::domain_error(err_str);
		}

		pself->assert_pos_bounds(pos_, title_);
	}
};


template <class TDerived>
class Index
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/**
	 * D-dimensional signed integer vector.
	 */
	using VecDi = Felt::VecDi<t_dims>;

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
		return Felt::index<t_dims>(pos_, pself->size(), pself->offset());
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
		return Felt::index<t_dims>(idx_, pself->size(), pself->offset());
	}

};


template <class TDerived>
class Size
{
private:
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	///The dimensions (size) of the grid.
	VecDi	m_size;
	/// The translational offset of the grid's zero coordinate.
	VecDi	m_offset;
	/// Cache for use in `inside`.
	VecDi	m_offset_plus_size;

protected:

	Size(const VecDi& size_, const VecDi& offset_) :
		m_size{size_}, m_offset{offset_}, m_offset_plus_size{offset_ + size_}
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
	 * @tparam Pos the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename T>
	bool inside (const Felt::VecDT<T, t_dims>& pos_) const
	{
		return inside(pos_, m_offset, m_offset_plus_size);
	}

	/**
	 * Test if a position is inside given bounds.
	 *
	 * @tparam Pos the type of position vector (i.e. float vs. int).
	 * @param pos_ position in grid to query.
	 * @param pos_min_ minimum allowed position.
	 * @param pos_max_ one more than the maximum allowed position.
	 * @return true if position lies inside the grid, false otherwise.
	 */
	template <typename T>
	static bool inside (
		const Felt::VecDT<T, t_dims>& pos_, const VecDi& pos_min_, const VecDi& pos_max_
	) {
		for (Dim i = 0; i < pos_.size(); i++)
		{
			if (pos_(i) >= static_cast<T>(pos_max_(i)))
				return false;
			if (pos_(i) < static_cast<T>(pos_min_(i)))
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


template <class TDerived>
class Resize : protected Size<TDerived>
{
private:
	using Base = Size<TDerived>;
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	Resize() : Base::Size(VecDi(), VecDi())
	{}

	void resize(const VecDi& size_, const VecDi& offset_)
	{
		this->m_size = size_;
		this->m_offset = offset_;
		this->m_offset_plus_size = this->m_offset + size_;
		this->m_offset_plus_one = this->m_offset + VecDi::Constant(1);
		this->m_offset_plus_size_minus_two = this->m_offset_plus_size - VecDi::Constant(2);
	}
};


namespace Access
{

template <class TDerived>
class Ref
{
private:
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits<TDerived>::Leaf;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:
	/**
	 * Get a reference to the underlying data in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return reference to grid value.
	 */
	Leaf& ref(const PosIdx pos_idx_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pself->index(pos_idx_), "ref: ");
		#endif
		return pself->data()[pos_idx_];
	}

	/**
	 * @copydoc Leaf& ref(const UINT)
	 */
	const Leaf& ref(const PosIdx pos_idx_) const
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
	Leaf& ref(const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = pself->index(pos_);
		return ref(idx);
	}

	/**
	 * @copydoc Leaf& ref(const VecDi&)
	 */
	const Leaf& ref(const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "ref: ");
		#endif
		const UINT idx = pself->index(pos_);
		return ref(idx);
	}
};


template <class TDerived>
class ByValue : protected Index<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Dimension of the grid.
	static const Dim t_dims = Traits<TDerived>::t_dims;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits<TDerived>::Leaf;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:

	/**
	 * Get a copy of the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	Leaf get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const PosIdx idx = this->index(pos_);
		return get(idx);
	}

	/**
	 * Get a copy of the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return internally stored value at given grid position
	 */
	Leaf get (const PosIdx pos_idx_) const
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
	void set (const VecDi& pos_, Leaf val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "set: ");
		pself->assert_is_active("set: ");
		#endif
		const PosIdx idx = this->index(pos_);
		set(idx, val_);
	}

	/**
	 * Set the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @param val_ value to copy into grid at pos_.
	 */
	void set (const PosIdx pos_idx_, Leaf val_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_idx_, "set: ");
		pself->assert_is_active("set: ");
		#endif
		pself->data()[pos_idx_] = val_;
	}
};


template <class TDerived>
class ByRef : protected Index<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits<TDerived>::Leaf;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;

protected:

	/**
	 * Get a reference to the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	Leaf& get (const VecDi& pos_)
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const PosIdx idx = this->index(pos_);
		return get(idx);
	}

	/**
	 * Get a const reference to the value stored in the grid.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	const Leaf& get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		const PosIdx idx = this->index(pos_);
		return get(idx);
	}

	/**
	 * Get a reference to the value stored in the grid.
	 *
	 * @param pos_idx_ data index of position to query.
	 * @return internally stored value at given grid position
	 */
	Leaf& get (const PosIdx pos_idx_)
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
	const Leaf& get (const PosIdx pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(this->index(pos_idx_), "get: ");
		#endif
		return pself->data()[pos_idx_];
	}

};


template <class TDerived>
class LazyByValue : protected ByValue<TDerived>
{
private:
	/// CRTP derived class.
	using Derived = TDerived;
	/// Dimension of the grid.
	static const UINT t_dims = Traits<TDerived>::t_dims;
	/// Type of data to store in grid nodes.
	using Leaf = typename Traits<TDerived>::Leaf;
	/// D-dimensional signed integer vector.
	using VecDi = Felt::VecDi<t_dims>;
	/// Base class, for subclass friending.
	using Base = ByValue<TDerived>;
protected:

	/**
	 * Get a copy of the value stored in the grid or background value if inactive.
	 *
	 * @param pos_ position in grid to query.
	 * @return internally stored value at given grid position
	 */
	Leaf get (const VecDi& pos_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pos_, "get: ");
		#endif
		if (pself->data().size())
		{
			return Base::get(pos_);
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
	Leaf get (const PosIdx pos_idx_) const
	{
		#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
		pself->assert_pos_bounds(pself->index(pos_idx_), "get: ");
		#endif
		if (pself->data().size())
		{
			return Base::get(pos_idx_);
		}
		else
		{
			return pself->m_background;
		}
	}
};

} // Access


} // Grid
} // Mixin
} // Impl
} // Felt

#endif // FELT_IMPL_GRID_HPP
