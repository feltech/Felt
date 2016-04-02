#ifndef INCLUDE_FELT_PARTITIONEDARRAY_HPP_
#define INCLUDE_FELT_PARTITIONEDARRAY_HPP_

#include "PartitionBase.hpp"
#include "AlignedArray.hpp"

namespace felt
{

/**
 * Base class for common features of PartitionedArray template specialisations.
 */
template <class Derived>
class PartitionedArrayBase : public PartitionBase<PartitionedArrayBase<Derived> >
{
public:
	/// This class.
	using ThisType = PartitionedArrayBase<Derived>;
	/// Base class.
	using Base = PartitionBase<ThisType>;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using Base::Traits;
	using LeafType = typename Base::Traits::LeafType;
protected:
	/// Spatial offset of 'imaginary' grid containing the list.
	VecDi	m_offset;
public:
	/**
	 * Set offset of 'imaginary' grid containing the list.
	 *
	 * @param offset_
	 */
	void offset (const VecDi& offset_)
	{
		m_offset = offset_;
		Base::offset(offset_);
	}

	/**
	 * Get spatial partition from leaf grid node in 'imaginary' grid.
	 *
	 * @param pos_leaf_
	 * @return location of spatial partition in children grid.
	 */
	const VecDi pos_child (const VecDi& pos_leaf_) const
	{
		return (
			(pos_leaf_ - m_offset).array() / this->m_isize_child.array()
		).matrix() + this->children().offset();
	}
};


/**
 * Spatially partitioned expandable lists.
 *
 * A specialised partitioned grid, where the child grids are simply expandable lists.
 *
 * @tparam T the type to store in elements of the list
 * @tparam D the dimension of the 'imaginary' grid.
 * @tparam N the dimension of the array(s).
 */
template <typename T, UINT D, UINT N=0>
class PartitionedArray
	: public PartitionedArrayBase<PartitionedArray<T, D, N> >
{
protected:
	/// This class.
	using ThisType = PartitionedArray<T, D, N>;
	/// Base class
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;

public:
	/// Explicitly defined default constructor.
	PartitionedArray () = default;

	/**
	 * Construct multiple spatially partitioned arrays contained in an 'imaginary' grid.
	 *
	 * @param size_ spatial size of 'imaginary' grid.
	 * @param offset_ spatial offset of 'imaginary' grid.
	 * @param size_partition_ spatial size of a single partition.
	 */
	PartitionedArray (const VecDu& size_, const VecDi& offset_, const VecDu& size_partition_) :
		Base()
	{
		this->init(size_, offset_, size_partition_);
	}

	/**
	 * Add val to list, placing in partition found from pos.
	 *
	 * @param pos_ position in 'imaginary' grid.
	 * @param val_ value to insert in list.
	 * @param arr_idx_ ID of list to insert into.
	 */
	void add(const VecDi& pos_, const T& val_, const UINT arr_idx_)
	{
		const VecDi& pos_child = this->pos_child(pos_);
		this->children().get(pos_child)[arr_idx_].push_back(val_);
		Base::add_child(pos_child, arr_idx_);
	}

	/**
	 * Loop all spatial partitions, resizing the given list to zero in each.
	 *
	 * @param arr_idx_ ID of list to reset.
	 */
	void reset(const UINT arr_idx_)
	{
		for (const VecDi& pos_child : this->children().list(arr_idx_))
			this->children().get(pos_child)[arr_idx_].clear();
		Base::reset(arr_idx_);
	}
};


/**
 * Spatially partitioned expandable list - 1D array specialisation.
 *
 * A specialised partitioned grid, where the child grids are simply
 * expandable lists.
 *
 * @tparam T type to store in elements of the list
 * @tparam D dimension of the 'imaginary' grid.
 */
template <typename T, UINT D>
class PartitionedArray<T, D, 0>
	: public PartitionedArrayBase<PartitionedArray<T, D, 0> >
{
protected:
	/// This class
	using ThisType = PartitionedArray<T, D, 0>;
	/// Base class.
	using Base = PartitionedArrayBase<ThisType>;

	using typename Base::Child;
	using typename Base::VecDu;
	using typename Base::VecDi;
	using typename Base::LeafType;

public:
	/// Explicitly defined default constructor.
	PartitionedArray () = default;

	/**
	 * Construct a single (1D) spatially partitioned array contained in an 'imaginary' grid.
	 *
	 * @param size_ spatial size of 'imaginary' grid.
	 * @param offset_ spatial offset of 'imaginary' grid.
	 * @param size_partition_ spatial size of a single partition.
	 */
	PartitionedArray (const VecDu& size_, const VecDi& offset_, const VecDu& size_partition_) :
		Base()
	{
		this->init(size_, offset_, size_partition_);
	}

	/**
	 * Add val to list, placing in partition found from pos.
	 *
	 * @param pos position in 'imaginary' grid.
	 * @param val value to insert in list.
	 */
	void add(const VecDi& pos, const T& val)
	{
		const VecDi& pos_child = this->pos_child(pos);
		this->children().get(pos_child).push_back(val);
		Base::add_child(pos_child);
	}

	/**
	 * Thread safely add val to list, placing in partition found from pos.
	 *
	 * @param pos position in 'imaginary' grid.
	 * @param val value to insert in list.
	 */
	void add_safe(const VecDi& pos, const T& val)
	{
		const VecDi& pos_child = this->pos_child(pos);
		Child& child = this->children().get(pos_child);
		Base::add_child(pos_child);
		std::lock_guard<std::mutex> lock(child.mutex());
		child.push_back(val);
	}

	/**
	 * Loop all spatial partitions, resizing their lists to zero.
	 */
	void reset()
	{
		for (const VecDi& pos_child : this->children().list())
			this->children().get(pos_child).clear();
		Base::reset();
	}
};


/**
 * Traits for PartitionedArrayBase.
 *
 * Just forward the traits defined for PartitionedArrayBase subclasses.
 */
template <class Derived>
struct GridTraits<PartitionedArrayBase<Derived> > : GridTraits<Derived>
{};


/**
 * Traits for PartitionedArray.
 *
 * @tparam T type to store in elements of the list
 * @tparam D dimension of the 'imaginary' grid.
 * @tparam N number of distinct arrays to partition.
 */
template <typename T, UINT D, UINT N>
struct GridTraits<PartitionedArray<T, D, N> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, N>;
	/// Child type to store in spatial partitions - in this case an array of lists.
	using ChildType = std::array<AlignedArray<T>, N>;
	/// Number of distinct arrays to partition - in this case a single list.
	static const UINT NumLists = N;
};


/**
 * Traits for PartitionedArray.
 *
 * @tparam T data type stored in array.
 * @tparam D dimension of the grid.
 */
template <typename T, UINT D>
struct GridTraits<PartitionedArray<T, D, 0> > : DefaultGridTraits<T, D>
{
	/// The class inheriting from the base.
	using ThisType = PartitionedArray<T, D, 0>;
	/// Child type to store in spatial partitions - in this case a single list.
	using ChildType = AlignedArray<T>;
	/// Number of distinct arrays to partition - in this case a single list.
	static const UINT NumLists = 1;
};

} // End namespace felt.

#endif /* INCLUDE_FELT_PARTITIONEDARRAY_HPP_ */
