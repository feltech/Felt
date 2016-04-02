#ifndef INCLUDE_FELT_ALIGNEDARRAY_HPP_
#define INCLUDE_FELT_ALIGNEDARRAY_HPP_

#include <mutex>
#include <vector>
#include <eigen3/Eigen/Core>

namespace felt
{

/// An SSE compatible aligned array of arbitrary type with a mutex member for external thread-safety.
template <typename T>
class AlignedArray :  public std::vector<T, Eigen::aligned_allocator<T> >
{
protected:
	/// Mutex for external locking.
	std::mutex m_mutex;
public:
	/// Base class.
	using Base = std::vector<T, Eigen::aligned_allocator<T> >;

	/// Use constructor of std::vector.
	using Base::vector;

	AlignedArray() = default;

	AlignedArray(const AlignedArray&& other) : Base(other)
	{}

	AlignedArray(const AlignedArray& other) : Base(other)
	{}

	void operator=(const AlignedArray& other)
	{
		Base::operator=(other);
	}

	/// Get mutex member.
	std::mutex& mutex()
	{
		return m_mutex;
	}
};

}
#endif /* INCLUDE_FELT_ALIGNEDARRAY_HPP_ */
