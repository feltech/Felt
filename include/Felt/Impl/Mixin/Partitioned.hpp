#ifndef INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_
#define INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_

#include <Felt/Impl/Tracked.hpp>

namespace Felt
{
namespace Impl
{
namespace Mixin
{
namespace Partitioned
{

template <class Derived>
class Children
{
private:
	/// Traits of derived class.
	using TraitsType = Traits<Derived>;
	/// Child grid type.
	using ChildType = typename TraitsType::ChildType;
	/// Number of tracking lists.
	static constexpr UINT NumLists = TraitsType::NumLists;
	/// Dimension of grid.
	static constexpr UINT Dims = TraitsType::Dims;
	/// D-dimensional integer vector.
	using VecDi = Felt::VecDi<Dims>;
protected:
	/// Grid of partitions with tracking list(s) of active partitions.
	using ChildrenGrid = Impl::Tracked::MultiByRef<ChildType, Dims, NumLists>;
	/// Grid of child grids.
	ChildrenGrid m_children;


protected:
	Children(
		const VecDi& size_, const VecDi& offset_, const VecDi& child_size_,
		const ChildType& background_
	) :
		m_children(
			calc_children_size(size_, child_size_),
			(offset_.array() / child_size_.array()).matrix(),
			background_
		)
	{}


	const ChildrenGrid& children() const
	{
		return m_children;
	}


	ChildrenGrid& children()
	{
		return m_children;
	}

private:
	VecDi calc_children_size(const VecDi& size_, const VecDi& child_size_)
	{
		VecDi children_size = (size_.array() / child_size_.array()).matrix();

		if ((children_size.array() * child_size_.array()).matrix() != size_)
		{
			children_size += VecDi::Constant(1);
		}

		return children_size;
	}
};


} // Partitioned.
} // Mixin.
} // Impl.
} // Felt.


#endif /* INCLUDE_FELT_IMPL_MIXIN_PARTITIONED_HPP_ */
