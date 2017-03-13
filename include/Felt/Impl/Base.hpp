#ifndef INCLUDE_FELT_IMPL_BASE_HPP_
#define INCLUDE_FELT_IMPL_BASE_HPP_
namespace Felt
{
namespace Impl
{

template <class Derived> struct Traits {};

template <class Derived>
class Base
{
protected:
	Derived* self()
	{
		return (Derived*)(this);
	}
	const Derived* self() const
	{
		return (const Derived*)(this);
	}
};

}
}
#endif /* INCLUDE_FELT_IMPL_BASE_HPP_ */
