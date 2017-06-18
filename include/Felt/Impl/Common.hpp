#ifndef INCLUDE_FELT_IMPL_COMMON_HPP_
#define INCLUDE_FELT_IMPL_COMMON_HPP_

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>

namespace Felt
{
namespace Impl
{
template <class Derived> struct Traits {};

/**
 * Helper to unwrap parenthesised argument.
 */
#define _FELT_UNWRAP(...) __VA_ARGS__

/**
 * Helper to add private inheritance of mixin to derived, comma-separated.
 */
#define _FELT_MIXIN_INHERIT(r, derived, idx, mixin)\
	BOOST_PP_COMMA_IF(idx) private Impl::Mixin::mixin<_FELT_UNWRAP derived>

/**
 * Helper to make derived a `friend` of mixin.
 */
#define _FELT_MIXIN_GRANT(r, derived, mixin)\
	friend Impl::Mixin::mixin<_FELT_UNWRAP derived>;

/**
 * Add an upcast method for use by mixins.
 */
#define _FELT_ENABLE_UPCAST(derived)\
	public:\
	template <class Base>\
	static _FELT_UNWRAP derived * upcast(Base* base_)\
	{\
		return static_cast<_FELT_UNWRAP derived *>(base_);\
	}\
	template <class Base>\
	static const _FELT_UNWRAP derived * upcast(const Base* base_)\
	{\
		return static_cast<_FELT_UNWRAP derived *>(base_);\
	}\
	private:

/**
 * Add inheritance of a sequence of derived classes.
 *
 * Sequence is defined as in BOOST_PP_SEQ style: (arg1)(arg2)(arg3)...
 *
 * Variadic part allows additional `friend`s - useful for derived classes of mixins.
 */
#define FELT_MIXINS(derived, mixins, ...)\
	BOOST_PP_SEQ_FOR_EACH_I(_FELT_MIXIN_INHERIT, derived, mixins)\
	{\
	BOOST_PP_SEQ_FOR_EACH(_FELT_MIXIN_GRANT, derived, mixins)\
	BOOST_PP_SEQ_FOR_EACH(_FELT_MIXIN_GRANT, derived, __VA_ARGS__)\
	_FELT_ENABLE_UPCAST(derived)

/**
 * Add `friend` access to a sequence of derived classes.
 *
 * Sequence is defined as in BOOST_PP_SEQ style: (arg1)(arg2)(arg3)...
 */
#define FELT_MIXINS_GRANT(derived, mixins) BOOST_PP_SEQ_FOR_EACH(_FELT_MIXIN_GRANT, derived, mixins)

} // Impl
} // Felt
#endif /* INCLUDE_FELT_IMPL_COMMON_HPP_ */
