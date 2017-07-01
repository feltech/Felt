#ifndef INCLUDE_FELT_IMPL_COMMON_HPP_
#define INCLUDE_FELT_IMPL_COMMON_HPP_

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/punctuation/comma_if.hpp>
#include <eigen3/Eigen/Dense>


namespace Felt
{
/**
 * Use 32 bit float by default.
 */
using FLOAT = float;

/**
 * Use 32 bit int by default.
 */
using INT = int;

/**
 * Use 32 bit unsigned int by default.
 */
using UINT = unsigned;

/**
 * Integer (large) array index type.
 */
using Idx = UINT;

/**
 * Shorthand for D-dimensional vector with elements of T type.
 */
template <typename T, UINT D>
using VecDT = Eigen::Matrix<T, D, 1>;
/**
 * Shorthand for D-dimensional float vector.
 */
template <UINT D>
using VecDf = VecDT<FLOAT, D>;
/**
 * Shorthand for D-dimensional integer vector.
 */
template <UINT D>
using VecDi = VecDT<INT, D>;
/**
 * Shorthand for D-dimensional unsigned integer vector.
 */
template <UINT D>
using VecDu = VecDT<UINT, D>;

/**
 * Shorthand for 2D float vector.
 */
using Vec2f = VecDf<2u>;
/**
 * Shorthand for 2D unsigned integer vector.
 */
using Vec2u = VecDu<2u>;
/**
 * Shorthand for 2D integer vector.
 */
using Vec2i = VecDi<2u>;
/**
 * Shorthand for 3D float vector.
 */
using Vec3f = VecDf<3u>;
/**
 * Shorthand for 3D unsigned integer vector.
 */
using Vec3u = VecDu<3u>;
/**
 * Shorthand for 3D integer vector.
 */
using Vec3i = VecDi<3u>;

namespace Impl
{
template <class Derived> struct Traits {};
} // Impl
} // Felt


/**
 * Shorthand cast to CRTP derived type.
 */
#define pself Derived::upcast(this)

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
		return static_cast<const _FELT_UNWRAP derived *>(base_);\
	}\
	private:

/**
 * Add inheritance of a sequence of derived classes.
 *
 * Sequence is defined as in BOOST_PP_SEQ style: (arg1)(arg2)(arg3)...
 *
 * Variadic part allows additional `friend`s - useful for base classes of `mixins` that require
 * access to `derived`.
 */
#define FELT_MIXINS(derived, mixins, ...)\
	BOOST_PP_SEQ_FOR_EACH_I(_FELT_MIXIN_INHERIT, derived, mixins)\
	{\
	BOOST_PP_SEQ_FOR_EACH(_FELT_MIXIN_GRANT, derived, mixins)\
	BOOST_PP_SEQ_FOR_EACH(_FELT_MIXIN_GRANT, derived, __VA_ARGS__)\
	_FELT_ENABLE_UPCAST(derived)

#if defined(FELT_EXCEPTIONS) || !defined(NDEBUG)
#define FELT_DEBUG(...) __VA_ARGS__;
#else
#define FELT_DEBUG(...)
#endif


#endif /* INCLUDE_FELT_IMPL_COMMON_HPP_ */
