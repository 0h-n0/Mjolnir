#ifndef MEGINGJORD_SIMD_FUNCTORS
#define MEGINGJORD_SIMD_FUNCTORS
#include <cmath>

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE this header alone"
#endif

namespace megingjord
{
namespace simd
{

template<typename T> struct set_impl
{
    constexpr static inline
    T invoke(T x){return x;}
};

template<typename T> struct load_impl
{
    constexpr static inline
    T invoke(const T *x){return *x;}
};

template<typename T> struct broadcast_impl
{
    constexpr static inline
    T invoke(const T *x){return *x;}
};

template<typename T> struct store_impl
{
    constexpr static inline
    T invoke(T x){return x;}

    static inline
    void invoke(T* dst, T x){*dst = x;}
};

template<typename T> struct add_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs + rhs;}
};

template<typename T> struct sub_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs - rhs;}
};

template<typename T> struct mul_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs * rhs;}
};

template<typename T> struct div_impl
{
    constexpr static inline
    T invoke(T lhs, T rhs) {return lhs / rhs;}
};

template<typename T> struct rcp_impl
{
    constexpr static inline
    T invoke(T a) {return 1 / a;}
};

template<typename T> struct rsqrt_impl
{
    static inline
    T invoke(T a) {return 1 / std::sqrt(a);}
};

template<typename T> struct sqrt_impl
{
    static inline
    T invoke(T a) {return std::sqrt(a);}
};

template<typename T> struct floor_impl
{
    static inline
    T invoke(T a) {return std::floor(a);}
};

template<typename T> struct ceil_impl
{
    static inline
    T invoke(T a) {return std::ceil(a);}
};

template<typename T> struct fmadd_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return a * b + c;}
};

template<typename T> struct fnmadd_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return -a * b + c;}
};

template<typename T> struct fmsub_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return a * b - c;}
};

template<typename T> struct fnmsub_impl
{
    constexpr static inline
    T invoke(T a, T b, T c) {return -a * b - c;}
};

} // simd
} // megingjord
#endif // MEGINGJORD_SIMD_FUNCTORS
