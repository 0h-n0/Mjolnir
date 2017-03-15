#ifndef MEGINGJORD_SIMD_PACK
#define MEGINGJORD_SIMD_PACK
#include <megingjord/util/aligned_array.hpp>

namespace megingjord
{
namespace simd
{

template<typename, std::size_t> struct pack;
template<typename, std::size_t, typename> struct packed_array;

template<typename> struct is_simd : public std::false_type{};
template<typename> struct is_packable : public std::false_type{};
template<> struct is_packable<float>  : public std::true_type{};
template<> struct is_packable<double> : public std::true_type{};

template<typename> struct single_type_of;
template<typename T, std::size_t N>
struct single_type_of<pack<T, N>>
{
    typedef T type;
};
template<typename T, std::size_t N, std::size_t align>
struct single_type_of<aligned_array<T, N, align>>
{
    typedef T type;
};
template<typename T, std::size_t N, typename trait>
struct single_type_of<packed_array<T, N, trait>>
{
    typedef T type;
};

} // simd
} // megingjord

#include "functors.hpp"

#ifdef MJOLNIR_HAVE_AVX
#include <immintrin.h>
#include "avx/pack_avx.hpp"
#include "avx/functor_avx.hpp"
#define MEGINGJORD_DEFAULT_SIMD avx_traits
#endif

#ifndef MJOLNIR_HAVE_AVX
#include "plain/pack_plain.hpp"
#define MEGINGJORD_DEFAULT_SIMD plain_traits
#endif

#include "packable_array.hpp"
#include "operation_pack.hpp"
#include "operation_array.hpp"

#endif /* MEGINGJORD_SIMD_PACK */
