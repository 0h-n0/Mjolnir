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
// __m256, __m256d
#include <immintrin.h>
#include "avx/pack_avx.hpp"
#include "avx/functor_avx.hpp"
#endif// avx

#ifdef MJOLNIR_HAVE_AVX2
// __m256i and FMA operation
// TODO #include "avx2/pack_avx2.hpp"
//      #include "avx2/functor_avx2.hpp"
#endif// avx2

#ifdef MJOLNIR_HAVE_FMA
#include "fma/functor_fma.hpp"
#elif defined(MJOLNIR_HAVE_AVX)
// TODO #include "avx/fallback_fma.hpp"
#elif defined(MJOLNIR_HAVE_SSE4_2)
// TODO #include "sse/fallback_fma.hpp"
#endif

// XXX set default simd operation
// TODO: AVX512
#if defined(MJOLNIR_HAVE_AVX2) || defined(MJOLNIR_HAVE_AVX)
#    define MEGINGJORD_DEFAULT_SIMD avx_traits
// TODO: #elif defined(MJOLNIR_HAVE_SSE4_2)
//       #    define MEGINGJORD_DEFAULT_SIMD sse_traits
#else // fallback for no known SIMD
#    include "plain/pack_plain.hpp"
#    define MEGINGJORD_DEFAULT_SIMD plain_traits
#endif

#include "packable_array.hpp"
#include "operation_pack.hpp"
#include "operation_array.hpp"

#endif /* MEGINGJORD_SIMD_PACK */
