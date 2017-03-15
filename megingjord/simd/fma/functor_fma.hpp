#ifndef MEGINGJORD_SIMD_FUNCTOR_FMA
#define MEGINGJORD_SIMD_FUNCTOR_FMA

#ifndef MEGINGJORD_SIMD_PACK
#error "DO NOT USE __FILE__ alone"
#endif

#ifndef MJOLNIR_HAVE_AVX2
#error "no AVX2 on the architecture"
#endif

namespace megingjord
{
namespace simd
{

// fmadd
template<>
struct fmadd_impl<__m128>
{
    static inline
    __m128 invoke(__m128 a, __m128 b, __m128 c)
    {
        return _mm_fmadd_ps(a, b, c);
    }
};

template<>
struct fmadd_impl<__m128d>
{
    static inline
    __m128d invoke(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fmadd_pd(a, b, c);
    }
};

template<>
struct fmadd_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a, __m256 b, __m256 c)
    {
        return _mm256_fmadd_ps(a, b, c);
    }
};

template<>
struct fmadd_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fmadd_pd(a, b, c);
    }
};

// fnmadd

template<>
struct fnmadd_impl<__m128>
{
    static inline
    __m128 invoke(__m128 a, __m128 b, __m128 c)
    {
        return _mm_fnmadd_ps(a, b, c);
    }
};

template<>
struct fnmadd_impl<__m128d>
{
    static inline
    __m128d invoke(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fnmadd_pd(a, b, c);
    }
};

template<>
struct fnmadd_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a, __m256 b, __m256 c)
    {
        return _mm256_fnmadd_ps(a, b, c);
    }
};

template<>
struct fnmadd_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fnmadd_pd(a, b, c);
    }
};


// fmsub
template<>
struct fmsub_impl<__m128>
{
    static inline
    __m128 invoke(__m128 a, __m128 b, __m128 c)
    {
        return _mm_fmsub_ps(a, b, c);
    }
};

template<>
struct fmsub_impl<__m128d>
{
    static inline
    __m128d invoke(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fmsub_pd(a, b, c);
    }
};

template<>
struct fmsub_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a, __m256 b, __m256 c)
    {
        return _mm256_fmsub_ps(a, b, c);
    }
};

template<>
struct fmsub_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fmsub_pd(a, b, c);
    }
};

// fnmsub

template<>
struct fnmsub_impl<__m128>
{
    static inline
    __m128 invoke(__m128 a, __m128 b, __m128 c)
    {
        return _mm_fnmsub_ps(a, b, c);
    }
};

template<>
struct fnmsub_impl<__m128d>
{
    static inline
    __m128d invoke(__m128d a, __m128d b, __m128d c)
    {
        return _mm_fnmsub_pd(a, b, c);
    }
};

template<>
struct fnmsub_impl<__m256>
{
    static inline
    __m256 invoke(__m256 a, __m256 b, __m256 c)
    {
        return _mm256_fnmsub_ps(a, b, c);
    }
};

template<>
struct fnmsub_impl<__m256d>
{
    static inline
    __m256d invoke(__m256d a, __m256d b, __m256d c)
    {
        return _mm256_fnmsub_pd(a, b, c);
    }
};

} // simd
} // megingjord
#endif /* MEGINGJORD_SIMD_PACK_FMA */
