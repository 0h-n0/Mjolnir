#define BOOST_TEST_MODULE "test_excluded_volume_potential"

#ifdef UNITTEST_FRAMEWORK_LIBRARY_EXIST
#include <boost/test/unit_test.hpp>
#else
#define BOOST_TEST_NO_LIB
#include <boost/test/included/unit_test.hpp>
#endif

#include <mjolnir/potential/ExcludedVolumePotential.hpp>
#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/util/make_unique.hpp>

typedef mjolnir::DefaultTraits traits;

constexpr static unsigned int      seed = 32479327;
constexpr static std::size_t       N    = 10000;
constexpr static traits::real_type h    = 1e-6;

BOOST_AUTO_TEST_CASE(EXV_constructable)
{
    std::unique_ptr<mjolnir::GlobalPotentialBase<traits>> exv = 
        mjolnir::make_unique<mjolnir::ExcludedVolumePotential<traits>>();

    BOOST_CHECK(exv);
}

BOOST_AUTO_TEST_CASE(EXV_derivative)
{
    auto exv = mjolnir::make_unique<mjolnir::ExcludedVolumePotential<traits>>();

    const traits::real_type sigma   = 3.0;
    const traits::real_type epsilon = 1.0;
    exv->epsilon() = epsilon;
    exv->emplace(sigma);
    exv->emplace(sigma);

    const traits::real_type x_min = 0.8 * sigma;
    const traits::real_type x_max =
        mjolnir::ExcludedVolumePotential<traits>::cutoff_ratio * sigma;
    const traits::real_type dx = (x_max - x_min) / N;

    traits::real_type x = x_min;
    for(std::size_t i=0; i<N; ++i)
    {
        const traits::real_type pot1 = exv->potential(0, 1, x + h);
        const traits::real_type pot2 = exv->potential(0, 1, x - h);
        const traits::real_type dpot = (pot1 - pot2) / (2 * h);

        const traits::real_type deri = exv->derivative(0, 1, x);

        BOOST_CHECK_CLOSE_FRACTION(dpot, deri, h);
        x += dx;
    }
}
