#ifndef MJOLNIR_GO_10_12_CONTACT_POTENTIAL
#define MJOLNIR_GO_10_12_CONTACT_POTENTIAL
#include <mjolnir/core/LocalPotentialBase.hpp>

namespace mjolnir
{

/*! @brief Go contact 10-12 potential *
 *  V(r) = epsilon * (5 * (r0/r)^12 - 6 * (r0/r)^10)  *
 * dV/dr = 60 * epsilon * ((r0/r)^10 - (r0/r)^12) / r */
template<typename traitsT>
class Go1012ContactPotential : public LocalPotentialBase<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:
    Go1012ContactPotential(const real_type e, const real_type r0)
        : epsilon_(e), r0_(r0)
    {}
    ~Go1012ContactPotential() override = default;

    real_type potential(const real_type r) const override
    {
        const real_type rd   = this->r0_ / r;
        const real_type rd2  = rd  * rd;
        const real_type rd4  = rd2 * rd2;
        const real_type rd8  = rd4 * rd4;
        const real_type rd10 = rd8 * rd2;
        const real_type rd12 = rd8 * rd4;

        return this->epsilon_ * (5. * rd12 - 6. * rd10);
    }

    real_type derivative(const real_type r) const override
    {
        const real_type invr = 1. / r;
        const real_type rd   = this->r0_ * invr;
        const real_type rd2  = rd  * rd;
        const real_type rd4  = rd2 * rd2;
        const real_type rd8  = rd4 * rd4;
        const real_type rd10 = rd8 * rd2;
        const real_type rd12 = rd8 * rd4;
        return this->epsilon_ * 60. * (rd10 - rd12) * invr;
    }

  private:

    const real_type epsilon_;
    const real_type r0_;
};

}
#endif /* MJOLNIR_GO_10_12_CONTACT_POTENTIAL */
