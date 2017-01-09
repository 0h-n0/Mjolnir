#ifndef MJOLNIR_BOND_LENGTH_INTERACTION
#define MJOLNIR_BOND_LENGTH_INTERACTION
#include "Particle.hpp"
#include "LocalPotentialBase.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

/*! @brief calculate energy and force of Bond length type local interaction */
template<typename traitsT>
class BondLengthInteraction
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef LocalPotentialBase<traits_type> potential_type;

  public:

    BondLengthInteraction() = default;
    ~BondLengthInteraction() = default;

    void
    calc_force(particle_type& p1, particle_type& p2, const potential_type& pot) const;

    real_type
    calc_energy(const particle_type& p1, const particle_type& p2,
                const potential_type& pot) const;
};

template<typename traitsT>
inline void
BondLengthInteraction<traitsT>::calc_force(particle_type& p1, particle_type& p2,
        const potential_type& pot) const
{
    const coordinate_type dpos = p2.position - p1.position;
    const real_type lensq = length_sq(dpos);
    const real_type f = -1 * pot.derivative(std::sqrt(lensq));
    const coordinate_type force = dpos * (fast_inv_sqrt(lensq) * f);

#ifdef MJOLNIR_PARALLEL_THREAD
    std::lock_guard<std::mutex> lock1(p1.mtx);
    std::lock_guard<std::mutex> lock2(p2.mtx);
#endif

    p1.force -= force;
    p2.force += force;
    return;
}

template<typename traitsT>
inline typename BondLengthInteraction<traitsT>::real_type
BondLengthInteraction<traitsT>::calc_energy(
        const particle_type& p1, const particle_type& p2,
        const potential_type& pot) const
{
    return pot.potential(length(p1.position - p2.position));
}

} // mjolnir
#endif /* MJOLNIR_BOND_LENGTH_INTERACTION */
