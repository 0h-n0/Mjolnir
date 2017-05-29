#ifndef MJOLNIR_BOND_ANGLE_INTERACTION
#define MJOLNIR_BOND_ANGLE_INTERACTION
#include "LocalInteractionBase.hpp"
#include "BoundaryCondition.hpp"
#include "constants.hpp"
#include <mjolnir/math/fast_inv_sqrt.hpp>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename potentialT>
class BondAngleInteraction : public LocalInteractionBase<traitsT>
{
  public:
    typedef traitsT    traits_type;
    typedef potentialT potential_type;
    typedef LocalInteractionBase<traits_type>   base_type;
    typedef typename base_type::real_type       real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::system_type     system_type;
    typedef typename base_type::particle_type   particle_type;
    typedef std::array<std::size_t, 3>          indices_type;
    typedef std::pair<indices_type, potentialT> potential_index_pair;
    typedef std::vector<potential_index_pair>   container_type;

  public:

    BondAngleInteraction() = default;
    ~BondAngleInteraction() = default;

    BondAngleInteraction(const container_type& pot): potentials(pot){}
    BondAngleInteraction(container_type&& pot): potentials(std::move(pot)){}

    void      calc_force (system_type&)        const override;
    real_type calc_energy(const system_type& ) const override;

  private:
    container_type potentials;
};

template<typename traitsT, typename potentialT>
void
BondAngleInteraction<traitsT, potentialT>::calc_force(system_type& sys) const
{
    for(const auto& idxp : this->potentials)
    {
        const coordinate_type r_ij = sys.adjust_direction(
                sys[idxp.first[0]].position - sys[idxp.first[1]].position);
        const real_type       inv_len_r_ij = fast_inv_sqrt(length_sq(r_ij));
        const coordinate_type r_ij_reg     = r_ij * inv_len_r_ij;

        const coordinate_type r_kj = sys.adjust_direction(
                sys[idxp.first[2]].position - sys[idxp.first[1]].position);
        const real_type       inv_len_r_kj = fast_inv_sqrt(length_sq(r_kj));
        const coordinate_type r_kj_reg     = r_kj * inv_len_r_kj;

        const real_type dot_ijk = dot_product(r_ij_reg, r_kj_reg);
        const real_type cos_theta = (-1. <= dot_ijk && dot_ijk <= 1.)
            ? dot_ijk : std::copysign(1.0, dot_ijk);

        const real_type theta = std::acos(cos_theta);

        const real_type coef = -(idxp.second.derivative(theta));

        const real_type sin_theta = std::sin(theta);
        const real_type coef_inv_sin =
            (sin_theta > constants<real_type>::tolerance) ?
            coef / sin_theta : coef / constants<real_type>::tolerance;

        const coordinate_type Fi =
            (coef_inv_sin * inv_len_r_ij) * (cos_theta * r_ij_reg - r_kj_reg);

        const coordinate_type Fk =
            (coef_inv_sin * inv_len_r_kj) * (cos_theta * r_kj_reg - r_ij_reg);

        sys[idxp.first[0]].force += Fi;
        sys[idxp.first[1]].force -= (Fi + Fk);
        sys[idxp.first[2]].force += Fk;
    }
    return;
}

template<typename traitsT, typename potentialT>
typename BondAngleInteraction<traitsT, potentialT>::real_type
BondAngleInteraction<traitsT, potentialT>::calc_energy(
        const system_type& sys) const
{
    real_type E = 0.;
    for(const auto& idxp : this->potentials)
    {
        const coordinate_type v_2to1 = sys.adjust_direction(
                sys[idxp.first[0]].position - sys[idxp.first[1]].position);
        const coordinate_type v_2to3 = sys.adjust_direction(
                sys[idxp.first[2]].position - sys[idxp.first[1]].position);

        const real_type lensq_v21 = length_sq(v_2to1);
        const real_type lensq_v23 = length_sq(v_2to3);
        const real_type dot_v21_v23 = dot_product(v_2to1, v_2to3);

        const real_type dot_ijk =
            dot_v21_v23 * fast_inv_sqrt(lensq_v21 * lensq_v23);
        const real_type cos_theta = (-1. <= dot_ijk && dot_ijk <= 1.)
                                    ? dot_ijk : std::copysign(1.0, dot_ijk);
        const real_type theta = std::acos(cos_theta);

        E += idxp.second.potential(theta);
    }
    return E;
}

}// mjolnir
#endif /* MJOLNIR_BOND_ANGL_INTERACTION */
