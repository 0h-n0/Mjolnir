#ifndef MJOLNIR_FORCE_FIELD
#define MJOLNIR_FORCE_FIELD
#include "LocalForceField.hpp"
#include "GlobalForceField.hpp"

namespace mjolnir
{

template<typename traitsT>
class ForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type       real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef System<traits_type>           system_type;
    typedef LocalForceField<traits_type>  local_forcefield_type;
    typedef GlobalForceField<traits_type> global_forcefield_type;

  public:

    ForceField(local_force_field_type&& local, global_force_field_type&& global)
        : local_(std::move(local)), global_(std::move(global))
    {}
    ~ForceField() = default;

    ForceField(const ForceField&) = delete;
    ForceField(ForceField&&)      = default;
    ForceField& operator=(const ForceField&) = delete;
    ForceField& operator=(ForceField&&)      = default;

    void initialize(const system_type& sys, const time_type dt)
    {
        global_.initialize(sys, dt);
    }

    void calc_force(system_type&) const
    {
        local_.calc_force(sys);
        global_.calc_force(sys);
    }
    energy_type calc_energy(const system_type&) const
    {
        return local_.calc_energy(sys) + global_.calc_energy(sys);
    }

  private:

    local_force_field_type  local_;
    global_force_field_type global_;
};

} // mjolnir
#endif /* MJOLNIR_FORCE_FIELD */
