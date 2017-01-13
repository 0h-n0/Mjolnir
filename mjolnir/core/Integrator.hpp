#ifndef MJOLNIR_INTEGRATOR
#define MJOLNIR_INTEGRATOR
#include "ParticleContainer.hpp"
#include "ForceField.hpp"

namespace mjolnir
{
template<typename traitsT>
class Integrator
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;

  public:

    virtual ~Integrator() = default;

    virtual void initialize(const ParticleContainer<traitsT>& pcon) = 0;

    virtual time_type
    step(const time_type time, ParticleContainer<traitsT>& pcon,
         const ForceField<traitsT>& ff) = 0;

#ifdef MJOLNIR_PARALLEL_THREAD
    virtual time_type
    step(const time_type time, ParticleContainer<traitsT>& pcon,
         const ForceField<traitsT>& ff, const std::size_t num_threads) = 0;
#endif //MJOLNIR_PARALLEL_THREAD

    virtual time_type& delta_t()       = 0;
    virtual time_type  delta_t() const = 0;
};

} // mjolnir
#endif /* MJOLNIR_INTEGRATOR */
