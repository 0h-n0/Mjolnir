#ifndef MJOLNIR_VELOCITY_VERLET_INTEGRATOR
#define MJOLNIR_VELOCITY_VERLET_INTEGRATOR
#include "Integrator.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>

namespace mjolnir
{

template<typename traitsT>
class VelocityVerlet : public Integrator<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    VelocityVerlet(const time_type dt, const std::size_t number_of_particles)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5),
          acceleration_(number_of_particles)
    {}
    ~VelocityVerlet() override = default;

    void initialize(const ParticleContainer<traitsT>& pcon) override;

    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   const ForceField<traitsT>& ff) override;

    time_type& delta_t()       override {return dt_;}
    time_type  delta_t() const override {return dt_;}

#ifdef MJOLNIR_PARALLEL_THREAD
    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
       const ForceField<traitsT>& ff, const std::size_t num_threads) override;

  private:
    void step_pos_and_half_vel(ParticleContainer<traitsT>& pcon,
            std::size_t index, const std::size_t end);
    void step_acc_and_half_vel(ParticleContainer<traitsT>& pcon,
            std::size_t index, const std::size_t end);

#endif // MJOLNIR_PARALLEL_THREAD

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT>
void VelocityVerlet<traitsT>::initialize(const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(), acceleration_.begin());
            iter != make_zip(pcon.cend(), acceleration_.end()); ++iter)
    {
        *get<1>(iter) = get<0>(iter)->force / get<0>(iter)->mass;
    }
    return;
}

// at the initial step, acceleration_ must be initialized
template<typename traitsT>
typename VelocityVerlet<traitsT>::time_type
VelocityVerlet<traitsT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, const ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.cbegin());
            iter != make_zip(pcon.end(), acceleration_.cend()); ++iter)
    {
        get<0>(iter)->position += dt_ * (get<0>(iter)->velocity) +
                                  halfdt2_ * (*get<1>(iter));
        get<0>(iter)->velocity += halfdt_ * (*get<1>(iter));
    }

    // calc f(t+dt)
    ff.calc_force(pcon);

    // calc a(t+dt) and v(t+dt)
    for(auto iter = make_zip(pcon.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(), acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<1>(iter) = acc;
        get<0>(iter)->velocity += halfdt_ * acc;
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    return time + dt_;
}

#ifdef MJOLNIR_PARALLEL_THREAD
template<typename traitsT>
typename VelocityVerlet<traitsT>::time_type
VelocityVerlet<traitsT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, const ForceField<traitsT>& ff,
        const std::size_t num_threads)
{
    std::vector<std::future<void>> futures(num_threads-1);

    // calc r(t+dt)
    const std::size_t num_particle = pcon.size();
    const std::size_t particle_per_thread = num_particle / num_threads + 1;
    std::size_t begin = 0;
    std::size_t end = particle_per_thread;
    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        futures[t] = std::async(std::launch::async,
            [this, &pcon, begin, end](){step_pos_and_half_vel(pcon, begin, end);});

        begin += particle_per_thread;
        end   += particle_per_thread;
    }
    step_pos_and_half_vel(pcon, begin, num_particle);
    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        iter->get();

    // calc f(t+dt)
    ff.calc_force(pcon, num_threads);

    // calc a(t+dt) and v(t+dt)
    begin = 0;
    end = particle_per_thread;
    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        futures[t] = std::async(std::launch::async,
            [this, &pcon, begin, end](){step_acc_and_half_vel(pcon, begin, end);});

        begin += particle_per_thread;
        end   += particle_per_thread;
    }
    step_pos_and_half_vel(pcon, begin, num_particle);
    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        iter->get();

    return time + dt_;
}

template<typename traitsT>
void VelocityVerlet<traitsT>::step_pos_and_half_vel(
        ParticleContainer<traitsT>& pcon, std::size_t i, const std::size_t end)
{
    for(; i != end; ++i)
    {
        pcon[i].position += dt_ * pcon[i].velocity + halfdt2_ * acceleration_[i];
        pcon[i].velocity += halfdt_ * acceleration_[i];
    }
    return;
}

template<typename traitsT>
void VelocityVerlet<traitsT>::step_acc_and_half_vel(
        ParticleContainer<traitsT>& pcon, std::size_t i, const std::size_t end)
{
    for(; i != end; ++i)
    {
        const coordinate_type acc = pcon[i].force / pcon[i].mass;
        acceleration_[i] = acc;
        pcon[i].velocity += halfdt_ * acc;
        pcon[i].force = coordinate_type(0., 0., 0.);
    }
    return;
}

#endif

} // mjolnir
#endif /* MJOLNIR_VELOCITY_VERLET_INTEGRATOR */
