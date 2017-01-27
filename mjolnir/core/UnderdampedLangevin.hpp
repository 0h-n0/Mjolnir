#ifndef MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#define MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR
#include "Integrator.hpp"
#include "RandomNumberGenerator.hpp"
#include <mjolnir/util/zip_iterator.hpp>
#include <mjolnir/util/make_zip.hpp>
#include <memory>

#ifdef MJOLNIR_PARALLEL_THREAD
#include <future>
#endif

namespace mjolnir
{

template<typename traitsT>
class UnderdampedLangevin : public Integrator<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;

  public:

    UnderdampedLangevin(const time_type dt, const std::size_t number_of_particles,
            const real_type temperature, const real_type kB,
            std::vector<real_type>&& friction_constant,
            const std::shared_ptr<RandomNumberGenerator<traits_type>>& rng)
        : dt_(dt), halfdt_(dt * 0.5), halfdt2_(dt * dt * 0.5), rng_(rng),
          temperature_(temperature), kB_(kB),
          gamma_(std::forward<std::vector<real_type>>(friction_constant)),
          noise_(number_of_particles), acceleration_(number_of_particles)
    {}
    ~UnderdampedLangevin() override = default;

    void initialize(const ParticleContainer<traitsT>& pcon) override;

    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
                   const ForceField<traitsT>& ff) override;

#ifdef MJOLNIR_PARALLEL_THREAD
    time_type step(const time_type time, ParticleContainer<traitsT>& pcon,
       const ForceField<traitsT>& ff, const std::size_t num_threads) override;
#endif // MJOLNIR_PARALLEL_THREAD

    real_type  T() const {return temperature_;}
    real_type& T()       {return temperature_;}

    time_type  delta_t() const override {return dt_;}
    time_type& delta_t()       override {return dt_;}

    void generate_noise(const ParticleContainer<traitsT>& pcon);

  private:
    time_type dt_;      //!< dt
    time_type halfdt_;  //!< dt/2
    time_type halfdt2_; //!< dt^2/2
    real_type kB_;
    real_type temperature_;
    std::shared_ptr<RandomNumberGenerator<traits_type>> rng_;
    std::vector<real_type> gamma_;
    std::vector<coordinate_type> noise_;
    std::vector<coordinate_type> acceleration_;
};

template<typename traitsT>
void UnderdampedLangevin<traitsT>::initialize(
        const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(),  gamma_.cbegin(),
                             noise_.begin(), acceleration_.begin());
            iter != make_zip(pcon.cend(),    gamma_.cend(),
                             noise_.end(),   acceleration_.end()); ++iter)
    {
        // set a = f/m
        *get<3>(iter) = get<0>(iter)->force / get<0>(iter)->mass;

        // set random force
        *get<2>(iter) = rng_->underdamped_langevin(
                get<0>(iter)->mass, *get<1>(iter), dt_, temperature_, kB_);
    }
    return;
}

// at the initial step, acceleration_ and noise_ must be initialized
template<typename traitsT>
typename UnderdampedLangevin<traitsT>::time_type
UnderdampedLangevin<traitsT>::step(const time_type time,
        ParticleContainer<traitsT>& pcon, const ForceField<traitsT>& ff)
{
    // calc r(t+dt)
    for(auto iter = make_zip(pcon.begin(),    gamma_.cbegin(),
                             noise_.cbegin(), acceleration_.cbegin());
            iter != make_zip(pcon.end(),      gamma_.cend(),
                             noise_.cend(),   acceleration_.cend()); ++iter)
    {
        const real_type hgdt   = (*get<1>(iter) * halfdt_);
        const real_type o_hgdt = 1. - hgdt;

        const coordinate_type noisy_force = ((*get<2>(iter)) + (*get<3>(iter)));

        get<0>(iter)->position +=
            (dt_ * o_hgdt) * (get<0>(iter)->velocity) + halfdt2_ * noisy_force;

        get<0>(iter)->velocity *= o_hgdt * (o_hgdt * o_hgdt + hgdt);
        get<0>(iter)->velocity += (halfdt_ * o_hgdt) * noisy_force;
    }

    // calc f(t+dt)
    ff.calc_force(pcon);

    // calc a(t+dt) and v(t+dt), generate noise
    for(auto iter = make_zip(pcon.begin(),   gamma_.cbegin(),
                             noise_.begin(), acceleration_.begin());
            iter != make_zip(pcon.end(),     gamma_.cend(),
                             noise_.end(),   acceleration_.end()); ++iter)
    {
        // consider cash 1/m
        const coordinate_type acc = get<0>(iter)->force / (get<0>(iter)->mass);
        *get<3>(iter) = acc;

        const coordinate_type noise = rng_->underdamped_langevin(
                get<0>(iter)->mass, *get<1>(iter), dt_, temperature_, kB_);
        *get<2>(iter) = noise;

        const real_type gm = (*get<1>(iter));
        get<0>(iter)->velocity +=
            halfdt_ * (1. - gm * halfdt_) * (acc + noise);
        get<0>(iter)->force = coordinate_type(0., 0., 0.);
    }

    return time + dt_;
}

template<typename traitsT>
void UnderdampedLangevin<traitsT>::generate_noise(
        const ParticleContainer<traitsT>& pcon)
{
    for(auto iter = make_zip(pcon.cbegin(), gamma_.cbegin(), noise_.begin());
            iter != make_zip(pcon.cend(),   gamma_.cend(),   noise_.end());
            ++iter)
    {
        *get<2>(iter) = rng_->underdamped_langevin(
                get<0>(iter)->mass, *get<1>(iter), dt_, temperature_, kB_);
    }
    return;
}

#ifdef MJOLNIR_PARALLEL_THREAD
template<typename traitsT>
typename UnderdampedLangevin<traitsT>::time_type
UnderdampedLangevin<traitsT>::step(
        const time_type time, ParticleContainer<traitsT>& pcon,
        const ForceField<traitsT>& ff, const std::size_t num_threads)
{
    std::vector<std::future<void>> futures(num_threads-1);

    // calc r(t+dt)
    const std::size_t num_particle = pcon.size();
    const std::size_t particle_per_thread = num_particle / num_threads + 1;
    std::size_t begin = 0;
    std::size_t end = particle_per_thread;
    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
    {
        *iter = std::async(std::launch::async, [this, &pcon, begin, end]()
            {
                for(std::size_t i=begin; i != end; ++i)
                {
                    const real_type hgdt   = gamma_[i] * halfdt_;
                    const real_type o_hgdt = 1. - hgdt;
                    const coordinate_type noisy_force = noise_[i] + acceleration_[i];

                    pcon[i].position +=
                        (dt_ * o_hgdt) * (pcon[i].velocity) + halfdt2_ * noisy_force;
                    pcon[i].velocity *= o_hgdt * (o_hgdt * o_hgdt + hgdt);
                    pcon[i].velocity += (halfdt_ * o_hgdt) * noisy_force;
                }
                return;
            });

        begin += particle_per_thread;
        end   += particle_per_thread;
    }

    // do "this" thread's job
    for(std::size_t i=begin; i != num_particle; ++i)
    {
        const real_type hgdt   = gamma_[i] * halfdt_;
        const real_type o_hgdt = 1. - hgdt;
        const coordinate_type noisy_force = noise_[i] + acceleration_[i];

        pcon[i].position +=
            (dt_ * o_hgdt) * (pcon[i].velocity) + halfdt2_ * noisy_force;
        pcon[i].velocity *= o_hgdt * (o_hgdt * o_hgdt + hgdt);
        pcon[i].velocity += (halfdt_ * o_hgdt) * noisy_force;
    }

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        iter->wait(); // synchronize

    // calc f(t+dt)
    ff.calc_force(pcon, num_threads);

    // generate random numbers
    this->generate_noise(pcon);

    // calc a(t+dt) and v(t+dt)
    begin = 0;
    end = particle_per_thread;
    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
    {
        *iter = std::async(std::launch::async, [this, &pcon, begin, end]()
            {
                for(std::size_t i=begin; i != end; ++i)
                {
                    const coordinate_type acc = pcon[i].force / pcon[i].mass;
                    acceleration_[i] = acc;
 
                    pcon[i].velocity += halfdt_ * (1.-gamma_[i]*halfdt_) * (acc+noise_[i]);
                    pcon[i].force = coordinate_type(0., 0., 0.);
                }
                return;
            });

        begin += particle_per_thread;
        end   += particle_per_thread;
    }

    for(std::size_t i=begin; i != num_particle; ++i)
    {
        const coordinate_type acc = pcon[i].force / pcon[i].mass;
        acceleration_[i] = acc;
 
        pcon[i].velocity += halfdt_ * (1.-gamma_[i]*halfdt_) * (acc+noise_[i]);
        pcon[i].force = coordinate_type(0., 0., 0.);
    }

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        iter->wait(); // synchronize

    return time + dt_;
}

#endif // MJOLNIR_PARALLEL_THREAD

} // mjolnir
#endif /* MJOLNIR_UNDERDAMPED_LANGEVIN_INTEGRATOR */
