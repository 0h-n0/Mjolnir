#ifndef MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#define MJOLNIR_GLOBAL_DISTANCE_INTEARACTION
#include "GlobalInteractionBase.hpp" 
#include "SpatialPartition.hpp" 
#include <memory>

namespace mjolnir
{

template<typename traitsT>
class GlobalDistanceInteraction : public GlobalInteractionBase<traitsT>
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef GlobalPotentialBase<traits_type> potential_type;
    typedef SpatialPartition<traits_type> spatial_partitioning_type;

  public:
    GlobalDistanceInteraction() = default;
    GlobalDistanceInteraction(std::unique_ptr<spatial_partitioning_type>&& sp)
        : spatial_partition_(
                std::forward<std::unique_ptr<spatial_partitioning_type>>(sp))
    {}
    ~GlobalDistanceInteraction() = default;

    void initialize(const particle_container_type& pcon,
                    const time_type dt) override;

    void
    calc_force(particle_container_type& pcon, potential_type& pot) const override;

    real_type
    calc_energy(const particle_container_type& pcon,
                const potential_type& pot) const override;

    void set_spatial_partition(std::unique_ptr<spatial_partitioning_type>&& sp)
    {
        spatial_partition_ = std::forward<spatial_partitioning_type>(sp);
    }

#ifdef MJOLNIR_PARALLEL_THREAD
    void
    calc_force(particle_container_type& pcon, const potential_type& pot,
               const std::size_t num_threads) const override;
    real_type
    calc_energy(const particle_container_type& pcon, const potential_type& pot,
                const std::size_t num_threads) const override;

  private:

    void
    calc_force_(particle_container_type& pcon,
                std::size_t begin, const std::size_t end,
                const potential_type& pot) const;
    real_type
    calc_energy_(const particle_container_type& pcon,
                 std::size_t begin, const std::size_t end,
                 const potential_type& pot) const;
#endif

  private:
    std::unique_ptr<spatial_partitioning_type> spatial_partition_;
};

template<typename traitsT>
void GlobalDistanceInteraction<traitsT>::calc_force(
        particle_container_type& pcon, potential_type& pot) const
{
    spatial_partition_->update(pcon);
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partitioning_type::index_list const& partners =
            spatial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type       l   = length(rij);
            const coordinate_type f   = rij * (pot.derivative(i, j, l) / l);
            pcon[i].force += f;
            pcon[j].force -= f;
        }
    }
    return ;
}

template<typename traitsT>
typename GlobalDistanceInteraction<traitsT>::real_type
GlobalDistanceInteraction<traitsT>::calc_energy(
        const particle_container_type& pcon, const potential_type& pot) const
{
    real_type e = 0.0;
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        typename spatial_partitioning_type::index_list const& partners =
            spatial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type l = length(rij);
            e += pot.potential(i, j, l);
        }
    }
    return e;
}

template<typename traitsT>
inline void GlobalDistanceInteraction<traitsT>::initialize(
        const particle_container_type& pcon, const time_type dt)
{
    this->spatial_partition_->update(pcon, dt);
    return ;
}

#ifdef MJOLNIR_PARALLEL_THREAD

template<typename traitsT>
inline void GlobalDistanceInteraction<traitsT>::calc_force(
        particle_container_type& pcon, const potential_type& pot,
        const std::size_t num_threads) const
{
    spatial_partition_->update(pcon);
    std::vector<std::future<void>> futures(num_threads-1);
    const std::size_t particle_per_threads = pcon.size() / num_threads + 1;

    std::size_t begin = 0;
    std::size_t end = particle_per_threads;
    for(std::size_t t=0; t<num_threads-1; ++t)
    {
        futures.at(t) = std::async(std::launch::async,
            [this, &pcon, begin, end, &pot](){calc_force_(pcon, begin, end, pot);});
        begin += particle_per_threads;
        end   += particle_per_threads;
    }
    calc_force_(pcon, begin, pcon.size(), pot);

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        iter->get();

    return ;
}

template<typename traitsT>
inline typename GlobalDistanceInteraction<traitsT>::real_type
GlobalDistanceInteraction<traitsT>::calc_energy(
        const particle_container_type& pcon, const potential_type& pot,
        const std::size_t num_threads) const
{
    spatial_partition_->update(pcon);
    std::vector<std::future<real_type>> futures(num_threads-1);
    const std::size_t particle_per_threads = pcon.size() / num_threads + 1;

    std::size_t begin = 0;
    std::size_t end = particle_per_threads;
    for(std::size_t t=0; t<num_threads-1; ++t)
    {
        futures.at(t) = std::async(std::launch::async,
            [this, &pcon, begin, end, &pot](){return calc_energy_(pcon, begin, end, pot);});
        begin += particle_per_threads;
        end   += particle_per_threads;
    }
    real_type energy = calc_energy_(pcon, begin, pcon.size(), pot);

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
        energy += iter->get();

    return energy;
}

template<typename traitsT>
inline void GlobalDistanceInteraction<traitsT>::calc_force_(
        particle_container_type& pcon,
        std::size_t i, const std::size_t end, const potential_type& pot) const
{
    for(; i != end; ++i)
    {
        typename spatial_partitioning_type::index_list const& partners =
            spatial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type       l   = length(rij);
            const coordinate_type f   = rij * (pot.derivative(i, j, l) / l);

#ifdef MJOLNIR_PARALLEL_THREAD
    std::lock_guard<std::mutex> lock1(pcon[i].mtx);
    std::lock_guard<std::mutex> lock2(pcon[j].mtx);
#endif
            pcon[i].force += f;
            pcon[j].force -= f;
        }
    }
    return;
}

template<typename traitsT>
inline typename GlobalDistanceInteraction<traitsT>::real_type
GlobalDistanceInteraction<traitsT>::calc_energy_(
        const particle_container_type& pcon,
        std::size_t i, const std::size_t end, const potential_type& pot) const
{
    real_type e = 0.0;
    for(; i != end; ++i)
    {
        typename spatial_partitioning_type::index_list const& partners =
            spatial_partition_->partners(i);
        for(auto iter = partners.cbegin(); iter != partners.cend(); ++iter)
        {
            const std::size_t j = *iter;
            const coordinate_type rij = pcon[j].position - pcon[i].position;
            const real_type l = length(rij);
            e += pot.potential(i, j, l);
        }
    }
    return e;
}

#endif
} // mjolnir
#endif /* MJOLNIR_GLOBAL_DISTANCE_INTEARACTION */
