#ifndef MJOLNIR_LOCAL_FORCE_FIELD
#define MJOLNIR_LOCAL_FORCE_FIELD
#include "BondLengthInteraction.hpp"
#include "BondAngleInteraction.hpp"
#include "DihedralAngleInteraction.hpp"
#include "LocalPotentialBase.hpp"
#include "ParticleContainer.hpp"
#include <utility>
#include <vector>
#include <array>
#include <memory>

#ifdef MJOLNIR_PARALLEL_THREAD
#include <thread>
#include <future>
#include <mutex>
#include <atomic>
#endif

#ifdef MJOLNIR_DUMP_TIME
#include <iostream>
#include <mjolnir/util/time_stamp.hpp>
#endif // MJOLNIR_DUMP_TIME

namespace mjolnir
{

template<typename traitsT>
class LocalForceField
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::time_type time_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Particle<coordinate_type> particle_type;
    typedef ParticleContainer<traits_type> particle_container_type;
    typedef LocalPotentialBase<traits_type> potential_base;
    typedef std::unique_ptr<potential_base> potential_ptr;
    typedef std::array<std::size_t, 2> pid2_type;
    typedef std::array<std::size_t, 3> pid3_type;
    typedef std::array<std::size_t, 4> pid4_type;

  private:
    struct bond_potential_type;
    struct angle_potential_type;
    struct dihd_potential_type;
    typedef std::vector<bond_potential_type>  bond_container;
    typedef std::vector<angle_potential_type> angle_container;
    typedef std::vector<dihd_potential_type>  dihd_container;
  
  public:

    LocalForceField() = default;

    ~LocalForceField() = default;
    LocalForceField(LocalForceField const&) = delete;
    LocalForceField(LocalForceField&&)      = default;
    LocalForceField& operator=(LocalForceField const&) = delete;
    LocalForceField& operator=(LocalForceField&&)      = default;

    void emplace_bond(std::size_t i, std::size_t j, potential_ptr&& pot);
    void emplace_angle(std::size_t i, std::size_t j, std::size_t k,
            potential_ptr&& pot);
    void emplace_dihedral(std::size_t i, std::size_t j, std::size_t k,
            std::size_t l, potential_ptr&& pot);

    void      calc_force(particle_container_type& pcon) const;
    real_type calc_energy(const particle_container_type& pcon) const;

#ifdef MJOLNIR_PARALLEL_THREAD
    void calc_force(particle_container_type& pcon,
                    const std::size_t num_threads) const;
    real_type calc_energy(const particle_container_type& pcon,
                          const std::size_t num_threads) const;
#endif

    void
    calc_bond_force(particle_container_type& pcon,
                    typename bond_container::const_iterator iter,
                    const typename bond_container::const_iterator end) const;
    void
    calc_angle_force(particle_container_type& pcon,
                     typename angle_container::const_iterator iter,
                     const typename angle_container::const_iterator end) const;
    void
    calc_dihd_force(particle_container_type& pcon,
                    typename dihd_container::const_iterator iter,
                    const typename dihd_container::const_iterator end) const;

    real_type
    calc_bond_energy(const particle_container_type& pcon,
                     typename bond_container::const_iterator iter,
                     const typename bond_container::const_iterator end) const;
    real_type
    calc_angle_energy(const particle_container_type& pcon,
                      typename angle_container::const_iterator iter,
                      const typename angle_container::const_iterator end) const;
    real_type
    calc_dihd_energy(const particle_container_type& pcon,
                     typename dihd_container::const_iterator iter,
                     const typename dihd_container::const_iterator end) const;

  private:

    BondLengthInteraction<traitsT>    bond;
    BondAngleInteraction<traitsT>     angle;
    DihedralAngleInteraction<traitsT> dihd;

    struct bond_potential_type
    {
        bond_potential_type() = default;
        bond_potential_type(std::size_t i_, std::size_t j_, potential_ptr&& pot_)
           : i(i_), j(j_), pot(std::forward<potential_ptr>(pot_))
        {}
        ~bond_potential_type() = default;
        bond_potential_type(const bond_potential_type&) = delete;
        bond_potential_type(bond_potential_type&&)      = default;
        bond_potential_type& operator=(const bond_potential_type&) = delete;
        bond_potential_type& operator=(bond_potential_type&&)      = default;

        std::size_t i, j;
        potential_ptr pot;
    };

    struct angle_potential_type
    {
        angle_potential_type() = default;
        angle_potential_type(std::size_t i_, std::size_t j_, std::size_t k_,
                             potential_ptr&& pot_)
           : i(i_), j(j_), k(k_), pot(std::forward<potential_ptr>(pot_))
        {}
        ~angle_potential_type() = default;
        angle_potential_type(const angle_potential_type&) = delete;
        angle_potential_type(angle_potential_type&&)      = default;
        angle_potential_type& operator=(const angle_potential_type&) = delete;
        angle_potential_type& operator=(angle_potential_type&&)      = default;

        std::size_t i, j, k;
        potential_ptr pot;
    };

    struct dihd_potential_type
    {
        dihd_potential_type() = default;
        dihd_potential_type(std::size_t i_, std::size_t j_, std::size_t k_,
                            std::size_t l_, potential_ptr&& pot_)
           : i(i_), j(j_), k(k_), l(l_), pot(std::forward<potential_ptr>(pot_))
        {}
        ~dihd_potential_type() = default;
        dihd_potential_type(const dihd_potential_type&) = delete;
        dihd_potential_type(dihd_potential_type&&)      = default;
        dihd_potential_type& operator=(const dihd_potential_type&) = delete;
        dihd_potential_type& operator=(dihd_potential_type&&)      = default;

        std::size_t i, j, k, l;
        potential_ptr pot;
    };

    bond_container  bond_potentials;
    angle_container angle_potentials;
    dihd_container  dihd_potentials;
};

template<typename traitsT>
void LocalForceField<traitsT>::emplace_bond(
        std::size_t i, std::size_t j, potential_ptr&& pot)
{
    bond_potentials.emplace_back(i, j, std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::emplace_angle(
        std::size_t i, std::size_t j, std::size_t k, potential_ptr&& pot)
{
    angle_potentials.emplace_back(i, j, k, std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::emplace_dihedral(
        std::size_t i, std::size_t j, std::size_t k, std::size_t l, potential_ptr&& pot)
{
    dihd_potentials.emplace_back(i, j, k, l, std::forward<potential_ptr>(pot));
    return ;
}

template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon) const
{
    calc_bond_force(pcon, bond_potentials.cbegin(),  bond_potentials.cend());
    calc_angle_force(pcon, angle_potentials.cbegin(), angle_potentials.cend());
    calc_dihd_force(pcon, dihd_potentials.cbegin(),  dihd_potentials.cend());
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon) const
{
    return calc_bond_energy(
            pcon, bond_potentials.cbegin(),  bond_potentials.cend()) +
        calc_angle_energy(
            pcon, angle_potentials.cbegin(), angle_potentials.cend())+
        calc_dihd_energy(
            pcon, dihd_potentials.cbegin(),  dihd_potentials.cend());
}

template<typename traitsT>
inline void LocalForceField<traitsT>::calc_bond_force(
        particle_container_type& pcon,
        typename bond_container::const_iterator iter,
        const typename bond_container::const_iterator end) const
{
    for(; iter != end; ++iter)
        bond.calc_force(pcon[iter->i], pcon[iter->j], *(iter->pot));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::calc_angle_force(
        particle_container_type& pcon,
        typename angle_container::const_iterator iter,
        const typename angle_container::const_iterator end) const
{
    for(; iter != end; ++iter)
        angle.calc_force(pcon[iter->i], pcon[iter->j], pcon[iter->k],
                         *(iter->pot));
    return;
}

template<typename traitsT>
inline void LocalForceField<traitsT>::calc_dihd_force(
        particle_container_type& pcon,
        typename dihd_container::const_iterator iter,
        const typename dihd_container::const_iterator end) const
{
    for(; iter != end; ++iter)
        dihd.calc_force(pcon[iter->i], pcon[iter->j], pcon[iter->k],
                        pcon[iter->l], *(iter->pot));
    return;
}

template<typename traitsT>
inline typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_bond_energy(
        const particle_container_type& pcon,
        typename bond_container::const_iterator iter,
        const typename bond_container::const_iterator end) const
{
    real_type energy = 0.;
    for(; iter != end; ++iter)
        energy += bond.calc_energy(
                pcon[iter->i], pcon[iter->j], *(iter->pot));
    return energy;
}

template<typename traitsT>
inline typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_angle_energy(
        const particle_container_type& pcon,
        typename angle_container::const_iterator iter,
        const typename angle_container::const_iterator end) const
{
    real_type energy = 0.;
    for(; iter != end; ++iter)
        energy += angle.calc_energy(
                pcon[iter->i], pcon[iter->j], pcon[iter->k], *(iter->pot));
    return energy;
}

template<typename traitsT>
inline typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_dihd_energy(
        const particle_container_type& pcon,
        typename dihd_container::const_iterator iter,
        const typename dihd_container::const_iterator end) const
{
    real_type energy = 0.;
    for(; iter != end; ++iter)
        energy += dihd.calc_energy(
                pcon[iter->i], pcon[iter->j], pcon[iter->k], pcon[iter->l],
                *(iter->pot));
    return energy;
}

// parallel impl----------------------------------------------------------------
#ifdef MJOLNIR_PARALLEL_THREAD
template<typename traitsT>
void LocalForceField<traitsT>::calc_force(particle_container_type& pcon,
        const std::size_t num_threads) const
{
#ifdef MJOLNIR_DUMP_TIME
    std::vector<std::future<timestamp<void>>> futures(num_threads-1);
#else
    std::vector<std::future<void>> futures(num_threads-1);
#endif // dump time

// bond ------------------------------------------------------------------------
    {
#ifdef MJOLNIR_DUMP_TIME
    auto start = std::chrono::system_clock::now();
    std::cerr << "calculation of bond starts at : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      start.time_since_epoch()).count() << std::endl;
#endif // dump time

    const std::size_t tasks_per_thread =
        bond_potentials.size() / num_threads + 1;

    auto begin = bond_potentials.cbegin();
    auto end   = bond_potentials.cbegin() + tasks_per_thread;

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread creating loop starts at: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto b = this->bond;
        futures.at(t) = std::async(std::launch::async,
            [b, &pcon, begin, end]()
            {
#ifdef MJOLNIR_DUMP_TIME
                timestamp<void> ts; ts.start();
#endif // dump time
                for(auto iter = begin; iter != end; ++iter)
                    b.calc_force(pcon.at(iter->i), pcon.at(iter->j),
                                 *(iter->pot));
#ifdef MJOLNIR_DUMP_TIME
                ts.stop();
                return ts;
#else
                return ;
#endif // dump time
            });
        begin += tasks_per_thread;
        end   += tasks_per_thread;
    }

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread creating loop ends at  : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    // last one is executed "this" thread.
    end = bond_potentials.cend();
    for(; begin != end; ++begin)
        bond.calc_force(pcon.at(begin->i), pcon.at(begin->j), *(begin->pot));

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
    {
#ifdef MJOLNIR_DUMP_TIME
        auto ts = iter->get();
        std::cerr << "-----------------------------------" << std::endl;
        std::cerr << "thread starts at              : "
                  << ts.start_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "thread ends at                : "
                  << ts.stop_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "elapsed per thread            : "
                  << ts.elapsed<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
#else  // normal one
        iter->wait();
#endif // dump time
    }


#ifdef MJOLNIR_DUMP_TIME
    auto end_t = std::chrono::system_clock::now();
    std::cerr << "-----------------------------------" << std::endl;
    std::cerr << "calc of bond length ends at   : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t.time_since_epoch()).count() << " [usec]" << std::endl;
    std::cerr << "elapsed for bond length/step  : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t - start).count() << " [usec]" << std::endl;
    std::cerr << "==================================="<< std::endl;
    std::cerr << std::endl;
#endif // dump time
    }
// angle -----------------------------------------------------------------------
    {
#ifdef MJOLNIR_DUMP_TIME
    auto start = std::chrono::system_clock::now();
    std::cerr << "calculation of angle starts at: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      start.time_since_epoch()).count() << std::endl;
#endif // dump time

    const std::size_t tasks_per_thread =
        angle_potentials.size() / num_threads + 1;

    auto begin = angle_potentials.cbegin();
    auto end   = angle_potentials.cbegin() + tasks_per_thread;

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread creating loop starts at: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto a = this->angle;
        futures.at(t) = std::async(std::launch::async,
            [a, &pcon, begin, end]()
            {
#ifdef MJOLNIR_DUMP_TIME
                timestamp<void> ts; ts.start();
#endif // dump time
                for(auto iter = begin; iter != end; ++iter)
                    a.calc_force(pcon.at(iter->i), pcon.at(iter->j),
                                 pcon.at(iter->k), *(iter->pot));
#ifdef MJOLNIR_DUMP_TIME
                ts.stop();
                return ts;
#else
                return ;
#endif // dump time
            });
        begin += tasks_per_thread;
        end   += tasks_per_thread;
    }

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread making loop ends at    : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    // last one is executed "this" thread.
    end = angle_potentials.cend();
    for(; begin != end; ++begin)
        angle.calc_force(pcon.at(begin->i), pcon.at(begin->j),
                         pcon.at(begin->k), *(begin->pot));
    

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
    {
#ifdef MJOLNIR_DUMP_TIME
        auto ts = iter->get();
        std::cerr << "-----------------------------------" << std::endl;
        std::cerr << "thread starts at              : "
                  << ts.start_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "thread ends at                : "
                  << ts.stop_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "elapsed per thread            : "
                  << ts.elapsed<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
#else
        iter->wait();
#endif // dump time
    }

#ifdef MJOLNIR_DUMP_TIME
    auto end_t = std::chrono::system_clock::now();
    std::cerr << "-----------------------------------" << std::endl;
    std::cerr << "calc of bond angle ends at    : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t.time_since_epoch()).count() << " [usec]" << std::endl;
    std::cerr << "elapsed for bond angle/step   : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t - start).count() << " [usec]" << std::endl;
    std::cerr << "==================================="<< std::endl;
    std::cerr << std::endl;
#endif // dump time
    }
// dihd ------------------------------------------------------------------------
    {
#ifdef MJOLNIR_DUMP_TIME
    auto start = std::chrono::system_clock::now();
    std::cerr << "calculation of dihd starts at : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      start.time_since_epoch()).count() << std::endl;
#endif // dump time

    const std::size_t tasks_per_thread =
        dihd_potentials.size() / num_threads + 1;
    auto begin = dihd_potentials.cbegin();
    auto end   = dihd_potentials.cbegin() + tasks_per_thread;

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread creating loop starts at: "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto d = this->dihd;
        futures.at(t) = std::async(std::launch::async,
            [d, &pcon, begin, end]()
            {
#ifdef MJOLNIR_DUMP_TIME
                timestamp<void> ts; ts.start();
#endif
                for(auto iter = begin; iter != end; ++iter)
                   d.calc_force(pcon.at(iter->i), pcon.at(iter->j),
                                pcon.at(iter->k), pcon.at(iter->l),
                                *(iter->pot));
#ifdef MJOLNIR_DUMP_TIME
                ts.stop();
                return ts;
#else
                return ;
#endif
            });
        begin += tasks_per_thread;
        end   += tasks_per_thread;
    }

#ifdef MJOLNIR_DUMP_TIME
    std::cerr << "thread creating loop ends at  : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                     (std::chrono::system_clock::now()).time_since_epoch()).count()
              << std::endl;
#endif // dump time

    // last one is executed "this" thread.
    end = dihd_potentials.cend();
    for(; begin != end; ++begin)
       dihd.calc_force(pcon.at(begin->i), pcon.at(begin->j),
                       pcon.at(begin->k), pcon.at(begin->l),
                       *(begin->pot));

    for(auto iter = futures.begin(); iter != futures.end(); ++iter)
    {
#ifdef MJOLNIR_DUMP_TIME
        auto ts = iter->get();
        std::cerr << "-----------------------------------" << std::endl;
        std::cerr << "thread starts at              : "
                  << ts.start_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "thread ends at                : "
                  << ts.stop_at<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
        std::cerr << "elapsed per thread            : "
                  << ts.elapsed<std::chrono::microseconds>().count()
                  << " [usec]" << std::endl;
#else
        iter->wait();
#endif
    }

#ifdef MJOLNIR_DUMP_TIME
    auto end_t = std::chrono::system_clock::now();
    std::cerr << "-----------------------------------" << std::endl;
    std::cerr << "calc of dihd angle ends at    :"
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t.time_since_epoch()).count()
              << std::endl;
 
    std::cerr << "elapsed for dihd angle/step   : "
              << std::chrono::duration_cast<std::chrono::microseconds>(
                      end_t - start).count() << " [usec]" << std::endl;
    std::cerr << "==================================="<< std::endl;
    std::cerr << std::endl;
#endif
    }
    return;
}

template<typename traitsT>
typename LocalForceField<traitsT>::real_type
LocalForceField<traitsT>::calc_energy(const particle_container_type& pcon,
        const std::size_t num_threads) const
{
    real_type energy = 0.;
    std::vector<std::future<real_type>> futures(num_threads-1);

// bond ------------------------------------------------------------------------
    const std::size_t bsize = bond_potentials.size();
    const std::size_t bonds_per_thread = bsize / num_threads + 1;
    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto begin = bond_potentials.cbegin() + t * bonds_per_thread;
        auto end   = bond_potentials.cbegin() + (t+1) * bonds_per_thread;
        futures.at(t) = std::async(std::launch::async,
            [this, &pcon, begin, end](){return calc_bond_energy(pcon, begin, end);});
    }
    {// last one is executed "this" thread.
        auto begin = bond_potentials.cbegin() + (num_threads - 1) * bonds_per_thread;
        auto end = bond_potentials.cend();
        energy += calc_bond_energy(pcon, begin, end);
    }
    for(std::size_t t = 0; t < num_threads-1; ++t)
        energy += futures.at(t).get();

// angle -----------------------------------------------------------------------
    const std::size_t asize = angle_potentials.size();
    const std::size_t angles_per_thread = asize / num_threads + 1;
    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto begin = angle_potentials.cbegin() + t     * angles_per_thread;
        auto end   = angle_potentials.cbegin() + (t+1) * angles_per_thread;
        futures.at(t) = std::async(std::launch::async,
            [this, &pcon, begin, end](){return calc_angle_energy(pcon, begin, end);});
    }
    {// last one is executed "this" thread.
        auto begin = angle_potentials.cbegin() + (num_threads - 1) * angles_per_thread;
        auto end   = angle_potentials.cend();
        energy += calc_angle_energy(pcon, begin, end);
    }
    for(std::size_t t = 0; t < num_threads-1; ++t)
        energy += futures.at(t).get();

// dihd ------------------------------------------------------------------------
    const std::size_t dsize = dihd_potentials.size();
    const std::size_t dihds_per_thread = dsize / num_threads + 1;
    for(std::size_t t=0; t < num_threads-1; ++t)
    {
        auto begin = dihd_potentials.cbegin() + t     * dihds_per_thread;
        auto end   = dihd_potentials.cbegin() + (t+1) * dihds_per_thread;
        futures.at(t) = std::async(std::launch::async,
            [this, &pcon, begin, end](){return calc_dihd_energy(pcon, begin, end);});
    }
    {// last one is executed "this" thread.
        auto begin = dihd_potentials.cbegin() + (num_threads - 1) * dihds_per_thread;
        auto end   = dihd_potentials.cend();
        energy += calc_dihd_energy(pcon, begin, end);
    }
    for(std::size_t t = 0; t < num_threads-1; ++t)
        energy += futures.at(t).get();

    return energy;
}
#endif

} // mjolnir
#endif /* MJOLNIR_LOCAL_FORCE_FIELD */
