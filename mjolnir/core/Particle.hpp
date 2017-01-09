#ifndef MJOLNIR_PARTICLE
#define MJOLNIR_PARTICLE
#include <mjolnir/util/scalar_type_of.hpp>

#ifdef MJOLNIR_PARALLEL_THREAD
#include <mutex>
#endif

namespace mjolnir
{

template<typename coordT>
struct Particle
{
    scalar_type<coordT> mass;
    coordT position;
    coordT velocity;
    coordT force;

#ifdef MJOLNIR_PARALLEL_THREAD
    Particle() = default;
    Particle(const scalar_type<coordT> m, const coordT& pos,
             const coordT& vel, const coordT& f)
        : mass(m), position(pos), velocity(vel), force(f), mtx()
    {}
    ~Particle() = default;
    Particle(const Particle& p) = delete;
    Particle(Particle&& p)      = default;
    Particle& operator=(const Particle& p) = delete;
    Particle& operator=(Particle&& p)      = default;
    std::mutex mtx;
#endif

};

} // mjolnir
#endif /* MJOLNIR_PARTICLE */
