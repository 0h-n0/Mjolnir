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
    std::mutex mtx;
#endif

};

template<typename coordT>
constexpr inline Particle<coordT>
make_particle(const scalar_type<coordT> mass, const coordT& pos,
              const coordT& vel, const coordT& f)
{
    return Particle<coordT>{mass, pos, vel, f};
}


} // mjolnir
#endif /* MJOLNIR_PARTICLE */
