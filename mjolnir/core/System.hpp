#ifndef MJOLNIR_SYSTEM
#define MJOLNIR_SYSTEM
#include "Particle.hpp"
#include <vector>

namespace mjolnir
{

template<typename traitsT>
class System
{
  public:
    typedef traitsT  traits_type;
    typedef typename traits_type::real_type         real_type;
    typedef typename traits_type::coordinate_type   coordinate_type;
    typedef typename traits_type::boundary_type     boundary_type;
    typedef Particle<coordinate_type>               particle_type;
    typedef std::vector<particle_type>              container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

    System()  = default;
    ~System() = default;

    System(container_type&& pcon, boundary_type&& bound)
        : boundary_(bound), particles_(pcon), temperature_(0.), pressure_(0.),
          ionic_strength_(0.)
    {}

    std::size_t size() const noexcept {return particles_.size();}

    real_type& pressure()             noexcept {return pressure_;}
    real_type  pressure()       const noexcept {return pressure_;}
    real_type& temperature()          noexcept {return temperature_;}
    real_type  temperature()    const noexcept {return temperature_;}
    real_type& ionic_strength()       noexcept {return ionic_strength_;}
    real_type  ionic_strength() const noexcept {return ionic_strength_;}

    coordinate_type adjust_direction(coordinate_type dr) const noexcept
    {return boundary_.adjust_direction(dr);}
    coordinate_type  adjust_position(coordinate_type dr) const noexcept
    {return boundary_.adjust_position(dr);}

    particle_type &      operator[](std::size_t i)       noexcept
    {return particles_[i];}
    particle_type const& operator[](std::size_t i) const noexcept
    {return particles_[i];}
    particle_type &      at(std::size_t i)       noexcept
    {return particles_.at(i);}
    particle_type const& at(std::size_t i) const noexcept
    {return particles_.at(i);}

    iterator       begin()        noexcept {return particles_.begin();}
    iterator       end()          noexcept {return particles_.end();}
    const_iterator begin()  const noexcept {return particles_.cbegin();}
    const_iterator end()    const noexcept {return particles_.cend();}
    const_iterator cbegin() const noexcept {return particles_.cbegin();}
    const_iterator cend()   const noexcept {return particles_.cend();}

    boundary_type const& boundary() const noexcept {return boundary_;}

    real_type& max_speed()       noexcept {return max_speed_;}
    real_type  max_speed() const noexcept {return max_speed_;}

  private:

    real_type      temperature_;
    real_type      pressure_;
    real_type      ionic_strength_;
    real_type      max_speed_;
    boundary_type  boundary_;
    container_type particles_;
};




} // mjolnir
#endif// MJOLNIR_SYSTEM
