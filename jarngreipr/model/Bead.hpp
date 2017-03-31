#ifndef JARNGREIPR_BEAD
#define JARNGREIPR_BEAD
#include <jarngreipr/pdb/PDBAtom.hpp>
#include <ostream>
#include <vector>
#include <string>

namespace jarngreipr
{

template<typename traitsT>
class Bead
{
  public:

    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBAtom<traits_type> atom_type;
    typedef std::vector<atom_type> container_type;

  public:

    Bead() = default;
    virtual ~Bead() = default;
    Bead(const Bead&) = default;
    Bead(Bead&&)      = default;
    Bead& operator=(const Bead&) = default;
    Bead& operator=(Bead&&)      = default;

    explicit Bead(const container_type& atoms) : atoms_(atoms){}
    explicit Bead(const std::string& name) : name_(name){}
    Bead(const container_type& atoms, const std::string& name)
        : atoms_(atoms), name_(name)
    {}
    Bead(container_type&& atoms, std::string&& name)
        : atoms_(std::forward<container_type>(atoms)),
          name_(std::forward<std::string>(name))
    {}

    virtual coordinate_type position() const = 0;

    virtual std::ostream& print(std::ostream& os) const = 0;

    container_type const& atoms() const {return atoms_;}
    container_type &      atoms()       {return atoms_;}

    std::string const& name() const {return name_;}
    std::string &      name()       {return name_;}

    std::size_t const& bead_id() const {return bead_id_;}
    std::size_t &      bead_id()       {return bead_id_;}

  protected:

    std::size_t     bead_id_;
    std::string     name_;
    container_type  atoms_;
};

}//jarngreipr
#endif /* JARNGREIPR_BEAD */
