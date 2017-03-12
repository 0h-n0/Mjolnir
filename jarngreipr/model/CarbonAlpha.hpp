#ifndef JARNGREIPR_CARBON_ALPHA
#define JARNGREIPR_CARBON_ALPHA
#include "Bead.hpp"
#include <algorithm>
#include <stdexcept>
#include <string>

namespace jarngreipr
{

/*! @brief carbon alpha 1 beads model */
template<typename traitsT>
class CarbonAlpha : public Bead<traitsT>
{
  public:

    typedef Bead<traitsT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::container_type container_type;
    typedef typename base_type::atom_type atom_type;
    typedef PDBResidue<traitsT> residue_type;

  public:

    CarbonAlpha() = default;
    CarbonAlpha(const residue_type& residue) : base_type(residue.atoms()){}
    CarbonAlpha(const container_type& atoms) : base_type(atoms){}
    CarbonAlpha(const residue_type& residue, const std::string& name)
        : base_type(residue.atoms(), name)
    {}
    CarbonAlpha(const container_type& atoms, const std::string& name)
        : base_type(atoms, name)
    {}
    ~CarbonAlpha() = default;

    coordinate_type position() const override;
};

template<typename traitsT>
typename CarbonAlpha<traitsT>::coordinate_type
CarbonAlpha<traitsT>::position() const
{
    auto finder = [](const atom_type& a){return a.atom_name == "CA";};

    const std::size_t num_ca =
        std::count_if(this->atoms_.cbegin(), this->atoms_.cend(), finder);

    if(num_ca == 0)
        throw std::runtime_error(
                "jarngreipr::Ca::position: no C-alpha atom in this residue");
    else if(num_ca > 1)
        throw std::out_of_range(
                "jarngreipr::Ca::position: too many C-alphas in the bead");

    const auto ca =
        std::find_if(this->atoms_.cbegin(), this->atoms_.cend(), finder);
    return ca->position;
}

}//jarngreipr
#endif /*JARNGREIPR_CARBON_ALPHA*/
