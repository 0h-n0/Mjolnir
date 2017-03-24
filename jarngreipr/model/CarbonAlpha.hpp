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
    CarbonAlpha(const residue_type& residue)
        : base_type(residue.atoms())
    {
        this->init_();
    }
    CarbonAlpha(const container_type& atoms)
        : base_type(atoms)
    {
        this->init_();
    }
    CarbonAlpha(const residue_type& residue, const std::string& name)
        : base_type(residue.atoms(), name)
    {
        this->init_();
    }
    CarbonAlpha(const container_type& atoms, const std::string& name)
        : base_type(atoms, name)
    {
        this->init_();
    }
    ~CarbonAlpha() = default;

    CarbonAlpha(const CarbonAlpha&) = default;
    CarbonAlpha(CarbonAlpha&&)      = default;
    CarbonAlpha& operator=(const CarbonAlpha&) = default;
    CarbonAlpha& operator=(CarbonAlpha&&)      = default;

    coordinate_type const& position() const override {return position_;}

    template<typename charT, typename ctrait>
    std::basic_ostream<charT, ctrait>&
    print(std::basic_ostream<charT, ctrait>& os) const override;

  private:

    void init_();

  private:

    coordinate_type position_;
};

template<typename traitsT>
void CarbonAlpha<traitsT>::init_()
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

    this->position_ =
        std::find_if(this->atoms_.cbegin(), this->atoms_.cend(),
                     finder)->position;
    return ;
}

template<typename traitsT>
template<typename charT, typename ctrait>
std::basic_ostream<charT, ctrait>&
CarbonAlpha<traitsT>::print(std::basic_ostream<charT, ctrait>& os) const
{
    auto finder = [](const atom_type& a){return a.atom_name == "CA";};
    os << *(std::find_if(this->atoms_.cbegin(), this->atoms_.cend(), finder));
    return os;
}



}//jarngreipr
#endif /*JARNGREIPR_CARBON_ALPHA*/
