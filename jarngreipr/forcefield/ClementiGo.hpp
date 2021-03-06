#ifndef JARNGREIPR_CLEMENTI_GO
#define JARNGREIPR_CLEMENTI_GO
#include <jarngreipr/forcefield/ForceFieldGenerator.hpp>
#include <jarngreipr/geometry/distance.hpp>
#include <jarngreipr/geometry/angle.hpp>
#include <jarngreipr/geometry/dihedral.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <iterator>
#include <iostream>
#include <vector>

namespace jarngreipr
{

template<typename realT>
class ClementiGo final : public IntraChainForceFieldGenerator<realT>
{
  public:
    typedef IntraChainForceFieldGenerator<realT> base_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::bead_type bead_type;

  public:

    ClementiGo(const toml::Table& para): parameters_(para){}
    ~ClementiGo() override = default;

    // generate parameters and write out to `ostrm`.
    void generate(toml::Table& out,
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

    bool check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const override;

  private:

    real_type min_distance_sq(const std::unique_ptr<bead_type>& bead1,
                              const std::unique_ptr<bead_type>& bead2) const
    {
        real_type min_dist_sq = std::numeric_limits<real_type>::max();
        for(const auto& atom1 : bead1->atoms())
        {
            if(atom1.element == " H"){continue;}
            for(const auto& atom2 : bead2->atoms())
            {
                if(atom2.element == " H") {continue;}
                const auto dsq = distance_sq(atom1.position, atom2.position);
                min_dist_sq = std::min(dsq, min_dist_sq);
            }
        }
        return min_dist_sq;
    }

  private:

    toml::Table parameters_;
};

template<typename realT>
void ClementiGo<realT>::generate(toml::Table& ff,
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    if(!this->check_beads_kind(chain))
    {
        std::cerr << "ClementiGo: stop generating forcefield..." << std::endl;
        return ;
    }

    if(ff.count("local") == 0)
    {
        ff["local"] = toml::Array();
    }

    /* bond-length */ {
        toml::Table bond_length;
        bond_length["interaction"] = toml::String("BondLength");
        bond_length["potential"]   = toml::String("Harmonic");
        bond_length["topology"]    = toml::String("bond");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 1; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            toml::Table para;
            para["indices"] = toml::value{i1, i2};
            para["eq"     ] = distance(bead1->position(), bead2->position());
            para["k"      ] = toml::get<toml::Float>(mjolnir::toml_value_at(
                    this->parameters_, "coef_bond", "[parameter.ClementiGo]"));
            params.push_back(std::move(para));
        }
        bond_length["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_length));
    }
    /* bond-angle */{
        toml::Table bond_angle;
        bond_angle["interaction"] = toml::String("BondAngle");
        bond_angle["potential"]   = toml::String("Harmonic");
        bond_angle["topology"]    = toml::String("none");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 2; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const auto& bead3 = chain.at(i+2);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            const std::size_t i3 = bead3->index();

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3};
            para["eq"     ] = angle(bead1->position(), bead2->position(),
                                    bead3->position());
            para["k"      ] = toml::get<toml::Float>(mjolnir::toml_value_at(
                    this->parameters_, "coef_angle", "[parameter.ClementiGo]"));
            params.push_back(std::move(para));
        }
        bond_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(bond_angle));
    }
    /* dihedral-angle */{
        toml::Table dihd_angle;
        dihd_angle["interaction"] = toml::String("DihedralAngle");
        dihd_angle["potential"]   = toml::String("ClementiDihedral");
        dihd_angle["topology"]    = toml::String("none");

        toml::Array params;
        for(std::size_t i=0, sz = chain.size() - 3; i<sz; ++i)
        {
            const auto& bead1 = chain.at(i);
            const auto& bead2 = chain.at(i+1);
            const auto& bead3 = chain.at(i+2);
            const auto& bead4 = chain.at(i+3);
            const std::size_t i1 = bead1->index();
            const std::size_t i2 = bead2->index();
            const std::size_t i3 = bead3->index();
            const std::size_t i4 = bead4->index();

            toml::Table para;
            para["indices"] = toml::value{i1, i2, i3, i4};
            para["eq"     ] = dihedral_angle(bead1->position(),
                    bead2->position(), bead3->position(), bead4->position());
            para["k1"     ] = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, "coef_dihedral_1", "[parameter.ClementiGo]"));
            para["k3"     ] = toml::get<toml::Float>(mjolnir::toml_value_at(
                this->parameters_, "coef_dihedral_3", "[parameter.ClementiGo]"));
            params.push_back(std::move(para));
        }
        dihd_angle["parameters"] = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(dihd_angle));
    }

    /* intra-chain-go-contacts */{
        const toml::Float threshold = toml::get<toml::Float>(
            mjolnir::toml_value_at(this->parameters_, "contact_threshold",
                                   "[parameter.ClementiGo]"));
        const real_type th2 = threshold * threshold;

        toml::Table go_contact;
        go_contact["interaction"] = toml::String("BondLength");
        go_contact["potential"  ] = toml::String("Go1012Contact");
        go_contact["topology"   ] = toml::String("contact");

        toml::Array params;
        for(std::size_t i=0, sz_i = chain.size()-4; i<sz_i; ++i)
        {
            for(std::size_t j=i+4, sz_j = chain.size(); j<sz_j; ++j)
            {
                if(this->min_distance_sq(chain.at(i), chain.at(j)) < th2)
                {
                    const auto& bead1 = chain.at(i);
                    const auto& bead2 = chain.at(j);
                    const std::size_t i1 = bead1->index();
                    const std::size_t i2 = bead2->index();

                    toml::Table para;
                    para["indices"] = toml::value{i1, i2};
                    para["eq"     ] = distance(bead1->position(), bead2->position());
                    para["k"      ] = toml::get<toml::Float>(mjolnir::toml_value_at(
                        this->parameters_, "coef_contact", "[parameter.ClementiGo]"));
                    params.push_back(std::move(para));
                }
            }
        }
        go_contact["parameters"]  = std::move(params);
        ff["local"].cast<toml::value_t::Array>().push_back(std::move(go_contact));
    }
    return;
}

template<typename realT>
bool ClementiGo<realT>::check_beads_kind(
        const std::vector<std::unique_ptr<bead_type>>& chain) const
{
    for(const auto& bead : chain)
    {
        if(bead->kind() != "CarbonAlpha")
        {
            std::cerr << "ClementiGo: invalid coarse-grained bead kind: "
                      << bead->kind() << '\n';
            std::cerr << "it allows only CarbonAlpha beads.\n";
            return false;
        }
    }
    return true;
}

}//jarngreipr
#endif /* JARNGREIPR_CLEMENTI_GO */
