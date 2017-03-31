#ifndef JARNGREIPR_CARBON_ALPHA_MODEL
#define JARNGREIPR_CARBON_ALPHA_MODEL
#include "CarbonAlpha.hpp"
#include "CGChain.hpp"
#include "Model.hpp"

namespace jarngreipr
{

template<typename traitsT>
class CarbonAlphaModel : public Model<traitsT>
{
  public:
    typedef traitsT traits_type;
    typedef Model<traitsT> base_type;
    typedef typename base_type::real_type           real_type;
    typedef typename base_type::coordinate_type     coordinate_type;
    typedef typename base_type::chain_type          chain_type;
    typedef typename base_type::bead_type           bead_type;
    typedef typename base_type::bead_ptr            bead_ptr;
    typedef typename base_type::bead_container_type bead_container_type;

  public:

    CarbonAlphaModel() = default;
    ~CarbonAlphaModel() = default;

    static
    CGChain<traitsT> apply(const chain_type& chain);
};

template<typename traitsT>
CGChain<traitsT> CarbonAlphaModel<traitsT>::apply(const chain_type& chain)
{
    CGChain<traitsT> beads;
    beads.chainID() = chain.chain_id();
    for(auto iter = chain.cbegin(); iter != chain.cend(); ++iter)
        beads.emplace_back(new CarbonAlpha<traitsT>(*iter, "CA"));
    return beads;
}

}//jarngreipr
#endif// JARNGREIPR_CARBON_ALPHA_MODEL
