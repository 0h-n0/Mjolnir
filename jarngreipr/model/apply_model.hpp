#ifndef JARNGREIPR_APPLY_MODEL
#define JARNGREIPR_APPLY_MODEL
#include "CarbonAlphaModel.hpp"

namespace jarngreipr
{

template<typename traitsT>
typename Model<traitsT>::bead_container_type
apply_model(const std::string& model,
            const typename Model<traitsT>::chain_type& chain)
{
    if(model == "CarbonAlpha")
    {
        return CarbonAlphaModel<traitsT>::apply(chain);
    }
    else
    {
        throw std::invalid_argument("Unknown Molecular Model: " + model);
    }
}

}//jarngreipr
#endif// JARNGREIPR_APPLY_MODEL
