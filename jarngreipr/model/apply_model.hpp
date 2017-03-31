#ifndef JARNGREIPR_APPLY_MODEL
#define JARNGREIPR_APPLY_MODEL
#include "CGChain.hpp"
#include "CarbonAlphaModel.hpp"

namespace jarngreipr
{

template<typename traitsT>
CGChain<traitsT>
apply_model(const std::string& model, const PDBChain<traitsT>& chain)
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
