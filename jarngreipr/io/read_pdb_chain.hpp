#ifndef JARNGREIPR_IO_READ_PDB
#define JARNGREIPR_IO_READ_PDB
#include "PDBReader.hpp"

namespace jarngreipr
{

template<typename traitsT>
PDBChain<traits_type>
read_pdb_chain(const std::string& fname, const std::string& chainID)
{
    PDBReader<traitsT> reader;
    for(auto&& chain : reader.parse(reader.read(fname)))
    {
        if(chain.chain_id() == chainID) return chain;
    }
    throw std::invalid_argument(std::string("no chain ") + chainID +
                                std::string(" in ") + fname);
}

}//jarngreipr
#endif// JARNGREIPR_IO_READ_PDB
