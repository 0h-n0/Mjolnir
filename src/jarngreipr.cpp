#include <jarngreipr/pdb/unpack_chainID.hpp>
#include <jarngreipr/pdb/read_pdb_chain.hpp>
#include <jarngreipr/model/apply_model.hpp>

#include <mjolnir/core/DefaultTraits.hpp>
#include <mjolnir/core/RandomNumberGenerator.hpp>
#include <toml/toml.hpp>
#include <map>

typedef mjolnir::DefaultTraits traits;

int main(int argc, char **argv)
{
    if(argc != 2)
    {
        std::cerr << "./jarngreipr [file.toml]" << std::endl;
        return 1;
    }

    std::ifstream tomlfile(argv[1]);
    if(!tomlfile.good())
    {
        std::cerr << "file open error: " << argv[1] << std::endl;
        return 1;
    }

    const auto input = toml::parse(tomlfile);
    const auto& general = toml::get<toml::Table>(input.at("general"));
    const std::string file_name =
        toml::get<toml::String>(general.at("output_path")) +
        toml::get<toml::String>(general.at("file_name"));
    const std::int64_t seed =
        toml::get<toml::Integer>(general.at("seed"));
    const std::string parameter_path =
        toml::get<toml::String>(general.at("parameter_path"));

    std::ifstream parameter_file(parameter_path + std::string("parameter.toml"));
    if(!parameter_file.good())
    {
        std::cerr << "file open error: " << parameter_path << "parameter.toml"
                  << std::endl;
        return 1;
    }
    const auto params = toml::parse(parameter_file);

    // generate coarse-grained structures
    const auto& structures = toml::get<toml::Table>(input.at("structures"));

    std::vector<jarngreipr::PDBChain<traits>> aa_chains;
    std::vector<jarngreipr::CGChain<traits>>  cg_chains;
    for(auto& item : structures)
    {
        const auto& ref_model = toml::get<toml::Table>(item.second);
        const auto& reference = toml::get<toml::String>(ref_model.at("reference"));
        const auto& model     = toml::get<toml::String>(ref_model.at("model"));
        for(auto chain : jarngreipr::unpack_chainID(item.first))
        {// chain == "A", "B", ..."Z", "AA", "AB", ...
            aa_chains.emplace_back(jarngreipr::read_pdb_chain<traits>(reference, chain));
            cg_chains.emplace_back(jarngreipr::apply_model<traits>(model, aa_chains.back()));
        }
    }
    {// update index
        std::size_t bead_idx    = 0;
        for(auto& cg_chain : cg_chains)
        {
            for(auto& cg_bead : cg_chain)
            {
                for(auto& atom : cg_bead->atoms())
                {
                    atom.residue_id = bead_idx; // XXX
                }
                cg_bead->bead_id() = bead_idx;
                ++bead_idx;
            }
        }
    }
    {// output coarse-grained structure to .cgs file
        const std::string cgsfile(file_name + ".cgs");
        std::ofstream ofs(cgsfile);
        if(!ofs.good())
        {
            std::cerr << "cannot write Coarse-Grained Structure file: "
                      << cgsfile << std::endl;
            return 1;
        }
        for(auto& cg_chain : cg_chains)
        {
            for(auto& cg_bead : cg_chain)
            {
                cg_bead->print(ofs);
                ofs << std::endl;
            }
        }
    }

    // generate forcefield parameters
    const auto forcefields =
        toml::get<toml::Array<toml::Table>>(input.at("forcefields"));

//     std::vector<jarngreipr::ForceFieldParameterTable<traits>>
//         ffpts(forcefields.size());
//
//     for(auto iter = forcefields.cbegin(); iter != forcefields.cend(); ++iter)
//     {
//         std::size_t index;
//         try
//         {
//             index = toml::get<toml::Integer>(iter->at("index"));
//             if(ffpts.size() <= index)
//             {
//                 std::cerr << "[Error] forcefields.index exceeds number of table"
//                           << std::endl;
//                 return 1;
//             }
//         }
//         catch(std::out_of_range&)
//         {
//             index = 0;
//             std::cerr << "[WARNING] no 'index' field in [[forcefields]] table."
//                       << "set 0." << std::endl;
//         }
//
//         const auto local = toml::get<toml::Table>(iter->at("local"));
//         for(auto chains = local.cbegin(); chains != local.cend(); ++chains)
//         {
//             for(auto chainID : jarngreipr::unpack_chainID(chain.first))
//             {
//                 const auto chain = std::find_if(
//                         cg_chains.cbegin(), cg_chains.cend(), []()
//                         );
//                 apply_forcefield();
//             }
//         }
//
//
//
//
//         const auto global = toml::get<toml::Table>(iter->at("global"));
//
//     }
//
//     // simulator setup
//     const auto simulators =
//         toml::get<toml::Array<toml::Table>>(input.at("simulators"));
//
//     for(auto iter = simulators.cbegin(); iter != simulators.cend(); ++iter)
//     {
//         try
//         {
//         const std::size_t index = toml::get<toml::Integer>(iter->at("index"));
//         if(ffpts.size() >= index)
//             ffpts.resize(index+1);
//         }
//         catch(std::out_of_range&)
//         {
//         std::cerr << "[WARNING] no 'index' field in [[simulators]] table"
//                   << std::endl;
//         }
//     }
    return 0;
}
