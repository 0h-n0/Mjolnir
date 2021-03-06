#ifndef MJOLNIR_READ_FORCEFIELD
#define MJOLNIR_READ_FORCEFIELD
#include <extlib/toml/toml.hpp>
#include <mjolnir/core/ForceField.hpp>
#include <mjolnir/util/get_toml_value.hpp>
#include <mjolnir/input/read_interaction.hpp>

namespace mjolnir
{

template<typename traitsT>
LocalForceField<traitsT>
read_local_forcefield(std::vector<toml::Table> interactions)
{
    LocalForceField<traitsT> lff;
    for(const auto& interaction : interactions)
    {
        lff.emplace(read_local_interaction<traitsT>(interaction));
    }
    return lff;
}

template<typename traitsT>
GlobalForceField<traitsT>
read_global_forcefield(std::vector<toml::Table> interactions)
{
    GlobalForceField<traitsT> gff;
    for(const auto& interaction : interactions)
    {
        gff.emplace(read_global_interaction<traitsT>(interaction));
    }
    return gff;
}

template<typename traitsT>
ExternalForceField<traitsT>
read_external_forcefield(std::vector<toml::Table> interactions)
{
    ExternalForceField<traitsT> eff;
    for(const auto& interaction : interactions)
    {
        eff.emplace(read_external_interaction<traitsT>(interaction));
    }
    return eff;
}

template<typename traitsT>
ForceField<traitsT>
read_forcefield(const toml::Table& data, std::size_t N)
{
    const auto& ffs = toml_value_at(data, "forcefields", "<root>"
            ).cast<toml::value_t::Array>();

    if(ffs.size() <= N)
    {
        throw std::out_of_range("no enough forcefields: " + std::to_string(N));
    }
    const auto& ff = ffs.at(N).cast<toml::value_t::Table>();

    std::vector<toml::Table> fflocal, ffglobal, ffexternal;

    if(ff.count("local") == 1)
    {
        fflocal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "local", "[forcefields]"));
    }
    if(ff.count("global") == 1)
    {
        ffglobal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "global", "[forcefields]"));
    }
    if(ff.count("external") == 1)
    {
        ffexternal = toml::get<std::vector<toml::Table>>(
                toml_value_at(ff, "external", "[forcefields]"));
    }
    return ForceField<traitsT>(
            read_local_forcefield<traitsT>(std::move(fflocal)),
            read_global_forcefield<traitsT>(std::move(ffglobal)),
            read_external_forcefield<traitsT>(std::move(ffexternal)));
}


} // mjolnir
#endif// MJOLNIR_READ_FORCEFIELD
