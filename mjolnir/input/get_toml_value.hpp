#ifndef MJOLNIR_INPUT_GET_TOML_VALUE
#define MJOLNIR_INPUT_GET_TOML_VALUE
#include <extlib/toml/toml.hpp>
#include <string>
#include <stdexcept>

namespace mjolnir
{
namespace detail
{

inline toml::value const&
value_at(const toml::Table& tab, const std::string& key,
         const std::string& tablename = "<couldn't provided>")
{
    try
    {
        return tab.at(key);
    }
    catch(const std::out_of_range& oor)
    {
        throw std::out_of_range("mjolnir: while reading toml file: key " + key +
                std::string("in table ") + tablename + std::string(" missing"));
    }
}

} // detail
} // mjolnir
#endif// MJOLNIR_INPUT_GET_TOML_VALUE