#ifndef JARNGREIPR_OUTPUT_PARAMETERS
#define JARNGREIPR_OUTPUT_PARAMETERS
#include <ostream>
#include <tuple>

namespace jarngreipr
{
namespace detail
{

template<typename charT, typename char_traits, std::size_t I, typename ... Ts>
struct output_params_impl
{
    static std::basic_ostream<charT, char_traits>&
    invoke(std::basic_ostream<charT, char_traits>& os,
           const std::tuple<std::pair<std::string, Ts>...>& paras)
    {
        os << std::get<sizeof...(Ts) - I>(paras).first << '='
           << std::get<sizeof...(Ts) - I>(paras).second << ',';
        return output_params_impl<charT, char_traits, I-1, Ts...>::invoke(os, paras);
    }
};

template<typename charT, typename char_traits, typename ... Ts>
struct output_params_impl<charT, char_traits, traitsT, 1, Ts...>
{
    static std::basic_ostream<charT, char_traits>&
    invoke(std::basic_ostream<charT, char_traits>& os,
           const std::tuple<std::pair<std::string, Ts>...>& paras)
    {
        os << std::get<sizeof...(Ts)-1>(paras).first << '='
           << std::get<sizeof...(Ts)-1>(paras).second;
        return os;
    }
};

}//detail
}// jarngreipr
#endif //JARNGREIPR_OUTPUT_PARAMETERS
