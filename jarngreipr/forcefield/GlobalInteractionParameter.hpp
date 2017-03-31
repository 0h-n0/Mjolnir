#ifndef JARNGREIPR_MODEL_GLOBAL_INTERACTION_PARAMETER
#define JARNGREIPR_MODEL_GLOBAL_INTERACTION_PARAMETER
#include "output_parameters.hpp"
#include <tuple>
#include <utility>
#include <string>
#include <array>
#include <ostream>

namespace jarngreipr
{

template<typename traitsT, typename ... paramTs>
struct GlobalInteractionParameter
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::tuple<std::pair<std::string, paramTs>...> params_type;

  public:
    GlobalInteractionParameter()  = default;
    virtual ~GlobalInteractionParameter() = default;

    std::size_t& index()       {return index;}
    std::size_t  index() const {return index;}

    template<std::size_t I>
    typename std::tuple_element<I, params_type>::type&
    parameter_at()
    {
        return std::get<I>(parameters_);
    }

    template<std::size_t I>
    typename std::tuple_element<I, params_type>::type const&
    parameter_at() const
    {
        return std::get<I>(parameters_);
    }

    params_type const& parameters() const {return paramters_;}

  private:
    std::size_t index_;
    params_type parameters_;
};

template<typename charT, typename char_traits,
         typename traitsT, std::size_t N, typename ... Ts>
std::basic_ostream<charT, char_traits>&
operator<<(std::basic_ostream<charT, char_traits>& os,
           const GlobalInteractionParameter<traitsT, N, Ts...>& interaction)
{
    os << "{index=" << interaction.index() << ",";
    detail::output_params_impl<charT, char_traits, traitsT,
        sizeof...(Ts), Ts...>::invoke(os, interaction.parameters());
    os << '}';
    return os;
}

}// jarngreipr
#endif //JARNGREIPR_MODEL_GLOBAL_INTERACTION_PARAMETER
