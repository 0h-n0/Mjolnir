#ifndef JARNGREIPR_MODEL_INTERACTION
#define JARNGREIPR_MODEL_INTERACTION
#include "output_parameters.hpp"
#include <tuple>
#include <utility>
#include <string>
#include <array>

namespace jarngreipr
{

template<typename traitsT, std::size_t N, typename ... paramTs>
struct LocalInteractionParameter
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::array<std::size_t, N> indices_type;
    typedef std::tuple<std::pair<std::string, paramTs>...> params_type;

  public:
    LocalInteractionParameter()  = default;
    virtual ~LocalInteractionParameter() = default;

    std::size_t& index_at(std::size_t i)       {return indices_.at(i);}
    std::size_t  index_at(std::size_t i) const {return indices_.at(i);}

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

    indices_type const& indices()    const {return indices_;}
    params_type  const& parameters() const {return paramters_;}

  private:
    indices_type indices_;
    params_type parameters_;
};

template<typename charT, typename char_traits,
         typename traitsT, std::size_t N, typename ... Ts>
std::basic_ostream<charT, char_traits>&
operator<<(std::basic_ostream<charT, char_traits>& os,
           const LocalInteractionParameter<traitsT, N, Ts...>& interaction)
{
    os << "{indices=[";
    for(auto iter : interaction.indices())
        os << item << ',';
    os << "],";
    detail::output_params_impl<charT, char_traits, traitsT,
        sizeof...(Ts), Ts...>::invoke(os, interaction.parameters());
    os << '}';
    return os;
}

}// jarngreipr
#endif //JARNGREIPR_MODEL_INTERACTION
