#ifndef JARNGREIPR_MODEL_INTERACTION
#define JARNGREIPR_MODEL_INTERACTION
#include <tuple>
#include <utility>
#include <string>
#include <ostream>

namespace jarngreipr
{

template<typename traitsT, typename ... paramTs>
struct Interaction
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef std::vector<std::size_t> indices_type;
    typedef std::tuple<std::pair<std::string, paramTs>...> params_type;

  public:
    InteractionBase()  = default;
    virtual ~InteractionBase() = default;

    void add_index(std::size_t i){indices_.emplace_back(i);}

    template<std::size_t I>
    void add_parameter(const std::tuple_element<I, params_type>::type& para)
    {
        std::get<I>(parameters_) = para;
    }

    indices_type const& indices()    const {return indices_;}
    params_type  const& parameters() const {return paramters_;}

  private:
    indices_type indices_;
    params_type parameters_;
};

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

template<typename charT, typename char_traits, typename traitsT, typename ... Ts>
std::basic_ostream<charT, char_traits>&
operator<<(std::basic_ostream<charT, char_traits>& os,
           const Interaction<traitsT, Ts...>& interaction)
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
