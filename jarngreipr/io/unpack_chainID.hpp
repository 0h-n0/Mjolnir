#ifndef JARNGREIPR_IO_UNPACK_CHAIN_ID
#define JARNGREIPR_IO_UNPACK_CHAIN_ID
#include <megingjord/util/split.hpp>
#include "chainID_iterator.hpp"

namespace jarngreipr
{

/*! @brief unpack string "A-D" to array of {"A", "B", "C", "D"}  *
 * pairof{"A-D", value} ->                                       *
 *     arrayof{pairof{"A", value}, pairof{"B", value}, ...}      */
template<typename T>
inline std::vector<std::pair<std::string, T>>
unpack_chainID(const std::pair<std::string, T>& value)
{
    const auto splitted = megingjord::split(value.first, '-');
    if(splitted.size() > 2)
        throw std::invalid_argument("too much chains");

    std::vector<std::pair<std::string, T>> retval;
    for(chainID_iterator iter(splitted.front()), last(splitted.back());
            iter != std::next(last); ++iter)
    {
        retval.emplace_back(*iter, value.second);
    }
    return retval;
}

}// jarngreipr
#endif // JARNGREIPR_IO_UNPACK_CHAIN_ID
