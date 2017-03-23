#ifndef MEGINGJORD_UTIL_STRING_HPP
#define MEGINGJORD_UTIL_STRING_HPP
#include <string>
#include <vector>

namespace megingjord
{

template<typename charT, typename ctrait, typename alloc>
std::vector<std::basic_string<charT, ctrait, alloc>>
split(const std::basic_string<charT, ctrait, alloc>& str, charT delim)
{
    std::vector<std::basic_string<charT, ctrait, alloc>> retval;
    std::size_t offset = 0;
    while(offset < str.size())
    {
        const std::size_t next = str.find_first_of(delim, offset);
        if(next == std::basic_string<charT, ctrait, alloc>::npos)
        {
            retval.emplace_back(str.cbegin()+offset, str.cend());
            break;
        }
        retval.emplace_back(str.cbegin()+offset, str.cbegin() + next);
        offset = next+1;
    }
    return retval;
}


} // megingjord
#endif// MEGINGJORD_UTIL_STRING_HPP
