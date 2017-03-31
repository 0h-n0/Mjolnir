#ifndef JARNGREIPR_IO_CHAIN_ID_ITERATOR
#define JARNGREIPR_IO_CHAIN_ID_ITERATOR
#include <string>
#include <iterator>
#include <stdexcept>

namespace jarngreipr
{

struct chainID_iterator
{
    // represent only range [A, ZZ]
    typedef std::size_t        difference_type;
    typedef std::string        value_type;
    typedef std::string const* pointer;
    typedef std::string const& reference;
    typedef std::bidirectional_iterator_tag iterator_category;

    chainID_iterator() : id("A"){}
    chainID_iterator(std::string const& str) : id(str){}
    ~chainID_iterator() = default;
    chainID_iterator(chainID_iterator const&) = default;
    chainID_iterator(chainID_iterator &&)     = default;
    chainID_iterator& operator=(chainID_iterator const&) = default;
    chainID_iterator& operator=(chainID_iterator &&)     = default;

    chainID_iterator& operator++()
    {
        if(id.size() == 1)
        {
            if(id.front() == 'Z')
                id = "AA";
            else
                id.front() += 1;
        }
        else if(id.back() == 'Z')
        {
            if(id.front() == 'Z')
                throw std::out_of_range("too many chains");
            else
                id.front() += 1;
            id.back() = 'A';
        }
        else
        {
            id.back() += 1;
        }
        return *this;
    }

    chainID_iterator& operator--()
    {
        if(id.size() == 1)
        {
            if(id.front() == 'A')
                throw std::out_of_range("negative chains");
            else
                id.front() -= 1;
        }
        else if(id.back() == 'A')
        {
            if(id.front() == 'A')
                id = "Z";
            else
                id.front() -= 1;
            id.back() = 'Z';
        }
        else
        {
            id.back() -= 1;
        }
        return *this;
    }

    chainID_iterator  operator++(int)
    {
        chainID_iterator retval = *this;
        ++(*this);
        return retval;
    }

    chainID_iterator  operator--(int)
    {
        chainID_iterator retval = *this;
        --(*this);
        return retval;
    }

    reference operator*()  const {return id;}
    pointer   operator->() const {return &id;}

    std::string id;
};

inline bool operator==(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    return lhs.id == rhs.id;
}

inline bool operator!=(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    return !(lhs == rhs);
}

inline bool operator<(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    if(lhs.id.size() != rhs.id.size()) return lhs.id.size() < rhs.id.size();
    else if(lhs.id.front() != rhs.id.front()) return lhs.id.front() < rhs.id.front();
    else return lhs.id.back() < rhs.id.back();
}

inline bool operator<=(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    return (lhs < rhs || lhs == rhs);
}

inline bool operator>(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    return !(lhs <= rhs);
}

inline bool operator>=(chainID_iterator const& lhs, chainID_iterator const& rhs)
{
    return !(lhs < rhs);
}

}// jarngreipr
#endif //JARNGREIPR_IO_CHAIN_ID_ITERATOR
