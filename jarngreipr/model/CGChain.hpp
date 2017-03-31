#ifndef JARNGREIPR_CG_CHAIN
#define JARNGREIPR_CG_CHAIN
#include "Bead.hpp"
#include <memory>

namespace jarngreipr
{

template<typename traitsT>
class CGChain
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef Bead<traits_type> bead_type;
    typedef typename bead_type::atom_type atom_type;
    typedef std::unique_ptr<bead_type> bead_ptr;
    typedef std::vector<bead_ptr> container_type;
    typedef typename container_type::iterator       iterator;
    typedef typename container_type::const_iterator const_iterator;

  public:

    CGChain()  = default;
    ~CGChain() = default;
    CGChain(CGChain&&)             = default;
    CGChain& operator=(CGChain&&)  = default;
    CGChain(const CGChain&)             = delete;
    CGChain& operator=(const CGChain&)  = delete;

    void push_back(bead_ptr&& bead)
    {
        beads_.push_back(std::forward<bead_ptr>(bead));
    }

    template<typename ...Ts>
    void emplace_back(Ts&& ... args)
    {
        beads_.emplace_back(std::forward<Ts>(args)...);
    }

    std::string&       chainID()       {return id_;}
    std::string const& chainID() const {return id_;}

    bead_ptr&       operator[](std::size_t i)       {return beads_[i];}
    bead_ptr const& operator[](std::size_t i) const {return beads_[i];}
    bead_ptr&       at(std::size_t i)       {return beads_.at(i);}
    bead_ptr const& at(std::size_t i) const {return beads_.at(i);}

    iterator begin() {return beads_.begin();}
    iterator end()   {return beads_.end();}
    const_iterator begin() const {return beads_.cbegin();}
    const_iterator end()   const {return beads_.cend();}
    const_iterator cbegin() const {return beads_.cbegin();}
    const_iterator cend()   const {return beads_.cend();}

  private:

    std::string    id_;
    container_type beads_;
};

template<typename T>
inline std::ostream& operator<<(std::ostream& os, const CGChain<T>& chain)
{
    for(auto bead : chain)
    {
        bead->print(os);
        os << std::endl;
    }
    return os;
}

}//jarngreipr
#endif// JARNGREIPR_CG_CHAIN
