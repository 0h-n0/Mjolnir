#ifndef JARNGREIPR_MODEL
#define JARNGREIPR_MODEL
#include "Bead.hpp"
#include <utility>
#include <memory>
#include <array>

namespace jarngreipr
{

template<typename traitsT>
class Model
{
  public:
    typedef traitsT traits_type;
    typedef typename traits_type::real_type real_type;
    typedef typename traits_type::coordinate_type coordinate_type;
    typedef PDBChain<traits_type> chain_type;
    typedef Bead<traits_type> bead_type;
    typedef std::unique_ptr<bead_type> bead_ptr;
    typedef std::vector<bead_ptr> bead_container_type;
    typedef typename bead_container_type::iterator       iterator;
    typedef typename bead_container_type::const_iterator const_iterator;

  public:

    Model() = default;
    ~Model() = default;
};

}//jarngreipr
#endif /* JARNGREIPR_MODEL */
