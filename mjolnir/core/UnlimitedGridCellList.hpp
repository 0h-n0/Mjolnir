#ifndef MJOLNIR_UNLIMITED_GRID_CELL_LIST
#define MJOLNIR_UNLIMITED_GRID_CELL_LIST
#include "BoundaryCondition.hpp"
#include "SpatialPartition.hpp"
#include <mjolnir/util/logger.hpp>
#include <functional>
#include <algorithm>
#include <limits>
#include <array>
#include <cmath>

namespace mjolnir
{

template<typename traitsT, typename boundaryT = UnlimitedBoundary<traitsT>>
class UnlimitedGridCellList : public SpatialPartition<traitsT>
{
  public:

    typedef boundaryT boundary_type;
    typedef SpatialPartition<traitsT> base_type;
    typedef typename base_type::time_type time_type;
    typedef typename base_type::real_type real_type;
    typedef typename base_type::coordinate_type coordinate_type;
    typedef typename base_type::coordinate_type position_type;
    typedef typename base_type::particle_container_type particle_container_type;
    typedef typename base_type::index_type index_type;
    typedef typename base_type::index_list index_list;
    constexpr static std::size_t dim_size   = 16;
    constexpr static int         dim        = dim_size;
    constexpr static std::size_t total_size = dim_size * dim_size * dim_size;
    constexpr static real_type mesh_epsilon = 1e-6;
    typedef std::array<int, 3>               cell_index_type;
    typedef std::array<std::size_t, 26>      neighbor_cell_idx;
    typedef std::pair<index_list, neighbor_cell_idx> unit_cell_type;
    typedef std::vector<index_list>          verlet_list_type;
    typedef std::vector<index_list>          except_list_type;
    typedef std::array<unit_cell_type, total_size> cell_list_type;

  public:

    UnlimitedGridCellList() = default;
    ~UnlimitedGridCellList() override = default;

    UnlimitedGridCellList(const real_type cutoff, const real_type mergin)
        : dt_(0.), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {
        initialize();
    }
    UnlimitedGridCellList(const real_type cutoff, const real_type mergin,
                          const time_type dt)
        : dt_(dt), cutoff_(cutoff), mergin_(mergin), current_mergin_(-1.),
          inv_cell_size_(1. / (cutoff + mergin + mesh_epsilon))
    {
        initialize();
    }

    bool valid() const noexcept override
    {
        return current_mergin_ >= 0. || dt_ == 0.;
    }

    void initialize();
    void make  (const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon) override;
    void update(const particle_container_type& pcon,
                const time_type dt) override;

    real_type cutoff() const {return this->cutoff_;}
    real_type mergin() const {return this->mergin_;}

    void set_cutoff(const real_type c);
    void set_mergin(const real_type m);

  private:


    std::size_t index(const position_type& pos)   const;
    std::size_t index(const cell_index_type& pos) const;
    cell_index_type add(const int x, const int y, const int z,
                        const cell_index_type&) const;

  private:

    real_type dt_;
    real_type cutoff_;
    real_type mergin_;
    real_type current_mergin_;
    real_type inv_cell_size_;

    cell_list_type cell_list_;

    static Logger& logger_;
};

template<typename traitsT, typename boundaryT>
Logger& UnlimitedGridCellList<traitsT, boundaryT>::logger_ =
        LoggerManager<char>::get_logger("UnlimitedGridCellList");

template<typename traitsT, typename boundaryT>
inline void
UnlimitedGridCellList<traitsT, boundaryT>::set_cutoff(const real_type c)
{
    this->cutoff_ = c;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
inline void
UnlimitedGridCellList<traitsT, boundaryT>::set_mergin(const real_type m)
{
    this->mergin_ = m;
    this->inv_cell_size_ = 1. / (cutoff_ + mergin_ + mesh_epsilon);
    return;
}

template<typename traitsT, typename boundaryT>
void UnlimitedGridCellList<traitsT, boundaryT>::make(
        const particle_container_type& pcon)
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::make CALLED");

    for(auto iter = this->list_.begin(); iter != this->list_.end(); ++iter)
        iter->clear();
    this->list_.resize(pcon.size());

    for(auto iter = cell_list_.begin(); iter != cell_list_.end(); ++iter)
        iter->first.clear(); // DON'T clear iter->second

    MJOLNIR_LOG_DEBUG("cell_list and verlet_list are cleared");

    std::size_t idx = 0;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
    {
        MJOLNIR_LOG_DEBUG("set", idx, "-th particle at", index(iter->position));
        cell_list_[index(iter->position)].first.push_back(idx);
        ++idx;
    }
    MJOLNIR_LOG_DEBUG("cell list is updated");

    const real_type r_c = cutoff_ + mergin_;
    const real_type r_c2 = r_c * r_c;

    MJOLNIR_LOG_DEBUG("except list size", this->except_.size());

    MJOLNIR_LOG_DEBUG("exception list is not empty.");
    MJOLNIR_LOG_DEBUG("lookup particles and also except list.");
    for(std::size_t i=0; i<pcon.size(); ++i)
    {
        const coordinate_type ri   = pcon[i].position;
        const unit_cell_type& cell = cell_list_.at(index(ri));

        MJOLNIR_LOG_DEBUG("particle position", pcon[i].position);
        MJOLNIR_LOG_DEBUG("cell index", index(ri));
        MJOLNIR_LOG_DEBUG("making verlet list for index", i);
        MJOLNIR_LOG_DEBUG("except list for ", i, "-th value has",
                          this->except_.at(i).size(), "particles");

        const auto cbeg = this->except_.at(i).cbegin();
        const auto cend = this->except_.at(i).cend();

        for(std::size_t j=0; j < cell.first.size(); ++j)
        {
            MJOLNIR_LOG_DEBUG("looking", j, "-th particle in the cell.",
                              "its index is", cell.first.at(j));

            const std::size_t k = cell.first.at(j);
            if(k <= i || std::find(cbeg, cend, k) != cend)
                continue;

            if(length_sq(boundary_type::adjust_direction(
                            pcon.at(k).position - ri)) < r_c2)
            {
                MJOLNIR_LOG_DEBUG("add index", k, "to verlet list", i);
                this->list_.at(i).push_back(k);
            }
        }
        MJOLNIR_LOG_DEBUG("end looking in the same cell");

        // see neighbor cells
        for(auto iter = cell.second.cbegin(); iter != cell.second.cend(); ++iter)
        {
            MJOLNIR_LOG_DEBUG("see neighboring cell at", std::distance(
                              cell.second.cbegin(), iter));
            MJOLNIR_LOG_DEBUG("neighboring cell index", *iter);

            const unit_cell_type& neighbor = cell_list_[*iter];
            for(std::size_t j=0; j < neighbor.first.size(); ++j)
            {
                MJOLNIR_LOG_DEBUG("looking", j, "-th particle in the cell.",
                                  "its index is", neighbor.first.at(j));
                const std::size_t k = neighbor.first.at(j);
                if(k <= i || std::find(cbeg, cend, k) != cend)
                    continue;

                if(length_sq(boundary_type::adjust_direction(
                                pcon.at(k).position - ri)) < r_c2)
                {
                    MJOLNIR_LOG_DEBUG("add index", k, "to verlet list", i);
                    this->list_.at(i).push_back(k);
                }
            }
        }
    }

    for(auto iter = this->list_.begin(); iter != this->list_.end(); ++iter)
    {
        std::sort(iter->begin(), iter->end());
    }

    this->current_mergin_ = mergin_;

    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList::make() RETURNED");
    return ;
}

template<typename traitsT, typename boundaryT>
void UnlimitedGridCellList<traitsT, boundaryT>::update(
        const particle_container_type& pcon)
{
    real_type max_speed = 0.;
    for(auto iter = pcon.cbegin(); iter != pcon.cend(); ++iter)
        max_speed = std::max(max_speed, length_sq(iter->velocity));

    this->current_mergin_ -= std::sqrt(max_speed) * dt_ * 2.;
    if(this->current_mergin_ < 0.)
        this->make(pcon);

    return ;
}

template<typename traitsT, typename boundaryT>
inline void UnlimitedGridCellList<traitsT, boundaryT>::update(
        const particle_container_type& pcon, const time_type dt)
{
    this->dt_ = dt;
    this->update(pcon);
    return ;
}

template<typename traitsT, typename boundaryT>
inline std::size_t
UnlimitedGridCellList<traitsT, boundaryT>::index(const position_type& pos) const
{
    const int x = static_cast<int>(std::floor(pos[0]*inv_cell_size_)) % dim;
    const int y = static_cast<int>(std::floor(pos[1]*inv_cell_size_)) % dim;
    const int z = static_cast<int>(std::floor(pos[2]*inv_cell_size_)) % dim;
    return index(cell_index_type{{(x < 0) ? x + dim : x,
                                  (y < 0) ? y + dim : y,
                                  (z < 0) ? z + dim : z}});
}

template<typename traitsT, typename boundaryT>
inline std::size_t
UnlimitedGridCellList<traitsT, boundaryT>::index(const cell_index_type& idx) const
{
    return idx[0] + dim_size * idx[1] + dim_size * dim_size * idx[2];
}


template<typename traitsT, typename boundaryT>
inline typename UnlimitedGridCellList<traitsT, boundaryT>::cell_index_type
UnlimitedGridCellList<traitsT, boundaryT>::add(
        const int x, const int y, const int z, const cell_index_type& idx) const
{
    int ret_x = (idx[0] + x);
    int ret_y = (idx[1] + y);
    int ret_z = (idx[2] + z);
    if(ret_x < 0) ret_x += dim; else if(ret_x >= dim) ret_x -= dim;
    if(ret_y < 0) ret_y += dim; else if(ret_y >= dim) ret_y -= dim;
    if(ret_z < 0) ret_z += dim; else if(ret_z >= dim) ret_z -= dim;
    return cell_index_type{{ret_x, ret_y, ret_z}};
}

template<typename traitsT, typename boundaryT>
void UnlimitedGridCellList<traitsT, boundaryT>::initialize()
{
    MJOLNIR_LOG_DEBUG("UnlimitedGridCellList<traitsT>::initialize CALLED");
    for(auto iter = cell_list_.begin(); iter != cell_list_.end(); ++iter)
    {
        MJOLNIR_LOG_DEBUG("filled by null", std::distance(cell_list_.begin(), iter));
        iter->first = index_list{};
        iter->second.fill(std::numeric_limits<std::size_t>::max());
    }

    for(int x = 0; x < dim; ++x)
    for(int y = 0; y < dim; ++y)
    for(int z = 0; z < dim; ++z)
    {
        const cell_index_type idx{{x, y, z}};
        auto& cell = this->cell_list_[index(idx)];

        MJOLNIR_LOG_DEBUG("cell", x, y, z, "index", index(idx));
        cell.second[ 0] = index(add( 1,  0,  0, idx));
        cell.second[ 1] = index(add( 0,  1,  0, idx));
        cell.second[ 2] = index(add( 0,  0,  1, idx));
        cell.second[ 3] = index(add(-1,  0,  0, idx));
        cell.second[ 4] = index(add( 0, -1,  0, idx));
        cell.second[ 5] = index(add( 0,  0, -1, idx));
        cell.second[ 6] = index(add( 1,  1,  0, idx));
        cell.second[ 7] = index(add( 0,  1,  1, idx));
        cell.second[ 8] = index(add( 1,  0,  1, idx));
        cell.second[ 9] = index(add(-1, -1,  0, idx));
        cell.second[10] = index(add( 0, -1, -1, idx));
        cell.second[11] = index(add(-1,  0, -1, idx));
        cell.second[12] = index(add( 1, -1,  0, idx));
        cell.second[13] = index(add( 0,  1, -1, idx));
        cell.second[14] = index(add(-1,  0,  1, idx));
        cell.second[15] = index(add(-1,  1,  0, idx));
        cell.second[16] = index(add( 0, -1,  1, idx));
        cell.second[17] = index(add( 1,  0, -1, idx));
        cell.second[18] = index(add(-1,  1,  1, idx));
        cell.second[19] = index(add( 1, -1,  1, idx));
        cell.second[20] = index(add( 1,  1, -1, idx));
        cell.second[21] = index(add(-1, -1,  1, idx));
        cell.second[22] = index(add( 1, -1, -1, idx));
        cell.second[23] = index(add(-1,  1, -1, idx));
        cell.second[24] = index(add( 1,  1,  1, idx));
        cell.second[25] = index(add(-1, -1, -1, idx));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  0,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1,  0, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 0, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  0, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1, -1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1,  1, -1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add( 1,  1,  1, idx)));
        MJOLNIR_LOG_DEBUG("neighbors", index(add(-1, -1, -1, idx)));
    }
    return;
}



} // mjolnir
#endif/* MJOLNIR_UNLIMITED_GRID_CELL_LIST */
