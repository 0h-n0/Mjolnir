#ifndef MJOLNIR_UTIL_TIMESTAMP
#define MJOLNIR_UTIL_TIMESTAMP
#include <chrono>

namespace mjolnir
{

template<typename T>
class timestamp
{
  public:
    typedef T value_type;
    typedef std::chrono::system_clock::time_point time_point_type;

  public:

    timestamp() = default;
    ~timestamp() = default;

    void start(){start_ = std::chrono::system_clock::now();}
    void stop() {stop_  = std::chrono::system_clock::now();}

    value_type const& value() const {return value_;}
    value_type&       value()       {return value_;}

    template<typename durationT>
    durationT start_at() const
    {
        return std::chrono::duration_cast<durationT>(start_.time_since_epoch());
    }

    template<typename durationT>
    durationT stop_at() const
    {
        return std::chrono::duration_cast<durationT>(stop_.time_since_epoch());
    }

    template<typename durationT>
    durationT elapsed() const
    {
        return std::chrono::duration_cast<durationT>(stop_ - start_);
    }


    time_point_type const& get_start() const {return start_;}
    time_point_type const& get_stop()  const {return stop_;}
 
  private:

    value_type      value_;
    time_point_type start_;
    time_point_type stop_;
};

template<>
class timestamp<void>
{
  public:
    typedef std::chrono::system_clock::time_point time_point_type;

  public:

    timestamp() = default;
    ~timestamp() = default;

    void start(){start_ = std::chrono::system_clock::now();}
    void stop() {stop_  = std::chrono::system_clock::now();}

    template<typename durationT>
    durationT start_at() const
    {
        return std::chrono::duration_cast<durationT>(start_.time_since_epoch());
    }

    template<typename durationT>
    durationT stop_at() const
    {
        return std::chrono::duration_cast<durationT>(stop_.time_since_epoch());
    }

    template<typename durationT>
    durationT elapsed() const
    {
        return std::chrono::duration_cast<durationT>(stop_ - start_);
    }


    time_point_type const& get_start() const {return start_;}
    time_point_type const& get_stop()  const {return stop_;}
 
  private:

    time_point_type start_;
    time_point_type stop_;
};


} // mjolnir
#endif/* MJOLNIR_UTIL_TIMESTAMP */
