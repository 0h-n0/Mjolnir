#ifndef MJOLNIR_SPIN_LOCK
#define MJOLNIR_SPIN_LOCK
#include <atomic>

namespace mjolnir
{

class spinlock
{
  public:
    spinlock(): state_(ATOMIC_FLAG_INIT){}
    ~spinlock() = default;

    void lock()
    {
        while(state_.test_and_set(std::memory_order_acquire)){/* spin! */;}
        return;
    }

    void unlock()
    {
        state_.clear(std::memory_order_release);
        return;
    }

  private:
    std::atomic_flag state_;
};

} // mjolnir
#endif /* MJOLNIR_SPIN_LOCK */
