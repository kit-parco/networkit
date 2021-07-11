#ifndef NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_
#define NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_

#include <atomic>

namespace Aux {

class Spinlock {
public:
    void lock() {
        while (spinner.test_and_set(std::memory_order_acquire)) {
            /* spin */
        }
    }
    void unlock() { spinner.clear(std::memory_order_release); }

private:
    std::atomic_flag spinner = ATOMIC_FLAG_INIT;
};

} // namespace Aux

#endif // NETWORKIT_AUXILIARY_SPIN_LOCK_HPP_
