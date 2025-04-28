#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>
#include <optional>

// Thread-safe blocking queue
// T: type of stored items
// push(): add an item and notify one waiter
// pop(): block until an item is available or queue is closed, then return optional<T>
// close(): wake all waiters and prevent further push/pop
// empty(): check if queue is empty

template <typename T>
class BlockingQueue {
private:
    std::queue<T>              queue_;       // underlying FIFO container
    mutable std::mutex         mutex_;       // guards queue_
    std::condition_variable    condVar_;     // signals push/close events
    bool                        closed_ = false;

public:
    // Add an item and wake one waiting thread
    void push(const T& item) {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            queue_.push(item);
        }
        condVar_.notify_one();
    }

    // Block until item available or closed. Returns nullopt if closed and empty.
    std::optional<T> pop() {
        std::unique_lock<std::mutex> lock(mutex_);
        condVar_.wait(lock, [this] { return closed_ || !queue_.empty(); });
        if (queue_.empty())
            return std::nullopt;
        T value = std::move(queue_.front());
        queue_.pop();
        return value;
    }

    // Prevent further push/pop and wake all waiters
    void close() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            closed_ = true;
        }
        condVar_.notify_all();
    }

    // Check if queue has no items
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }
};
