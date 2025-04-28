#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>

// A thread-safe blocking queue matching push/pop/finish API
template<typename T>
class BlockingQueue {
private:
    std::queue<T>           queue_;
    mutable std::mutex      mutex_;
    std::condition_variable cv_;
    bool                    finished_ = false;

public:
    // Enqueue a new item
    void push(const T& item) {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            queue_.push(item);
        }
        cv_.notify_one();
    }

    // Block until an item is available or finish() called; returns false if done and empty
    bool pop(T& out) {
        std::unique_lock<std::mutex> lk(mutex_);
        cv_.wait(lk, [this]{ return finished_ || !queue_.empty(); });
        if (queue_.empty()) return false;
        out = std::move(queue_.front());
        queue_.pop();
        return true;
    }

    // Signal that no more items will be added, unblocking all waiting threads
    void finish() {
        {
            std::lock_guard<std::mutex> lk(mutex_);
            finished_ = true;
        }
        cv_.notify_all();
    }

    // Check whether the queue is empty (non-blocking)
    bool empty() const {
        std::lock_guard<std::mutex> lk(mutex_);
        return queue_.empty();
    }
};
