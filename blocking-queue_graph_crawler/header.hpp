#pragma once
#include <queue>
#include <mutex>
#include <condition_variable>
#include <optional>

template<typename T>
class BlockingQueue {
private:
    std::queue<T>           queue_;
    mutable std::mutex      mutex_;
    std::condition_variable cond_var_;
    bool                    finished_ = false;

public:
    // enqueue a new item
    void enqueue(const T& item) {
        std::lock_guard<std::mutex> lock(mutex_);
        queue_.push(item);
        cond_var_.notify_one();
    }

    // try to dequeue; returns false if the queue is empty and closed
    bool try_dequeue(T& out) {
        std::unique_lock<std::mutex> lock(mutex_);
        cond_var_.wait(lock, [this]{ return finished_ || !queue_.empty(); });
        if (queue_.empty()) return false;
        out = queue_.front();
        queue_.pop();
        return true;
    }

    // signal that no more items will be added
    void close() {
        {
            std::lock_guard<std::mutex> lock(mutex_);
            finished_ = true;
        }
        cond_var_.notify_all();
    }

    // check if empty (at the moment of the call)
    bool empty() const {
        std::lock_guard<std::mutex> lock(mutex_);
        return queue_.empty();
    }
};
