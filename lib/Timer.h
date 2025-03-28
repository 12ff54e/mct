#ifndef ZQ_TIMER
#define ZQ_TIMER

#include <chrono>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

class Timer {
   public:
    static Timer& get_timer();
    void start(std::string func_name);
    void pause_last_and_start_next(std::string func_name);
    void pause();
    void pause(std::string func_name);
    void reset();
    void print();

   private:
    Timer() = default;
    Timer(const Timer&) = delete;
    Timer(Timer&&) = delete;
    Timer& operator=(const Timer&) = delete;
    Timer& operator=(Timer&&) = delete;

    std::vector<std::string> entries;
    std::unordered_map<
        std::string,
        std::pair<std::chrono::high_resolution_clock::duration,
                  std::chrono::high_resolution_clock::time_point>>
        time_consuming;

    std::string current_func_name;
};

#ifdef ZQ_TIMER_IMPLEMENTATION

#include <iomanip>
#include <iostream>

using namespace std::chrono;

void Timer::start(std::string func_name) {
    current_func_name = func_name;
    auto emplace_reuslt = time_consuming.emplace(
        func_name, std::make_pair(high_resolution_clock::duration::zero(),
                                  high_resolution_clock::now()));
    if (emplace_reuslt.second) {
        entries.push_back(func_name);
    } else {
        emplace_reuslt.first->second.second = high_resolution_clock::now();
    }
}

void Timer::pause_last_and_start_next(std::string func_name) {
    pause();
    start(func_name);
}

void Timer::pause() {
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = end_time - time_consuming.at(current_func_name).second;
    time_consuming.at(current_func_name).first += elapsed_time;
}

void Timer::pause(std::string func_name) {
    auto end_time = high_resolution_clock::now();
    auto elapsed_time = end_time - time_consuming.at(func_name).second;
    time_consuming.at(func_name).first += elapsed_time;
}

void Timer::reset() {
    entries.clear();
    time_consuming.clear();
}

void Timer::print() {
    std::size_t max_length = 0;
    for (auto& name : entries) {
        max_length = max_length < name.size() ? name.size() : max_length;
    }

    std::cout << '+' << std::setfill('-') << std::setw(max_length + 16) << '-'
              << std::setfill(' ');
    std::cout << "+\n"
              << std::left << std::setw(static_cast<int>(max_length + 17))
              << "| Time consumption";
    std::cout << "|\n";
    std::cout << '+' << std::setfill('-') << std::setw(max_length + 2) << '-'
              << std::setfill(' ');
    std::cout << "+-------------+\n";
    for (const auto& name : entries) {
        auto& time_used = time_consuming.at(name).first;
        std::cout << std::left << "| "
                  << std::setw(static_cast<int>(max_length)) << name << " | "
                  << std::right << std::setw(9)
                  << duration<double, milliseconds::period>(time_used).count()
                  << "ms |\n";
    }
    std::cout << '+' << std::setfill('-') << std::setw(max_length + 2) << '-'
              << std::setfill(' ');
    std::cout << "+-------------+\n";
}

Timer& Timer::get_timer() {
    static Timer timer{};
    return timer;
}

#endif  // ZQ_TIMER_IMPLEMENTATION
#endif  // ZQ_TIMER
