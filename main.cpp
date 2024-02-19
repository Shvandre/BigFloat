#include <chrono>
#include <thread>
#include <mutex>
#include <iostream>
#include "BigFloat.h"

//https://habr.com/ru/articles/443998/
void CalcPi(BigFloat &pi, const int k_start, const int k_finish, const BigFloat &bs) {
    static BigFloat one = 1._bf;
    static BigFloat two = 2._bf;
    static BigFloat four = 4._bf;
    static std::mutex m;
    BigFloat base = bs;
    BigFloat res = 0.0_bf;
    for(int i = k_start; i < k_finish; ++i) {
        res = res + ((four / BigFloat(8 * i + 1)) -
                        (two / BigFloat(8 * i + 4)) -
                        (one / BigFloat(8 * i + 5)) -
                        (one / BigFloat(8 * i + 6))) / base;
        base = base * BigFloat(16);
    }
    m.lock();
    pi = pi + res;
    m.unlock();
}

int main(int argc, char**argv) {
    std::cout << "Enter the precision of pi you want to get (number of digits after dot)" << std::endl;
    int n;
    std::cin >> n;
    int signs = n / 16;
    auto start = std::chrono::high_resolution_clock::now();
    BigFloat pi = 0._bf;
    BigFloat curBs = 1._bf;
    //Let's use 16 threads to calculate pi
    std::vector<std::thread> threads;

    for (int k = 0; k <= n - 15; k++) {
        if (k % signs == 0)
            threads.emplace_back(CalcPi, std::ref(pi), k, k+signs, curBs);
        curBs = curBs * BigFloat(16);
    }
    for (auto &thread: threads)
        thread.join();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Pi with needed precision: ";
    display(pi, n);
    std::cout << "Total time (in ms) " << duration.count() << '\n';
    return 0;
}