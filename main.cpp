#include <chrono>
#include <thread>
#include <iostream>
#include "BigFloat.h"

//https://habr.com/ru/articles/443998/
void CalcPi(BigFloat &pi, const int k, const BigFloat &bs) {
    static BigFloat one = 1._bf;
    static BigFloat two = 2._bf;
    static BigFloat four = 4._bf;
    BigFloat res = ((four / BigFloat(8*k+1)) -
                    (two / BigFloat(8*k+4)) -
                    (one / BigFloat(8*k+5)) -
                    (one / BigFloat(8*k+6))) / bs;
    pi = pi + res;
}

int main(int argc, char**argv) {
    std::cout << "Enter the precision of pi you want to get (number of digits after dot)" << std::endl;
    int n;
    std::cin >> n;
    auto start = std::chrono::high_resolution_clock::now();
    BigFloat pi = 0._bf;
    BigFloat curBs = 1._bf;
    std::vector<std::thread> threads(80);

    for (int k = 0; k < 80; k++) {
        threads[k] = std::thread(CalcPi, std::ref(pi), k, curBs);
        curBs = curBs * BigFloat(16);
    }
    for (auto &thread : threads)
        thread.join();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(std::chrono::high_resolution_clock::now() - start);
    std::cout << "Pi with needed precision: ";
    display(pi, n);
    std::cout << "Total time (in ms) " << duration.count() << '\n';
    return 0;
}