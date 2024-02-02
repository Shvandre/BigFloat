#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <thread>
#include "catch2/catch_session.hpp"
#include <catch2/catch_test_macros.hpp>

#define all(x) x.begin(), x.end()
const int sizeOfFracPart = 128; //Must be a power of 2 because of Karatsuba algorithm
const int base = 10;

class BigFloat {
private:
    std::vector<int> integerPart, fractionalPart;
    char sign = 0; // 1 - negative, 0 - positive
    //int precision; //Means how many digits after dot we need to store
    void addToIntPart(const std::vector<int> &digits) {
        int carry = 0;
        for(int i = 0; i < std::max(integerPart.size(), digits.size()); ++i) {
            int anum = i < integerPart.size() && integerPart[i];
            int bnum = i < digits.size() && digits[i];
            integerPart[i] = (anum + bnum + carry) % base;
            carry = (anum + bnum + carry) / base;
        }
        if(carry) integerPart.push_back(carry);
    }
    //Following three functions are taken from https://habr.com/ru/articles/262705/
    static void extend_vec(std::vector<int>& v, size_t len) {
        while (len & (len - 1)) {
            ++len;
        }

        v.resize(len);
    }
    static std::vector<int> naive_mul(const std::vector<int>& x, const std::vector<int>& y) {
        auto len = x.size();
        std::vector<int> res(2 * len, 0);

        for (auto i = 0; i < len; ++i) {
            for (auto j = 0; j < len; ++j) {
                res[i + j] += x[i] * y[j];
            }
        }

        return res;
    }
    static std::vector<int> karatsuba_mul(const std::vector<int>& x, const std::vector<int>& y) {
        auto len = x.size();
        std::vector<int> res(2 * len);

        if (len <= 100) { //This constant means that we will use native mult for small numbers (because of better constant)
            return naive_mul(x, y);
        }

        auto k = len / 2;

        std::vector<int> Xr {x.begin(), x.begin() + k};
        std::vector<int> Xl {x.begin() + k, x.end()};
        std::vector<int> Yr {y.begin(), y.begin() + k};
        std::vector<int> Yl {y.begin() + k, y.end()};

        std::vector<int> P1 = karatsuba_mul(Xl, Yl);
        std::vector<int> P2 = karatsuba_mul(Xr, Yr);

        std::vector<int> Xlr(k);
        std::vector<int> Ylr(k);

        for (int i = 0; i < k; ++i) {
            Xlr[i] = Xl[i] + Xr[i];
            Ylr[i] = Yl[i] + Yr[i];
        }

        std::vector<int> P3 = karatsuba_mul(Xlr, Ylr);

        for (auto i = 0; i < len; ++i) {
            P3[i] -= P2[i] + P1[i];
        }

        for (auto i = 0; i < len; ++i) {
            res[i] = P2[i];
        }

        for (auto i = len; i < 2 * len; ++i) {
            res[i] = P1[i - len];
        }

        for (auto i = k; i < len + k; ++i) {
            res[i] += P3[i - k];
        }

        return res;
    }
    static void finalize(std::vector<int>& res) {
        int i = 0;
        for (; i < res.size() - 1; ++i) {
            res[i + 1] += res[i] / base;
            res[i] %= base;
        }
        while(res[i] >= base) {
            res.push_back(res[i] / base);
            res[i] %= base;
            ++i;
        }
    }
    static std::vector<int> mult(std::vector<int>& x, std::vector<int>& y) {
        extend_vec(x, y.size());
        extend_vec(y, x.size());

        std::vector<int> res = karatsuba_mul(x, y);
        finalize(res);
        return res;
    }
    static std::vector<int> sum(const std::vector<int> &x, const std::vector<int> &y) {
        std::vector<int> res(std::max(x.size(), y.size()));
        int carry = 0;
        for(int i = 0; i < std::max(x.size(), y.size()); ++i) {
            int anum = i < x.size() ? x[i] : 0;
            int bnum = i < y.size() ? y[i] : 0;
            res[i] = (anum + bnum + carry) % base;
            carry = (anum + bnum + carry) / base;
        }
        if(carry) res.push_back(carry);
        return res;
    }
    static std::vector<int> substract(const std::vector<int> &x, const std::vector<int> &y) {
        std::vector<int> res(std::max(x.size(), y.size()));
        int loan = 0;
        for(int i = 0; i < std::max(x.size(), y.size()); ++i) {
            int anum = i < x.size() ? x[i] : 0;
            int bnum = i < y.size() ? y[i] : 0;
            if(anum - loan >= bnum) {
                res[i] = anum - loan - bnum;
                loan = 0;
            } else {
                res[i] = anum - loan + base - bnum;
                loan = 1;
            }
        }
        if(loan) {
            return substract(y, x);
        }
        return res;
    }
    void addToFracPart(const std::vector<int> &digits) {
        int carry = 0;
        for(int i = 0; i < std::max(fractionalPart.size(), digits.size()); ++i) {
            int anum = i < fractionalPart.size() && fractionalPart[i];
            int bnum = i < digits.size() && digits[i];
            fractionalPart[i] = (anum + bnum + carry) % base;
            carry = (anum + bnum + carry) / base;
        }
        if(carry) fractionalPart.push_back(carry);
    }
    void normalise() {
        if(fractionalPart.size() > sizeOfFracPart) {
            addToIntPart(std::vector<int>(fractionalPart.begin() + sizeOfFracPart, fractionalPart.end()));
            fractionalPart.resize(sizeOfFracPart);
        }
    }
    static void divideByTwo(std::vector<int> &digits) {
        //https://ru.wikipedia.org/wiki/%D0%94%D0%B5%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BD%D0%B0_%D0%B4%D0%B2%D0%B0
        std::vector<int> res;
        for(int i = 0; i < digits.size() - 1; ++i) {
            if(digits[i+1] % 2 == 0) {
                res.push_back(digits[i] / 2);
            }
            else{
                res.push_back(5 + digits[i] / 2);
            }
        }
        res.push_back(digits.back() / 2);
        std::swap(digits, res);
    }
    void divideByTwo() {
        std::vector<int> res;
        bool divByTwo = true;
        if(integerPart[0] % 2) {
            divByTwo = false;

        }
        divideByTwo(integerPart);
        divideByTwo(fractionalPart);
        if(!divByTwo) {
            fractionalPart[sizeOfFracPart - 1] += 5;
        }
        if(fractionalPart[sizeOfFracPart - 1] >= base) {
            addToIntPart(std::vector<int>(1, fractionalPart[sizeOfFracPart - 1] / base));
            fractionalPart[sizeOfFracPart - 1] %= base;
        }
    }
public:
    void inverseSign() {
        sign = 1 - sign;
    }
    BigFloat abs(const BigFloat &x) {
        return BigFloat(x.integerPart, x.fractionalPart, 0);
    }
    explicit BigFloat(const std::vector<int> &intPart, const std::vector<int> &fracPart, char sign_) {
        integerPart = intPart;
        fractionalPart = fracPart;
        sign = sign_;
    }
    explicit BigFloat(int x) {
        while(x) {
            integerPart.push_back(x % base);
            x /= base;
        }
        fractionalPart = std::vector<int> (sizeOfFracPart, 0);
    }
    explicit BigFloat(const char* x) {
        char* s = const_cast<char *>(x);
        if(s[0] == '-') {
            sign = 1;
            ++s;
        } else {
            sign = 0;
        }
        size_t n = strlen(s);
        std::vector<int> digitsInt, digitsFrac;
        digitsInt.reserve(n);
        digitsFrac.reserve(n);
        bool meetDot = false;
        int digitsCnt = n - (strchr(x, '.') - x);
        digitsFrac.reserve(digitsCnt);
        for(int i = 0; i <= sizeOfFracPart - digitsCnt; ++i) digitsFrac.push_back(0);
        for(int i = n - 1; i >= 0; --i) {
            if(s[i] == '.') {
                meetDot = true;
                continue;
            }
            if(!meetDot) {
                digitsFrac.push_back(s[i] - '0');
            } else {
                digitsInt.push_back(s[i] - '0');
            }
        }
        integerPart = digitsInt;
        fractionalPart = digitsFrac;
    }
    friend BigFloat operator+(const BigFloat &a, const BigFloat &b) {
        if(a.sign == 0 && b.sign == 1) { //Полож +- отриц
            BigFloat newB = b;
            newB.inverseSign();
            return a - b;
        }
        if(a.sign == 1 && b.sign == 0) { //Отриц + полож
            BigFloat newA = a;
            newA.inverseSign();
            return b - a;
        }
        std::vector<int> resInt = sum(a.integerPart, b.integerPart);
        std::vector<int> resFrac = sum(a.fractionalPart, b.fractionalPart);
        if(resFrac.size() > sizeOfFracPart) {
            resInt = sum(resInt, std::vector<int>(resFrac.begin() + sizeOfFracPart, resFrac.end()));
            resFrac.resize(sizeOfFracPart);
        }
        return BigFloat(resInt, resFrac, a.sign);
    }
    friend BigFloat operator-(const BigFloat &a, const BigFloat &b) {
        if(a.sign != b.sign) {
            return a + b;
        }
        if(a.sign == 1 && b.sign == 1) {
            BigFloat newB = b;
            newB.inverseSign();
            return b - a;
        }
        std::vector<int> A = a.fractionalPart;
        A.insert(A.end(), a.integerPart.begin(), a.integerPart.end());
        std::vector<int> B = b.fractionalPart;
        B.insert(B.end(), b.integerPart.begin(), b.integerPart.end());
        std::vector<int> resVec;
        if(a < b) {
            std::swap(A, B);
            resVec = substract(A, B);
            return BigFloat(std::vector<int> (resVec.begin() + sizeOfFracPart, resVec.end()),
                            std::vector<int> (resVec.begin(), resVec.begin() + sizeOfFracPart),
                            1 - a.sign);
        }
        else {
            resVec = substract(A, B);
            return BigFloat(std::vector<int> (resVec.begin() + sizeOfFracPart, resVec.end()),
                            std::vector<int> (resVec.begin(), resVec.begin() + sizeOfFracPart),
                            a.sign);
        }

    }
    bool operator >= (const BigFloat& other) const
    {
        return (*this-other).sign == 0;
    }
    bool operator < (const BigFloat& other) const
    {
        if(sign == 1 && other.sign == 0) {
            return true;
        }
        if(sign == 0 && other.sign == 1) {
            return false;
        }
        if(integerPart.size() < other.integerPart.size()) {
            return true;
        }
        if(integerPart.size() > other.integerPart.size()) {
            return false;
        }
        for(int i = integerPart.size() - 1; i >= 0; --i) {
            if(integerPart[i] < other.integerPart[i]) {
                return true;
            }
            if(integerPart[i] > other.integerPart[i]) {
                return false;
            }
        }
        for(int i = sizeOfFracPart - 1; i >= 0; --i) {
            if(fractionalPart[i] < other.fractionalPart[i]) {
                return true;
            }
            if(fractionalPart[i] > other.fractionalPart[i]) {
                return false;
            }
        }
        return false;
    }
    bool operator == (const BigFloat& other) const
    {
        for(int i = 0; i < sizeOfFracPart; ++i) {
            if(fractionalPart[i] != other.fractionalPart[i]) {
                return false;
            }
        }
        if(integerPart.back() == 0 || other.integerPart.back() == 0) {
            std::cerr << "Leading zeroes";
            throw std::exception();
        }
        if(integerPart.size() != other.integerPart.size())
            return false;
        for(int i = 0; i < integerPart.size(); ++i) {
            if(integerPart[i] != other.integerPart[i]) {
                return false;
            }
        }
        return sign == other.sign || integerPart.empty() && other.integerPart.empty();
    }
    BigFloat& operator -() {
        sign = 1 - sign;
        return *this;
    }
    friend BigFloat operator*(const BigFloat &a, const BigFloat &b) {
        std::vector<int> A = a.fractionalPart;
        A.insert(A.end(), a.integerPart.begin(), a.integerPart.end());
        std::vector<int> B = b.fractionalPart;
        B.insert(B.end(), b.integerPart.begin(), b.integerPart.end());
        std::vector<int> resVec = mult(A, B);
        int firstNonZero = resVec.size() - 1;
        for(; firstNonZero >= 0 && resVec[firstNonZero] == 0; --firstNonZero);
        BigFloat res(std::vector<int>(resVec.begin() + 2 * sizeOfFracPart, resVec.begin() + std::max(2 * sizeOfFracPart + 1, firstNonZero + 1)),
                    std::vector<int>(resVec.begin() + sizeOfFracPart, resVec.begin() + 2 * sizeOfFracPart),
                    a.sign ^ b.sign);
        return res;
    }
    friend BigFloat operator/(const BigFloat &a, const BigFloat &b) {
        BigFloat L("0.");
        BigFloat R = a;
        BigFloat eps("0.0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001");
        while (R - L >= eps) {
            BigFloat M = (L + R);
            M.divideByTwo();
            if (M * b < a) {
                L = M;
            } else {
                R = M;
            }
        }
        return R;
    }
    friend void display(const BigFloat&x) {
        for(int i = x.integerPart.size() - 1; i >= 0; --i) {
            std::cout << x.integerPart[i];
        }
        std::cout << '.';
        for(int i = x.fractionalPart.size() - 1; i >= 0; --i) {
            std::cout << x.fractionalPart[i];
        }
        std::cout << '\n';
    }
};



BigFloat operator""_bf(const char *s) {
    return BigFloat{s};
}

TEST_CASE( "All", "[BigFloat]" ) {
    BigFloat a = 1.96_bf;
    BigFloat b = 2.785_bf;
    display(a);
    display(b);
    display(a / b);
    REQUIRE( (a / b - 0.7037701975216623395681381225585937500000_bf) < 0.004_bf);
    display(a * b);
    REQUIRE( a * b == 524.8586_bf );
    display(a + b);
    REQUIRE( a + b == 269.745_bf );
    display(a - b);
    REQUIRE( a - b == -265.825_bf );
}

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
 //   int result = Catch::Session().run( argc, argv );
    clock_t start = clock();
    BigFloat pi = 0._bf;
    BigFloat bs = 1._bf;
    std::vector<std::thread> threads;
    for (int k = 0; k < 80; k++) {
        std::thread calcThread(CalcPi, std::ref(pi), k, bs);
        bs = bs * BigFloat(16);
        threads.push_back(std::move(calcThread));
    }
    for (auto &thread : threads) {
        thread.join();
    }
    clock_t alltime = clock() - start;
    display(pi);
    std::cout << alltime * 1000 / CLOCKS_PER_SEC;
    return 0;
}
