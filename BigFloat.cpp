#include <iostream>
#include <vector>
#include <algorithm>
#include <cstring>
#include <numeric>
#include "BigFloat.h"


void BigFloat::addToIntPart(const std::vector<int> &digits) {
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
void BigFloat::extend_vec(std::vector<int>& v, size_t len) {
    while (len & (len - 1)) {
        ++len;
    }

    v.resize(len);
}
std::vector<int> BigFloat::naive_mul(const std::vector<int>& x, const std::vector<int>& y) {
    auto len = x.size();
    std::vector<int> res(2 * len, 0);

    for (auto i = 0; i < len; ++i) {
        for (auto j = 0; j < len; ++j) {
            res[i + j] += x[i] * y[j];
        }
    }

    return res;
}
std::vector<int> BigFloat::karatsuba_mul(const std::vector<int>& x, const std::vector<int>& y) {
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
void BigFloat::finalize(std::vector<int>& res) {
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
std::vector<int> BigFloat::mult(std::vector<int>& x, std::vector<int>& y) {
    extend_vec(x, y.size());
    extend_vec(y, x.size());

    std::vector<int> res = karatsuba_mul(x, y);
    finalize(res);
    return res;
}
std::vector<int> BigFloat::sum(const std::vector<int> &x, const std::vector<int> &y) {
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
std::vector<int> BigFloat::substract(const std::vector<int> &x, const std::vector<int> &y) {
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
void BigFloat::addToFracPart(const std::vector<int> &digits) {
    int carry = 0;
    for(int i = 0; i < std::max(fractionalPart.size(), digits.size()); ++i) {
        int anum = i < fractionalPart.size() && fractionalPart[i];
        int bnum = i < digits.size() && digits[i];
        fractionalPart[i] = (anum + bnum + carry) % base;
        carry = (anum + bnum + carry) / base;
    }
    if(carry) fractionalPart.push_back(carry);
}
void BigFloat::normalise() {
    if(fractionalPart.size() > BigFloat::sizeOfFracPart) {
        addToIntPart(std::vector<int>(fractionalPart.begin() + BigFloat::sizeOfFracPart, fractionalPart.end()));
        fractionalPart.resize(BigFloat::sizeOfFracPart);
    }
}
void BigFloat::divideByTwo(std::vector<int> &digits, bool remove_leading_zeroes) {
    //https://ru.wikipedia.org/wiki/%D0%94%D0%B5%D0%BB%D0%B5%D0%BD%D0%B8%D0%B5_%D0%BD%D0%B0_%D0%B4%D0%B2%D0%B0
    std::vector<int> res;
    int last_nonzero = 0;
    for(int i = 0; i < digits.size() - 1; ++i) {
        if(digits[i+1] % 2 == 0) {
            res.push_back(digits[i] / 2);
            last_nonzero = digits[i] / 2 != 0 ? i : last_nonzero;
        }
        else{
            res.push_back(5 + digits[i] / 2);
            last_nonzero = i;
        }
    }
    res.push_back(digits.back() / 2);
    if(digits[digits.size() - 1] / 2 != 0) {
        last_nonzero = digits.size() - 1;
    }
    if(remove_leading_zeroes)
        res.resize(last_nonzero + 1);
    std::swap(digits, res);
}
void BigFloat::divideByTwo() {
    std::vector<int> res;
    bool divByTwo = true;
    if(integerPart[0] % 2) {
        divByTwo = false;

    }
    divideByTwo(integerPart, true);
    divideByTwo(fractionalPart, false);
    if(!divByTwo) {
        fractionalPart[BigFloat::sizeOfFracPart - 1] += 5;
    }
    if(fractionalPart[BigFloat::sizeOfFracPart - 1] >= base) {
        addToIntPart(std::vector<int>(1, fractionalPart[BigFloat::sizeOfFracPart - 1] / base));
        fractionalPart[BigFloat::sizeOfFracPart - 1] %= base;
    }
}
BigFloat& BigFloat::inverseSign() {
    sign = 1 - sign;
    return *this;
}
BigFloat abs(const BigFloat &x) {
    return BigFloat(x.integerPart, x.fractionalPart, 0);
}
BigFloat::BigFloat(const std::vector<int> &intPart, const std::vector<int> &fracPart, char sign_) : integerPart(intPart), fractionalPart(fracPart), sign(sign_){
    normalise();
}
BigFloat::BigFloat(int x) {
    if(x < 0) {
        sign = 1;
        x = -x;
    }
    else{
        sign = 0;
    }
    if(x == 0) {
        integerPart.push_back(0);
    }
    while(x) {
        integerPart.push_back(x % base);
        x /= base;
    }
    fractionalPart = std::vector<int> (BigFloat::sizeOfFracPart, 0);

}
BigFloat::BigFloat(const char* x) {
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
    for(int i = 0; i <= BigFloat::sizeOfFracPart - digitsCnt; ++i) digitsFrac.push_back(0);
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
BigFloat operator+(const BigFloat &a, const BigFloat &b) {
    if(a.sign == 0 && b.sign == 1) { //Полож +- отриц
        BigFloat newB = b;
        newB.inverseSign();
        return a - newB;
    }
    if(a.sign == 1 && b.sign == 0) { //Отриц + полож
        BigFloat newA = a;
        newA.inverseSign();
        BigFloat res = b - newA;
        if(b > newA) {
            res.sign = 0;
            return res;
        }
        else{
            res.sign = 1;
            return res;
        }
    }
    std::vector<int> resInt = BigFloat::sum(a.integerPart, b.integerPart);
    std::vector<int> resFrac = BigFloat::sum(a.fractionalPart, b.fractionalPart);
    if(resFrac.size() > BigFloat::sizeOfFracPart) {
        resInt = BigFloat::sum(resInt, std::vector<int>(resFrac.begin() + BigFloat::sizeOfFracPart, resFrac.end()));
        resFrac.resize(BigFloat::sizeOfFracPart);
    }
    return BigFloat(resInt, resFrac, a.sign);
}
BigFloat operator-(const BigFloat &a, const BigFloat &b) {
    if(a.sign != b.sign) {
        if(a.sign == 0) {
            BigFloat newB = b;
            newB.inverseSign();
            return a + newB;
        }
        else if(b.sign == 0) {
            BigFloat newA = a;
            newA.inverseSign();
            return (newA + b).inverseSign();
        }
    }
    if(a.sign == 1 && b.sign == 1) {
        BigFloat newB = b;
        newB.inverseSign();
        return newB + a;
    }
    std::vector<int> A = a.fractionalPart;
    A.insert(A.end(), a.integerPart.begin(), a.integerPart.end());
    std::vector<int> B = b.fractionalPart;
    B.insert(B.end(), b.integerPart.begin(), b.integerPart.end());
    std::vector<int> resVec;
    if(a < b) {
        std::swap(A, B);
        resVec = BigFloat::substract(A, B);
        int zeroesCnt = resVec.size() - 1;
        for(; zeroesCnt >= 0 && resVec[zeroesCnt] == 0; --zeroesCnt);
        return BigFloat(std::vector<int> (resVec.begin() + BigFloat::sizeOfFracPart, resVec.begin() + std::max(BigFloat::sizeOfFracPart + 1, zeroesCnt + 1)),
                        std::vector<int> (resVec.begin(), resVec.begin() + BigFloat::sizeOfFracPart),
                        1 - a.sign);
    }
    else {
        resVec = BigFloat::substract(A, B);
        int zeroesCnt = resVec.size() - 1;
        for(; zeroesCnt >= 0 && resVec[zeroesCnt] == 0; --zeroesCnt);
        return BigFloat(std::vector<int> (resVec.begin() + BigFloat::sizeOfFracPart, resVec.begin() + std::max(BigFloat::sizeOfFracPart + 1, zeroesCnt + 1)),
                        std::vector<int> (resVec.begin(), resVec.begin() + BigFloat::sizeOfFracPart),
                        a.sign);
    }

}
bool BigFloat::operator >= (const BigFloat& other) const
{
    return !(*this < other);
}
bool BigFloat::operator <= (const BigFloat& other) const
{
    return *this < other || *this == other;
}
bool BigFloat::operator < (const BigFloat& other) const
{
    if(sign == 1 && other.sign == 0) {
        return true;
    }
    if(sign == 0 && other.sign == 1) {
        return false;
    }
    if(integerPart.size() < other.integerPart.size()) {
        return 1 - sign;
    }
    if(integerPart.size() > other.integerPart.size()) {
        return sign;
    }
    for(int i = integerPart.size() - 1; i >= 0; --i) {
        if(integerPart[i] < other.integerPart[i]) {
            return 1 - sign;
        }
        if(integerPart[i] > other.integerPart[i]) {
            return sign;
        }
    }
    for(int i = BigFloat::sizeOfFracPart - 1; i >= 0; --i) {
        if(fractionalPart[i] < other.fractionalPart[i]) {
            return 1 - sign;
        }
        if(fractionalPart[i] > other.fractionalPart[i]) {
            return sign;
        }
    }
    return false;
}
bool BigFloat::operator == (const BigFloat& other) const
{
    for(int i = 0; i < BigFloat::sizeOfFracPart; ++i) {
        if(fractionalPart[i] != other.fractionalPart[i]) {
            return false;
        }
    }
    if(integerPart.back() == 0  && integerPart.size() > 1 || other.integerPart.back() == 0 && other.integerPart.size() > 1) {
        std::cerr << "ERROR! - Leading zeroes";
        std::cerr << this->toString(BigFloat::sizeOfFracPart) << ' ' << other.toString(BigFloat::sizeOfFracPart) << '\n';
    }
    if(integerPart.size() != other.integerPart.size())
        return false;
    for(int i = 0; i < integerPart.size(); ++i) {
        if(integerPart[i] != other.integerPart[i]) {
            return false;
        }
    }
    return sign == other.sign || integerPart == std::vector{0} && other.integerPart == std::vector{0};
}

bool BigFloat::operator!=(const BigFloat &other) const {
    return !(*this == other);
}
bool BigFloat::operator>(const BigFloat &other) const {
    return !(*this <= other);
}
BigFloat& BigFloat::operator -() {
    sign = 1 - sign;
    return *this;
}
BigFloat operator*(const BigFloat &a, const BigFloat &b) {
    std::vector<int> A = a.fractionalPart;
    A.insert(A.end(), a.integerPart.begin(), a.integerPart.end());
    std::vector<int> B = b.fractionalPart;
    B.insert(B.end(), b.integerPart.begin(), b.integerPart.end());
    if(std::accumulate(A.begin(), A.end(), 0) == 0 || std::accumulate(B.begin(), B.end(), 0) == 0) {
        return BigFloat("0.");
    }
    std::vector<int> resVec = BigFloat::mult(A, B);
    int firstNonZero = resVec.size() - 1;
    for(; firstNonZero >= 0 && resVec[firstNonZero] == 0; --firstNonZero);
    BigFloat res(std::vector<int>(resVec.begin() + 2 * BigFloat::sizeOfFracPart, resVec.begin() + std::max(2 * BigFloat::sizeOfFracPart + 1, firstNonZero + 1)),
                 std::vector<int>(resVec.begin() + BigFloat::sizeOfFracPart, resVec.begin() + 2 * BigFloat::sizeOfFracPart),
                 a.sign ^ b.sign);
    return res;
}
BigFloat operator/(const BigFloat &x, const BigFloat &y) {
    if(std::accumulate(y.fractionalPart.begin(), y.fractionalPart.end(), 0) == 0 && std::accumulate(y.integerPart.begin(), y.integerPart.end(), 0) == 0) {
            std::cerr << "ERROR! - Division by zero";
            throw std::runtime_error("Division by zero");
    }
    if(std::accumulate(x.fractionalPart.begin(), x.fractionalPart.end(), 0) == 0 && std::accumulate(x.integerPart.begin(), x.integerPart.end(), 0) == 0) {
        return BigFloat("0.");
    }
    BigFloat a = x;
    BigFloat b = y;
    if(a.sign == 1) {
        a.inverseSign();
    }
    if(b.sign == 1) {
        b.inverseSign();
    }
    BigFloat L("0.");
    BigFloat R = a;
    BigFloat eps(std::vector<int>(1, 0), std::vector<int> (BigFloat::sizeOfFracPart, 0), 0);
    static_assert(BigFloat::sizeOfFracPart >= 103, "BigFloat::sizeOfFracPart must be at least 103"); //103 is a magic precision for 100-digit precise PI calculation
    //Precision (number a bit bigger than 1, e.x. 20)
    eps.fractionalPart[10] = 1;
    while (R - L >= eps) {
        BigFloat M = (L + R);
        M.divideByTwo();
        if (M * b < a) {
            L = M;
        } else {
            R = M;
        }
    }
    R.sign = x.sign ^ y.sign;
    return R;
}

std::string BigFloat::toString(size_t precision) const {
    std::string res;
    if(this->sign) {
        res += '-';
    }
    for(int i = this->integerPart.size() - 1; i >= 0; --i) {
        res += static_cast<char>('0' + this->integerPart[i]);
    }
    res += '.';
    for(int i = this->fractionalPart.size() - 1; i >= 0 && i >= this->fractionalPart.size() - std::min(precision, this->fractionalPart.size()); --i) {
        res += static_cast<char>('0' + this->fractionalPart[i]);
    }
    return res;
}
void display(const BigFloat &x, size_t precision) {
    std::cout << x.toString(precision) << '\n';
}
std::ostream& operator<<(std::ostream &out, BigFloat &x) {
    out << x.toString(BigFloat::sizeOfFracPart);
    return out;
}

BigFloat::operator double() const {
    return std::stod(toString(10));
}

BigFloat operator""_bf(const char *s) {
    return BigFloat{s};
}
