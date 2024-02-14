#pragma once
#include <vector>

class BigFloat {
private:
    std::vector<int> integerPart, fractionalPart;
    char sign;

    void addToIntPart(const std::vector<int> &digits);
    static void extend_vec(std::vector<int>& v, size_t len);
    static std::vector<int> naive_mul(const std::vector<int>& x, const std::vector<int>& y);
    static std::vector<int> karatsuba_mul(const std::vector<int>& x, const std::vector<int>& y);
    static void finalize(std::vector<int>& res);
    static std::vector<int> mult(std::vector<int>& x, std::vector<int>& y);
    static std::vector<int> sum(const std::vector<int> &x, const std::vector<int> &y);
    static std::vector<int> substract(const std::vector<int> &x, const std::vector<int> &y);
    void addToFracPart(const std::vector<int> &digits);
    void normalise();
    static void divideByTwo(std::vector<int> &digits, bool removeZeroes);
    void divideByTwo();

public:
    BigFloat& inverseSign();
    friend BigFloat abs(const BigFloat &x);
    explicit BigFloat(const std::vector<int> &intPart, const std::vector<int> &fracPart, char sign_);
    explicit BigFloat(int x);
    explicit BigFloat(const char* x);
    [[nodiscard]] std::string toString(size_t precision) const;

    friend BigFloat operator+(const BigFloat &a, const BigFloat &b);
    friend BigFloat operator-(const BigFloat &a, const BigFloat &b);
    bool operator >= (const BigFloat& other) const;
    bool operator < (const BigFloat& other) const;
    bool operator == (const BigFloat& other) const;
    bool operator <= (const BigFloat& other) const;
    bool operator != (const BigFloat& other) const;
    bool operator > (const BigFloat& other) const;
    BigFloat& operator -();
    friend BigFloat operator*(const BigFloat &a, const BigFloat &b);
    friend BigFloat operator/(const BigFloat &a, const BigFloat &b);
    friend std::ostream& operator << (std::ostream &in, BigFloat &x);
    friend void display(const BigFloat&x, size_t precision);
    explicit operator double() const;

};

BigFloat operator""_bf(const char *s);
