#pragma once

#include <cstring>

#include <algorithm>
#include <iostream>
#include <string>
#include <vector>

class BigInteger {
  private:
    enum Signs { MINUS = -1, PLUS = 1 };
    std::vector<long long> number_;
    static const int base_ = 1e9;
    static const int b_cnt_ = 9;

    int sign_;

    friend void swap(BigInteger& number1, BigInteger& number2);

    static void borrow(std::vector<long long>& number, int to);

    static void cut(std::vector<long long>& number);

    void fix(bool is_changed = true);

    static std::pair<int, BigInteger> DivBinsearch(
        std::vector<long long>& divisible,
        const std::vector<long long>& divider, int lower_bit, int higher_bit);
    void clear();

    void subtrOrAdd(const BigInteger& number, int operation);

    void subtrOrAddSameSign(const BigInteger& number);

    void subtrOrAddPosSign(const BigInteger& number);

    void subtrOrAddNegSign(const BigInteger& number);

    void modOrDiv(const BigInteger& number, bool to_mod);

  public:
    BigInteger(int number);
    BigInteger(const char* number, size_t sz = 0);
    BigInteger();
    explicit operator bool() const;

    std::string toString() const;
    const std::vector<long long>& getData() const;
    int getSize() const;
    int getSign() const;
    void setSign(int sign) {
        sign_ = sign;
    }
    static int getBase();

    BigInteger& operator+=(const BigInteger& number);
    BigInteger& operator-=(const BigInteger& number);
    BigInteger& operator*=(const BigInteger& number);
    BigInteger& operator/=(const BigInteger& number);
    BigInteger& operator%=(const BigInteger& number);
    BigInteger operator-() const;
    BigInteger operator+() const;
    BigInteger& operator++();
    BigInteger operator++(int);
    BigInteger& operator--();
    BigInteger operator--(int);
};

BigInteger operator""_bi(unsigned long long number);
BigInteger operator""_bi(const char* number, unsigned long sz);

bool operator<(const BigInteger& a, const BigInteger& b);
bool operator==(const BigInteger& a, const BigInteger& b);
bool operator>(const BigInteger& a, const BigInteger& b);
bool operator<=(const BigInteger& a, const BigInteger& b);
bool operator>=(const BigInteger& a, const BigInteger& b);
bool operator!=(const BigInteger& a, const BigInteger& b);
BigInteger operator+(const BigInteger& a, const BigInteger& b);
BigInteger operator-(const BigInteger& a, const BigInteger& b);
BigInteger operator*(const BigInteger& a, const BigInteger& b);
BigInteger operator/(const BigInteger& a, const BigInteger& b);
BigInteger operator%(const BigInteger& a, const BigInteger& b);

class Rational {
  private:
    BigInteger numerator_;
    BigInteger denominator_;

    void shorten();

    static BigInteger gcd(const BigInteger& a, const BigInteger& b);

    void swap(Rational& frac);

    void subtrOrAdd(const Rational& frac, bool to_substr);

  public:
    Rational(const BigInteger& num, const BigInteger& den = 1);
    Rational();
    Rational(int number);
    friend bool operator<(const Rational& a, const Rational& b);
    Rational operator-() const;
    Rational operator+() const;
    Rational& operator+=(const Rational& frac);
    Rational& operator*=(const Rational& frac);
    Rational& operator/=(const Rational& frac);
    Rational& operator-=(const Rational& frac);
    explicit operator double();

    std::pair<BigInteger, BigInteger> getData() const {
        return std::make_pair(numerator_, denominator_);
    }

    std::string toString() const;

    std::string asDecimal(size_t precision = 20) const;
};

bool operator==(const Rational& a, const Rational& b);
bool operator!=(const Rational& a, const Rational& b);
bool operator>(const Rational& a, const Rational& b);
bool operator>=(const Rational& a, const Rational& b);
bool operator<=(const Rational& a, const Rational& b);
Rational operator+(const Rational& a, const Rational& b);
Rational operator-(const Rational& a, const Rational& b);
Rational operator*(const Rational& a, const Rational& b);
Rational operator/(const Rational& a, const Rational& b);