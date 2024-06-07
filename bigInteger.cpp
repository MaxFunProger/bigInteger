#include "biginteger.h"

BigInteger::BigInteger(int number) : number_({}), sign_(PLUS) {
    if (number == 0) {
        number_ = {0};
        return;
    }
    if (number < 0) {
        sign_ = MINUS;
        number = -number;
    }
    while (number > 0) {
        number_.push_back(number % base_);
        number /= base_;
    }
    cut(number_);
}

BigInteger::BigInteger(const char* number, size_t sz)
    : number_({}), sign_(PLUS) {
    if (sz == 0 && (strcmp(number, "") != 0)) {
        sz = strlen(number);
    }
    std::string tmp, result;
    if (number[0] == '-') {
        sign_ = MINUS;
    } else {
        sign_ = PLUS;
    }
    for (int i = sz - 1; i >= 0;) {
        for (int j = 0; j < (b_cnt_ < i + 1 ? b_cnt_ : i + 1); ++j) {
            tmp.push_back(number[i - j]);
        }
        if (tmp.back() == '-' || tmp.back() == '+') {
            tmp.pop_back();
        }
        if (!tmp.empty()) {
            std::reverse(tmp.begin(), tmp.end());
            number_.push_back(std::stoi(tmp));
            tmp.clear();
        }
        i -= (b_cnt_ < i + 1 ? b_cnt_ : i + 1);
    }
    cut(number_);
}

BigInteger::BigInteger() : number_({}), sign_(PLUS) {}

void swap(BigInteger& number1, BigInteger& number2) {
    std::swap(number1.number_, number2.number_);
    std::swap(number1.sign_, number2.sign_);
}

BigInteger::operator bool() const {
    return !(number_.empty() || (number_.size() == 1 && number_[0] == 0));
}

std::string BigInteger::toString() const {
    std::string result;
    if (sign_ == MINUS) {
        result.push_back('-');
    }
    for (int i = number_.size() - 1; i >= 0; --i) {
        if (result.empty() && number_[i] == 0 && i > 0) {
            continue;
        }
        std::string tmp = std::to_string(number_[i]);
        if (!result.empty() && result != "-") {
            int sz = tmp.size();
            for (int i = 0; i < b_cnt_ - sz; ++i) {
                result.push_back('0');
            }
        }
        result += tmp;
    }
    return result;
}

const std::vector<long long>& BigInteger::getData() const {
    return number_;
}

int BigInteger::getSize() const {
    return number_.size();
}

int BigInteger::getSign() const {
    return sign_;
}

int BigInteger::getBase() {
    return base_;
}

/*BigInteger& BigInteger::operator=(BigInteger number) {
    swap(*this, number);
    return *this;
}*/

BigInteger operator""_bi(unsigned long long number) {
    return BigInteger(number);
}

BigInteger operator""_bi(const char* number, unsigned long sz) {
    return BigInteger(number, sz);
}

bool operator==(const BigInteger& a, const BigInteger& b) {
    return a.getData() == b.getData() && a.getSign() == b.getSign();
}

bool operator<(const BigInteger& a, const BigInteger& b) {
    if (a.getSign() != b.getSign()) {
        return a.getSign() < b.getSign();
    }
    bool flag = false;
    if (a.getSign() == 1) {
        if (a.getSize() != b.getSize()) {
            return a.getSize() < b.getSize();
        }
    } else {
        if (a.getSize() != b.getSize()) {
            return a.getSize() > b.getSize();
        }
        flag = true;
    }
    for (int i = a.getSize() - 1; i >= 0; --i) {
        if (a.getData()[i] > b.getData()[i]) {
            return flag;
        }
        if (a.getData()[i] < b.getData()[i]) {
            return !flag;
        }
    }
    return false;
}

bool operator>(const BigInteger& a, const BigInteger& b) {
    return b < a;
}

bool operator<=(const BigInteger& a, const BigInteger& b) {
    return !(a > b);
}

bool operator>=(const BigInteger& a, const BigInteger& b) {
    return !(a < b);
}

bool operator!=(const BigInteger& a, const BigInteger& b) {
    return !(a == b);
}

BigInteger BigInteger::operator-() const {
    BigInteger result = *this;
    if (!(result.number_.size() == 1 && result.number_[0] == 0)) {
        result.sign_ *= MINUS;
    }
    return result;
}

BigInteger BigInteger::operator+() const {
    return *this;
}

void BigInteger::borrow(std::vector<long long>& number, int to) {
    int index = to + 1;
    while (number[index] == 0) {
        ++index;
    }
    --number[index];
    --index;
    while (index > to) {
        number[index] += base_ - 1;
        --index;
    }
    number[index] += base_;
}

void BigInteger::cut(std::vector<long long>& number) {
    size_t size = number.size();
    while (size > 1) {
        if (number[size - 1] == 0) {
            number.pop_back();
            --size;
        } else {
            break;
        }
    }
}

void BigInteger::clear() {
    number_.clear();
    sign_ = PLUS;
}

std::ostream& operator<<(std::ostream& out, const BigInteger& number) {
    out << number.toString();
    return out;
}

std::istream& operator>>(std::istream& in, BigInteger& number) {
    std::string s;
    in >> s;
    number = BigInteger(s.data(), s.size());
    return in;
}

void BigInteger::fix(bool is_changed) {
    while (number_.back() == 0 && number_.size() > 1) {
        number_.pop_back();
    }
    if (number_.back() < 0) {
        sign_ = MINUS;
        for (size_t i = 0; i < number_.size(); ++i) {
            number_[i] *= -1;
        }
    } else if (!is_changed) {
        sign_ = PLUS;
    }
    for (size_t i = 0; i < number_.size(); ++i) {
        if (number_[i] < 0) {
            number_[i] += base_;
            --number_[i + 1];
        }
    }
    for (size_t i = 0; i < number_.size(); ++i) {
        if (number_[i] >= base_) {
            int tmp = number_[i] / base_;
            number_[i] %= base_;
            if (i == number_.size()) {
                number_.push_back(tmp);
                break;
            } else {
                number_[i + 1] += tmp;
            }
        }
    }
    while (number_.back() == 0 && number_.size() > 1) {
        number_.pop_back();
    }
}

void BigInteger::subtrOrAddSameSign(const BigInteger& number) {
    if (number_.size() >= number.number_.size()) {
        for (size_t i = 0; i < number.number_.size(); ++i) {
            number_[i] += number.number_[i];
        }
    } else {
        for (size_t i = 0; i < number_.size(); ++i) {
            number_[i] += number.number_[i];
        }
        for (size_t i = number_.size(); i < number.number_.size(); ++i) {
            number_.push_back(number.number_[i]);
        }
    }
}

void BigInteger::subtrOrAddPosSign(const BigInteger& number) {
    if (number_.size() >= number.number_.size()) {
        for (size_t i = 0; i < number.number_.size(); ++i) {
            number_[i] = number.number_[i] - number_[i];
        }
    } else {
        for (size_t i = 0; i < number_.size(); ++i) {
            number_[i] = number.number_[i] - number_[i];
        }
        for (size_t i = number_.size(); i < number.number_.size(); ++i) {
            number_.push_back(number.number_[i]);
        }
    }
}

void BigInteger::subtrOrAddNegSign(const BigInteger& number) {
    if (number_.size() >= number.number_.size()) {
        for (size_t i = 0; i < number.number_.size(); ++i) {
            number_[i] -= number.number_[i];
        }
    } else {
        for (size_t i = 0; i < number_.size(); ++i) {
            number_[i] -= number.number_[i];
        }
        for (size_t i = number_.size(); i < number.number_.size(); ++i) {
            number_.push_back(-number.number_[i]);
        }
    }
}

void BigInteger::subtrOrAdd(const BigInteger& number, int operation) {
    if (sign_ == number.sign_ * operation) {
        subtrOrAddSameSign(number);
    } else if (number.sign_ * operation == MINUS) {
        subtrOrAddNegSign(number);
    } else {
        subtrOrAddPosSign(number);
    }
    if (sign_ != number.sign_ * operation) {
        fix(false);
    } else {
        fix();
    }
}

BigInteger& BigInteger::operator+=(const BigInteger& number) {
    subtrOrAdd(number, 1);
    return *this;
}

BigInteger operator+(const BigInteger& a, const BigInteger& b) {
    BigInteger tmp = a;
    tmp += b;
    return tmp;
}

BigInteger operator-(const BigInteger& a, const BigInteger& b) {
    BigInteger output = a;
    return output -= b;
}

BigInteger& BigInteger::operator-=(const BigInteger& number) {
    subtrOrAdd(number, -1);
    return *this;
}

BigInteger& BigInteger::operator*=(const BigInteger& number) {
    BigInteger result(0);
    BigInteger term;
    int zeros_cnt = 1;
    long long next = 0, tmp;
    for (size_t i = 0; i < number.number_.size(); ++i) {
        for (size_t j = 0; j < number_.size(); ++j) {
            tmp = number.number_[i] * number_[j] + next;
            next = tmp / base_;
            term.number_.push_back(tmp % base_);
        }
        if (next != 0) {
            term.number_.push_back(next);
        }
        result += term;
        term.clear();
        for (int i = 0; i < zeros_cnt; ++i) {
            term.number_.push_back(0);  // increasing bit height
        }
        ++zeros_cnt;
        next = 0;
    }
    int sign = sign_ * number.sign_;
    *this = result;
    sign_ = sign;
    cut(number_);
    if (number_.size() == 1 && number_[0] == 0) {
        sign_ = PLUS;
    }
    return *this;
}

BigInteger operator*(const BigInteger& a, const BigInteger& b) {
    BigInteger tmp = a;
    tmp *= b;
    return tmp;
}

BigInteger& BigInteger::operator++() {
    *this += 1;
    return *this;
}

BigInteger BigInteger::operator++(int) {
    BigInteger tmp = *this;
    *this += 1;
    return tmp;
}

BigInteger& BigInteger::operator--() {
    *this -= 1;
    return *this;
}

BigInteger BigInteger::operator--(int) {
    BigInteger tmp = *this;
    *this -= 1;
    return tmp;
}

std::pair<int, BigInteger> BigInteger::DivBinsearch(
    std::vector<long long>& divisible, const std::vector<long long>& divider,
    int lower_bit, int higher_bit) {
    int left = 0, right = 1e9;
    BigInteger divis, divid;
    for (int i = lower_bit; i <= higher_bit; ++i) {
        divis.number_.push_back(divisible[i]);
    }
    cut(divis.number_);
    divis.sign_ = PLUS;
    for (size_t i = 0; i < divider.size(); ++i) {
        divid.number_.push_back(divider[i]);
    }
    cut(divid.number_);
    divid.sign_ = PLUS;

    while (right - left > 1) {
        int mid = (right + left) / 2;
        if (divid * mid > divis) {
            right = mid;
        } else {
            left = mid;
        }
    }
    return std::make_pair(left, divis - divid * left);
}

void BigInteger::modOrDiv(const BigInteger& number, bool to_mod) {
    int sign = sign_ * number.sign_;
    int module_sign = sign_;
    if (*this == 0) {
        return;
    }
    if (*this == number) {
        if (to_mod) {
            *this = 0;
            return;
        }
        *this = 1;
        return;
    }
    if (*this == -number) {
        if (to_mod) {
            *this = 0;
            return;
        }
        *this = -1;
        return;
    }
    sign_ = number.sign_;
    if (*this * sign_ < number * number.sign_) {
        if (to_mod) {
            sign_ = module_sign;
            return;
        }
        *this = 0;
        return;
    }

    BigInteger result;
    int higher_bit = -1, lower_bit = -1;
    std::pair<int, BigInteger> p;
    for (int i = static_cast<int>(number_.size() - number.number_.size());
         i >= 0;) {
        if (i == static_cast<int>(number_.size() - number.number_.size())) {
            for (size_t j = 0; j < number.number_.size(); ++j) {
                if (number_[number_.size() - 1 - j] >
                    number.number_[number.number_.size() - 1 - j]) {
                    break;
                }
                if (number_[number_.size() - 1 - j] <
                    number.number_[number.number_.size() - 1 - j]) {
                    --i;
                    break;
                }
            }
            p = DivBinsearch(number_, number.number_, i, number_.size() - 1);
            if (!to_mod) {
                result.number_.push_back(p.first);
            }
            higher_bit = lower_bit = i - 1;
            if (p.second != 0) {
                for (size_t j = i; j < i + p.second.number_.size(); ++j) {
                    number_[j] = p.second.number_[j - i];
                }
                higher_bit = i + p.second.number_.size() - 1;
            }
            --i;
        } else {
            p = DivBinsearch(number_, number.number_, lower_bit, higher_bit);
            if (!to_mod) {
                result.number_.push_back(p.first);
            }
            if (p.second != 0) {
                for (size_t j = lower_bit;
                     j < lower_bit + p.second.number_.size(); ++j) {
                    number_[j] = p.second.number_[j - lower_bit];
                }
                --lower_bit;
                i = lower_bit;
                higher_bit = lower_bit + p.second.number_.size();
            } else {
                --lower_bit;
                i = higher_bit = lower_bit;
            }
        }
    }
    if (to_mod) {
        result = p.second;
        result.sign_ = module_sign;
    } else {
        std::reverse(result.number_.begin(), result.number_.end());
        result.sign_ = sign;
    }
    *this = result;
    cut(number_);
}

BigInteger& BigInteger::operator/=(const BigInteger& number) {
    modOrDiv(number, false);
    return *this;
}

BigInteger operator/(const BigInteger& a, const BigInteger& b) {
    BigInteger tmp = a;
    tmp /= b;
    return tmp;
}

BigInteger& BigInteger::operator%=(const BigInteger& number) {
    modOrDiv(number, true);
    return *this;
}

BigInteger operator%(const BigInteger& a, const BigInteger& b) {
    BigInteger tmp = a;
    tmp %= b;
    return tmp;
}

Rational::Rational() : numerator_(0), denominator_(1) {}

Rational::Rational(const BigInteger& num, const BigInteger& den)
    : numerator_(num), denominator_(den) {
    shorten();
}

Rational::Rational(int number) : numerator_(number), denominator_(1) {}

void Rational::shorten() {
    if (numerator_.getSign() != denominator_.getSign()) {
        numerator_.setSign(-1);
        denominator_.setSign(1);
    } else {
        numerator_.setSign(1);
        denominator_.setSign(1);
    }
    BigInteger divider = gcd(numerator_, denominator_);
    numerator_ /= divider;
    denominator_ /= divider;
}

BigInteger Rational::gcd(const BigInteger& a, const BigInteger& b) {
    if (b == 0) {
        return a * a.getSign();
    }
    return gcd(b * b.getSign(), (a * a.getSign()) % (b * b.getSign()));
}

void Rational::swap(Rational& frac) {
    ::swap(numerator_, frac.numerator_);
    ::swap(denominator_, frac.denominator_);
}

Rational Rational::operator+() const {
    return *this;
}
Rational Rational::operator-() const {
    Rational tmp = *this;
    tmp.numerator_.setSign(tmp.numerator_.getSign() * (-1));
    return tmp;
}
bool operator==(const Rational& a, const Rational& b) {
    return a.getData() == b.getData();
}
bool operator<(const Rational& a, const Rational& b) {
    return a.numerator_ * b.denominator_ < b.numerator_ * a.denominator_;
}
bool operator!=(const Rational& a, const Rational& b) {
    return !(a == b);
}
bool operator>(const Rational& a, const Rational& b) {
    return b < a;
}
bool operator>=(const Rational& a, const Rational& b) {
    return !(a < b);
}
bool operator<=(const Rational& a, const Rational& b) {
    return !(a > b);
}

Rational::operator double() {
    return std::stod(asDecimal());
}

void Rational::subtrOrAdd(const Rational& frac, bool to_substr) {
    numerator_ *= frac.denominator_;
    if (to_substr) {
        numerator_ -= frac.numerator_ * denominator_;
    } else {
        numerator_ += frac.numerator_ * denominator_;
    }
    denominator_ *= frac.denominator_;
}

Rational& Rational::operator+=(const Rational& frac) {
    subtrOrAdd(frac, false);
    shorten();
    return *this;
}

Rational& Rational::operator*=(const Rational& frac) {
    numerator_ *= frac.numerator_;
    denominator_ *= frac.denominator_;
    shorten();
    return *this;
}

Rational& Rational::operator/=(const Rational& frac) {
    numerator_ *= frac.denominator_;
    denominator_ *= frac.numerator_;
    shorten();
    return *this;
}

Rational& Rational::operator-=(const Rational& frac) {
    subtrOrAdd(frac, true);
    shorten();
    return *this;
}

Rational operator+(const Rational& a, const Rational& b) {
    Rational tmp = a;
    return tmp += b;
}

Rational operator-(const Rational& a, const Rational& b) {
    Rational tmp = a;
    return tmp -= b;
}

Rational operator*(const Rational& a, const Rational& b) {
    Rational tmp = a;
    return tmp *= b;
}

Rational operator/(const Rational& a, const Rational& b) {
    Rational tmp = a;
    return tmp /= b;
}

std::string Rational::toString() const {
    std::string result;
    result += numerator_.toString();
    if (denominator_ != 1) {
        result += "/" + denominator_.toString();
    }
    return result;
}

std::string Rational::asDecimal(size_t precision) const {
    if (numerator_ == 0) {
        return "0";
    }
    if (denominator_ == 1) {
        return numerator_.toString();
    }
    int for_equals = 0;
    int sign = numerator_.getSign();
    int base = numerator_.getBase();
    BigInteger tmp_num = numerator_;
    tmp_num.setSign(1);
    BigInteger integer = tmp_num / denominator_;
    while (tmp_num < denominator_) {
        tmp_num *= 10;
        ++for_equals;
    }
    for (size_t i = 0; i < precision + 1; ++i) {
        tmp_num *= base;
    }
    tmp_num /= denominator_;
    std::string out = tmp_num.toString();
    int int_sz = integer.getData().size();
    int delta = out.size() - int_sz;
    while (static_cast<size_t>(delta) > precision) {
        out.pop_back();
        --delta;
    }
    while (static_cast<size_t>(delta) < precision) {
        out.push_back('0');
        ++delta;
    }
    std::string result;
    for (int i = 0; i < int_sz; ++i) {
        result.push_back(out[i]);
    }
    result.push_back('.');
    for (size_t i = int_sz; i < out.size(); ++i) {
        result.push_back(out[i]);
    }
    if (sign == -1) {
        result = "-" + result;
    }
    return result;
}