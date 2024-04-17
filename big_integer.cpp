#include "big_integer.h"

#include <algorithm>
#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstring>
#include <iostream>
#include <ostream>
#include <stdexcept>

template <class T>
void big_integer::createBIFromNum(T a) {
  sign_ = a < 0;
  while (a != 0) {
    number_.push_back(a % BASE);
    a /= BASE;
  }
  deleteNeutral();
}

big_integer::big_integer() : sign_(false){};

big_integer::big_integer(const big_integer& other) = default;

big_integer::big_integer(int a) : big_integer(static_cast<long long>(a)){};

big_integer::big_integer(long a) : big_integer(static_cast<long long>(a)){};

big_integer::big_integer(long long a) : sign_(false) {
  createBIFromNum(a);
}

big_integer::big_integer(unsigned int a) : big_integer(static_cast<unsigned long long>(a)){};

big_integer::big_integer(unsigned long a) : big_integer(static_cast<unsigned long long>(a)){};

big_integer::big_integer(unsigned long long a) : sign_(false) {
  createBIFromNum(a);
}

big_integer::big_integer(const std::string& str) : sign_(false) {
  if (str.empty() || str == "-" || str == "+") {
    throw std::invalid_argument("Invalid number " + str);
  }
  bool negative = str[0] == '-';
  bool is_sign = negative || str[0] == '+';
  for (size_t i = is_sign; i < str.size(); i += big_integer::str_BASE_pow) {
    size_t substr_size = std::min(static_cast<size_t>(big_integer::str_BASE_pow), str.size() - i);
    *this *= static_cast<std::uint32_t>(pow(big_integer::str_BASE, substr_size));
    std::uint32_t digit = 0;
    const char* substr_end = str.data() + i + substr_size;
    std::from_chars_result res = std::from_chars(str.data() + i, substr_end, digit);
    if (res.ec == std::errc::invalid_argument || res.ptr != substr_end) {
      throw std::invalid_argument("Invalid number " + str);
    }
    if (negative) {
      *this -= digit;
    } else {
      *this += digit;
    }
  }
  deleteNeutral();
}

big_integer::~big_integer() = default;

big_integer& big_integer::operator=(const big_integer& other) = default;

big_integer& big_integer::sumOrDiff(const big_integer& rhs, bool change_sign) {
  bool overflow = false;
  if (change_sign) {
    ++*this;
  }
  size_t old_size = number_.size();
  number_.resize(std::max(number_.size(), rhs.number_.size()) + 2);
  for (size_t i = 0; i < number_.size(); ++i) {
    std::uint64_t s = static_cast<std::uint64_t>(getTrueDigit(i, old_size)) +
                      (change_sign ? ~rhs.getDigit(i) : rhs.getDigit(i)) + overflow;
    overflow = (s >= big_integer::BASE);
    number_[i] = static_cast<std::uint32_t>(s % big_integer::BASE);
  }
  sign_ = number_.back() == big_integer::BASE - 1;
  deleteNeutral();
  return *this;
}

big_integer& big_integer::operator+=(const big_integer& rhs) {
  return sumOrDiff(rhs, false);
}

big_integer& big_integer::operator+=(std::int64_t rhs) {
  bool overflow = false;
  size_t old_size = number_.size();
  number_.resize(std::max(number_.size(), sizeof(rhs) / sizeof(BASE)) + 2);
  bool b_sign = rhs < 0;
  for (size_t i = 0; i < number_.size(); ++i) {
    std::uint32_t b_dig = rhs % big_integer::BASE;
    if (b_dig == 0) {
      b_dig = big_integer::getNeutral(b_sign);
    }
    std::uint64_t s = static_cast<std::uint64_t>(getTrueDigit(i, old_size)) + b_dig + overflow;
    rhs >>= big_integer::log2_BASE;
    overflow = (s >= big_integer::BASE);
    std::uint32_t new_dig = static_cast<std::uint32_t>(s % big_integer::BASE);
    number_[i] = new_dig;
  }
  sign_ = number_.back() == big_integer::BASE - 1;
  deleteNeutral();
  return *this;
}

big_integer& big_integer::operator-=(const big_integer& rhs) {
  return sumOrDiff(rhs, true);
}

big_integer& big_integer::operator-=(const std::int64_t rhs) {
  if (rhs == std::numeric_limits<std::int64_t>::min()) {
    *this += std::numeric_limits<std::int64_t>::max();
    return ++*this;
  }
  return *this += -rhs;
}

std::vector<std::uint32_t> multiplication_cycle(const std::vector<std::uint32_t>& number1,
                                                const std::vector<std::uint32_t>& number2) {
  std::vector<std::uint32_t> res;
  res.resize(number1.size() + number2.size());
  for (size_t i = 0; i < number2.size(); ++i) {
    std::uint32_t carry = 0;
    for (size_t j = 0; j < number1.size() || carry != 0; ++j) {
      std::uint64_t tmp =
          static_cast<std::uint64_t>(number2[i]) * (j < number1.size() ? number1[j] : 0) + res[i + j] + carry;
      res[i + j] = static_cast<std::uint32_t>(tmp % big_integer::BASE);
      carry = static_cast<std::uint32_t>(tmp / big_integer::BASE);
    }
  }
  return res;
}

big_integer& big_integer::operator*=(const big_integer& rhs) {
  if (*this == 0 || rhs == 0) {
    return *this = 0;
  }
  bool new_sign = sign_ ^ rhs.sign_;
  big_integer abs_this = abs();
  big_integer abs_rhs = rhs.abs();
  sign_ = false;
  number_ = multiplication_cycle(abs_this.number_, abs_rhs.number_);
  deleteNeutral();
  if (new_sign) {
    *this = -*this;
  }
  return *this;
}

big_integer& big_integer::operator*=(const std::int64_t rhs) {
  if (*this == 0 || rhs == 0) {
    return *this = 0;
  }
  bool new_sign = sign_ ^ (rhs < 0);
  big_integer abs_this = abs();
  std::uint64_t abs_rhs;
  if (rhs == std::numeric_limits<std::int64_t>::min()) {
    abs_rhs = std::numeric_limits<std::int64_t>::max();
    ++abs_rhs;
  } else {
    abs_rhs = rhs < 0 ? -rhs : rhs;
  }
  sign_ = false;
  number_ = multiplication_cycle(abs_this.number_, {static_cast<uint32_t>(abs_rhs % big_integer::BASE),
                                                    static_cast<uint32_t>(abs_rhs / big_integer::BASE)});
  deleteNeutral();
  if (new_sign) {
    *this = -*this;
  }
  return *this;
}

std::pair<big_integer, big_integer> divWithRem(const big_integer& a, const big_integer& b) {
  if (b == 0) {
    throw std::invalid_argument("division by zero");
  }
  if (a == 0) {
    return {a, 0};
  }
  big_integer copy_a = a.abs();
  big_integer copy_b = b.abs();
  if (copy_a < copy_b) {
    return {0, a};
  }
  int shift = std::max(static_cast<int>(big_integer::log2_BASE - log2(copy_b.number_.back()) - 1), static_cast<int>(0));
  copy_a <<= shift;
  copy_b <<= shift;
  while (copy_b.number_.back() < big_integer::BASE / 2) {
    copy_a <<= 1;
    copy_b <<= 1;
    shift++;
  }
  size_t diff = copy_a.number_.size() - copy_b.number_.size();
  size_t n = copy_b.number_.size();
  big_integer q;
  q.number_.resize(diff + 1);
  copy_b <<= big_integer::log2_BASE * diff;
  if (copy_a >= copy_b) {
    q.number_[diff] = 1;
    copy_a -= copy_b;
  }
  copy_b.rightShift();
  for (auto i = diff; i > 0; i--) {
    std::uint64_t tmp =
        (static_cast<std::uint64_t>(copy_a.getDigit(n + i - 1)) * big_integer::BASE + copy_a.getDigit(n + i - 2)) /
        copy_b.number_.back();
    q.number_[i - 1] = std::min(tmp, big_integer::BASE - 1);
    big_integer moved_digit = copy_b;
    copy_a -= moved_digit * big_integer(q.number_[i - 1]);
    while (copy_a.sign_) {
      q.number_[i - 1]--;
      copy_a += moved_digit;
    }
    copy_b.rightShift();
  }
  q.deleteNeutral();
  copy_a >>= shift;
  if (a.sign_ ^ b.sign_) {
    q = -q;
  }
  if (a.sign_) {
    copy_a = -copy_a;
  }
  return {q, copy_a};
}

std::pair<big_integer, big_integer> divWithRemShort(const big_integer& a, const std::int64_t b) {
  if (b == std::numeric_limits<std::int64_t>::min() || std::abs(b) > std::numeric_limits<std::uint32_t>::max()) {
    return divWithRem(a, b);
  }
  if (b == 0) {
    throw std::invalid_argument("division by zero");
  }
  if (a == 0) {
    return {a, 0};
  }
  big_integer q = a.abs();
  std::uint32_t abs_b = std::abs(b);
  if (q.number_.size() == 1 && q.getDigit(0) < abs_b) {
    return {0, a};
  }
  std::uint32_t carry = 0;
  for (size_t i = q.number_.size(); i > 0; --i) {
    std::uint64_t tmp = q.getDigit(i - 1) + static_cast<std::uint64_t>(big_integer::BASE) * carry;
    q.number_[i - 1] = tmp / abs_b;
    carry = tmp % abs_b;
  }
  q.deleteNeutral();
  if (a.sign_ ^ (b < 0)) {
    q = -q;
  }
  return {q, carry};
}

big_integer& big_integer::operator/=(const big_integer& rhs) {
  return *this = divWithRem(*this, rhs).first;
}

big_integer& big_integer::operator/=(const std::int64_t rhs) {
  return *this = divWithRemShort(*this, rhs).first;
}

big_integer& big_integer::operator%=(const big_integer& rhs) {
  return *this = divWithRem(*this, rhs).second;
}

big_integer& big_integer::operator%=(const std::int64_t rhs) {
  return *this = divWithRemShort(*this, rhs).second;
}

template <class Operation>
big_integer& big_integer::binaryOperation(const big_integer& rhs, Operation op) {
  size_t size = std::max(number_.size(), rhs.number_.size());
  big_integer res;
  res.sign_ = op(sign_, rhs.sign_);
  res.number_.resize(size);
  for (size_t i = 0; i < size; i++) {
    res.number_[i] = op(getDigit(i), rhs.getDigit(i));
  }
  res.deleteNeutral();
  res.swap(*this);
  return *this;
}

big_integer& big_integer::operator&=(const big_integer& rhs) {
  return binaryOperation(rhs, std::bit_and<std::uint32_t>{});
}

big_integer& big_integer::operator|=(const big_integer& rhs) {
  return binaryOperation(rhs, std::bit_or<std::uint32_t>{});
}

big_integer& big_integer::operator^=(const big_integer& rhs) {
  return binaryOperation(rhs, std::bit_xor<std::uint32_t>{});
}

big_integer& big_integer::operator<<=(int rhs) {
  if (rhs < 0) {
    return *this >>= -rhs;
  }
  std::uint32_t digit_shift = rhs % big_integer::log2_BASE;
  std::uint32_t size_diff = rhs / big_integer::log2_BASE;
  number_.insert(number_.begin(), size_diff, 0);
  *this *= (static_cast<uint32_t>(1) << digit_shift);
  return *this;
}

big_integer& big_integer::operator>>=(int rhs) {
  if (rhs < 0) {
    return *this <<= -rhs;
  }
  std::uint32_t digit_shift = rhs % big_integer::log2_BASE;
  std::uint32_t size_diff = rhs / big_integer::log2_BASE;
  number_.erase(number_.begin(), number_.begin() + size_diff);
  std::pair<big_integer, big_integer> res = divWithRemShort(*this, static_cast<uint32_t>(1) << digit_shift);
  swap(res.first);
  if (sign_ && rhs != 0 && res.second != 0) {
    --*this;
  }
  return *this;
}

big_integer& big_integer::rightShift() {
  if (number_.empty()) {
    return *this;
  }
  for (size_t i = 0; i < number_.size() - 1; i++) {
    number_[i] = number_[i + 1];
  }
  number_.pop_back();
  return *this;
}

big_integer big_integer::operator+() const {
  return *this;
}

big_integer big_integer::operator-() const {
  return (~*this) += 1;
}

big_integer big_integer::operator~() const {
  big_integer res;
  res.number_.resize(number_.size());
  for (size_t i = 0; i < res.number_.size(); i++) {
    res.number_[i] = ~number_[i];
  }
  res.sign_ = !sign_;
  res.deleteNeutral();
  return res;
}

big_integer& big_integer::operator++() {
  return *this += 1;
}

big_integer big_integer::operator++(int) {
  big_integer answer = *this;
  ++*this;
  return answer;
}

big_integer& big_integer::operator--() {
  return *this -= 1;
}

big_integer big_integer::operator--(int) {
  *this -= 1;
  return *this + 1;
}

big_integer operator+(const big_integer& a, const big_integer& b) {
  return big_integer(a) += b;
}

big_integer operator+(const big_integer& a, std::int64_t b) {
  return big_integer(a) += b;
}

big_integer operator+(std::int64_t b, const big_integer& a) {
  return big_integer(a) += b;
}

big_integer operator-(const big_integer& a, const big_integer& b) {
  return big_integer(a) -= b;
}

big_integer operator-(const big_integer& a, const std::int64_t b) {
  return big_integer(a) -= b;
}

big_integer operator*(const big_integer& a, const big_integer& b) {
  return big_integer(a) *= b;
}

big_integer operator*(const big_integer& a, const std::int64_t b) {
  return big_integer(a) *= b;
}

big_integer operator*(const std::int64_t b, const big_integer& a) {
  return big_integer(a) *= b;
}

big_integer operator/(const big_integer& a, const big_integer& b) {
  return divWithRem(a, b).first;
}

big_integer operator/(const big_integer& a, const std::int64_t b) {
  return divWithRemShort(a, b).first;
}

big_integer operator%(const big_integer& a, const big_integer& b) {
  return divWithRem(a, b).second;
}

big_integer operator%(const big_integer& a, const std::int64_t b) {
  return divWithRemShort(a, b).second;
}

big_integer operator&(const big_integer& a, const big_integer& b) {
  return big_integer(a) &= b;
}

big_integer operator|(const big_integer& a, const big_integer& b) {
  return big_integer(a) |= b;
}

big_integer operator^(const big_integer& a, const big_integer& b) {
  return big_integer(a) ^= b;
}

big_integer operator<<(const big_integer& a, int b) {
  return big_integer(a) <<= b;
}

big_integer operator>>(const big_integer& a, int b) {
  return big_integer(a) >>= b;
}

bool operator==(const big_integer& a, const big_integer& b) {
  if (a.sign_ != b.sign_ || a.number_.size() != b.number_.size()) {
    return false;
  }
  return std::equal(a.number_.begin(), a.number_.end(), b.number_.begin(), b.number_.end());
}

bool operator!=(const big_integer& a, const big_integer& b) {
  return !(a == b);
}

bool operator<(const big_integer& a, const big_integer& b) {
  if (a == b) {
    return false;
  }
  if (a.sign_ && !b.sign_) {
    return true;
  }
  if (b.sign_ && !a.sign_) {
    return false;
  }
  if (a.number_.size() != b.number_.size()) {
    return (a.number_.size() < b.number_.size()) ^ a.sign_;
  }
  return std::lexicographical_compare(a.number_.rbegin(), a.number_.rend(), b.number_.rbegin(), b.number_.rend());
}

bool operator>(const big_integer& a, const big_integer& b) {
  return b < a;
}

bool operator<=(const big_integer& a, const big_integer& b) {
  return !(a > b);
}

bool operator>=(const big_integer& a, const big_integer& b) {
  return !(a < b);
}

std::string to_string(const big_integer& a) {
  if (a.number_.empty() && !a.sign_) {
    return "0";
  }
  std::string res;
  std::pair<big_integer, big_integer> div = {a.abs(), 0};
  const int32_t str_base = pow(big_integer::str_BASE, big_integer::str_BASE_pow);

  while (div.first != 0) {
    div = divWithRemShort(div.first, str_base);
    std::uint32_t digit = div.second.getDigit(0);
    for (size_t i = 0; digit != 0 || (i < big_integer::str_BASE_pow && div.first != 0); ++i) {
      res += digit % big_integer::str_BASE + '0';
      digit /= big_integer::str_BASE;
    }
  }

  if (a.sign_) {
    res += '-';
  }
  std::reverse(res.begin(), res.end());
  return res;
}

std::ostream& operator<<(std::ostream& out, const big_integer& a) {
  return out << to_string(a);
}

void big_integer::swap(big_integer& other) {
  std::swap(number_, other.number_);
  std::swap(sign_, other.sign_);
}

std::uint32_t big_integer::getNeutral(bool sign) {
  return sign ? big_integer::BASE - 1 : 0;
}

std::uint32_t big_integer::getDigit(size_t index) const {
  return getTrueDigit(index, number_.size());
}

std::uint32_t big_integer::getTrueDigit(size_t index, size_t true_size) const {
  return true_size > index ? number_[index] : getNeutral(sign_);
}

void big_integer::deleteNeutral() {
  while (!number_.empty() && number_.back() == getNeutral(sign_)) {
    number_.pop_back();
  }
}

big_integer big_integer::abs() const {
  return *this < 0 ? -*this : *this;
}
