#pragma once

#include <iosfwd>
#include <limits>
#include <string>
#include <vector>

struct big_integer {

private:
  static const std::uint64_t BASE = static_cast<std::uint64_t>(std::numeric_limits<std::uint32_t>::max()) + 1;
  static const std::uint32_t log2_BASE = 32;
  static const std::uint32_t str_BASE = 10;
  static const std::uint32_t str_BASE_pow = 9;
  std::vector<std::uint32_t> number_;
  bool sign_;
  template <class Operation>
  big_integer& binaryOperation(const big_integer& rhs, Operation op);
  static std::uint32_t getNeutral(bool sign);
  std::uint32_t getDigit(size_t index) const;
  std::uint32_t getTrueDigit(size_t index, size_t true_size) const;
  void deleteNeutral();
  big_integer& rightShift();
  friend std::pair<big_integer, big_integer> divWithRem(const big_integer& copy_a, const big_integer& copy_b);
  friend std::pair<big_integer, big_integer> divWithRemShort(const big_integer& a, std::int64_t b);
  friend big_integer shift(const big_integer& a, int b, bool is_left);
  friend std::vector<std::uint32_t> multiplication_cycle(const std::vector<std::uint32_t>& number1,
                                                         const std::vector<std::uint32_t>& number2);
  template <class T>
  void createBIFromNum(T a);
  big_integer& sumOrDiff(const big_integer& rhs, bool change_sign);

public:
  big_integer();
  big_integer(const big_integer& other);
  big_integer(int a);
  big_integer(long a);
  big_integer(long long a);
  big_integer(unsigned int a);
  big_integer(unsigned long a);
  big_integer(unsigned long long a);

  explicit big_integer(const std::string& str);
  ~big_integer();

  big_integer abs() const;
  big_integer& operator=(const big_integer& other);

  big_integer& operator+=(const big_integer& rhs);
  big_integer& operator+=(std::int64_t rhs);
  big_integer& operator-=(const big_integer& rhs);
  big_integer& operator-=(std::int64_t rhs);
  big_integer& operator*=(const big_integer& rhs);
  big_integer& operator*=(std::int64_t rhs);
  big_integer& operator/=(const big_integer& rhs);
  big_integer& operator/=(std::int64_t rhs);
  big_integer& operator%=(const big_integer& rhs);
  big_integer& operator%=(std::int64_t rhs);

  big_integer& operator&=(const big_integer& rhs);
  big_integer& operator|=(const big_integer& rhs);
  big_integer& operator^=(const big_integer& rhs);

  big_integer& operator<<=(int rhs);
  big_integer& operator>>=(int rhs);

  big_integer operator+() const;
  big_integer operator-() const;
  big_integer operator~() const;

  big_integer& operator++();
  big_integer operator++(int);

  big_integer& operator--();
  big_integer operator--(int);

  friend bool operator==(const big_integer& a, const big_integer& b);
  friend bool operator!=(const big_integer& a, const big_integer& b);
  friend bool operator<(const big_integer& a, const big_integer& b);
  friend bool operator>(const big_integer& a, const big_integer& b);
  friend bool operator<=(const big_integer& a, const big_integer& b);
  friend bool operator>=(const big_integer& a, const big_integer& b);

  friend big_integer operator%(const big_integer& a, const big_integer& b);
  friend big_integer operator%(const big_integer& a, std::int64_t b);

  friend big_integer operator<<(const big_integer& a, int b);
  friend big_integer operator>>(const big_integer& a, int b);

  friend std::string to_string(const big_integer& a);
  void swap(big_integer& integer);
};

big_integer operator+(const big_integer& a, const big_integer& b);
big_integer operator+(const big_integer& a, std::int64_t b);
big_integer operator+(std::int64_t b, const big_integer& a);
big_integer operator-(const big_integer& a, const big_integer& b);
big_integer operator-(const big_integer& a, std::int64_t b);
big_integer operator*(const big_integer& a, const big_integer& b);
big_integer operator*(const big_integer& a, std::int64_t b);
big_integer operator*(std::int64_t b, const big_integer& a);
big_integer operator/(const big_integer& a, const big_integer& b);
big_integer operator/(const big_integer& a, std::int64_t b);
big_integer operator&(const big_integer& a, const big_integer& b);
big_integer operator|(const big_integer& a, const big_integer& b);
big_integer operator^(const big_integer& a, const big_integer& b);

bool operator==(const big_integer& a, const big_integer& b);
bool operator!=(const big_integer& a, const big_integer& b);
bool operator<(const big_integer& a, const big_integer& b);
bool operator>(const big_integer& a, const big_integer& b);
bool operator<=(const big_integer& a, const big_integer& b);
bool operator>=(const big_integer& a, const big_integer& b);

std::string to_string(const big_integer& a);
std::ostream& operator<<(std::ostream& out, const big_integer& a);
