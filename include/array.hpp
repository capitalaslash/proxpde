#pragma once

// array class with improved constexpr usability.
// based on sprout: https://github.com/bolero-MURAKAMI/Sprout

#include <cstddef>
#include <array>

// #define SPROUT_USE_INDEX_ITERATOR_IMPLEMENTATION 1

using index_t = std::ptrdiff_t;

template<index_t... Indexes>
using index_tuple = std::integer_sequence<index_t, Indexes...>;

template<index_t N>
using make_index_tuple = std::make_integer_sequence<index_t, N>;

template<typename T, std::size_t N>
class array {
public:
  typedef T value_type;
#if SPROUT_USE_INDEX_ITERATOR_IMPLEMENTATION
  typedef sprout::index_iterator<array&, true> iterator;
  typedef sprout::index_iterator<array const&, true> const_iterator;
#else
  typedef T* iterator;
  typedef T const* const_iterator;
#endif
  typedef T& reference;
  typedef T const& const_reference;
  typedef std::size_t size_type;
  typedef std::ptrdiff_t difference_type;
  typedef T* pointer;
  typedef T const* const_pointer;
  typedef std::reverse_iterator<iterator> reverse_iterator;
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
public:
  static constexpr size_type static_size = N;
private:
  template<index_t Index>
  static constexpr value_type const&
  dummy_get(value_type const& value) {
    return value;
  }
  template<index_t... Indexes>
  static constexpr array
  fill_impl(value_type const& value, index_tuple<Indexes...>) {
    return array{{dummy_get<Indexes>(value)...}};
  }
public:
  value_type elems[static_size ? static_size : 1];
private:
  template<index_t... Indexes>
  constexpr std::array<T, N>
  to_std_array(index_tuple<Indexes...>) const
  noexcept
  {
    return std::array<T, N>{{elems[Indexes]...}};
  }
public:
  // construct/copy/destroy:
  template<typename T2>
  constexpr array& operator=(array<T2, N> const& rhs) {
    std::copy(rhs.begin(), rhs.end(), begin());
    return *this;
  }
  template<typename T2>
  constexpr array& operator=(array<T2, N>&& rhs) {
    std::move(rhs.begin(), rhs.end(), begin());
    return *this;
  }
  // modifiers:
  constexpr void fill(const_reference value) {
    std::fill_n(begin(), size(), value);
  }
  // constexpr array fill(const_reference value) const {
  // 	return fill_impl(value, sprout::index_n<0, N>::make());
  // }
  constexpr void assign(const_reference value) {
    fill(value);
  }
  constexpr array assign(const_reference value) const {
    return fill(value);
  }
  constexpr void swap(array& other)
  noexcept
  {
    std::swap_ranges(other.begin(), other.end(), begin());
  }
  // iterators:
#if SPROUT_USE_INDEX_ITERATOR_IMPLEMENTATION
  constexpr iterator begin() noexcept {
    return iterator(*this, 0);
  }
  constexpr const_iterator begin() const noexcept {
    return const_iterator(*this, 0);
  }
  constexpr iterator end() noexcept {
    return iterator(*this, size());
  }
  constexpr const_iterator end() const noexcept {
    return const_iterator(*this, size());
  }
#else
  constexpr iterator begin() noexcept {
    return iterator(elems);
  }
  constexpr const_iterator begin() const noexcept {
    return iterator(elems);
  }
  constexpr iterator end() noexcept {
    return iterator(elems) + size();
  }
  constexpr const_iterator end() const noexcept {
    return iterator(elems) + size();
  }
#endif
  constexpr reverse_iterator rbegin() noexcept {
    return reverse_iterator(end());
  }
  constexpr const_reverse_iterator rbegin() const noexcept {
    return const_reverse_iterator(end());
  }
  constexpr reverse_iterator rend() noexcept {
    return reverse_iterator(begin());
  }
  constexpr const_reverse_iterator rend() const noexcept {
    return const_reverse_iterator(begin());
  }
#if SPROUT_USE_INDEX_ITERATOR_IMPLEMENTATION
  constexpr const_iterator cbegin() const noexcept {
    return const_iterator(*this, 0);
  }
  constexpr const_iterator cend() const noexcept {
    return const_iterator(*this, size());
  }
#else
  constexpr const_iterator cbegin() const noexcept {
    return const_iterator(elems);
  }
  constexpr const_iterator cend() const noexcept {
    return const_iterator(elems) + size();
  }
#endif
  constexpr const_reverse_iterator crbegin() const noexcept {
    return const_reverse_iterator(end());
  }
  constexpr const_reverse_iterator crend() const noexcept {
    return const_reverse_iterator(begin());
  }
  // capacity:
  constexpr size_type size() const noexcept {
    return static_size;
  }
  constexpr size_type max_size() const noexcept {
    return size();
  }
  constexpr bool empty() const noexcept {
    return size() == 0;
  }
  // element access:
  constexpr reference operator[](size_type i) {
    return elems[i];
  }
  constexpr const_reference operator[](size_type i) const {
    return elems[i];
  }
  constexpr reference at(size_type i) {
    return i < size() ? elems[i]
                        : (throw std::out_of_range("array<>: index out of range"), elems[i])
                        ;
  }
  constexpr const_reference at(size_type i) const {
    return i < size() ? elems[i]
                        : (throw std::out_of_range("array<>: index out of range"), elems[i])
                        ;
  }
  constexpr reference front() {
    return elems[0];
  }
  constexpr const_reference front() const {
    return elems[0];
  }
  constexpr reference back() {
    return elems[size() - 1];
  }
  constexpr const_reference back() const {
    return elems[size() - 1];
  }

  constexpr pointer data() noexcept {
    return elems;
  }
  constexpr const_pointer data() const noexcept {
    return elems;
  }
  constexpr pointer c_array() noexcept {
    return data();
  }
  constexpr const_pointer c_array() const noexcept {
    return data();
  }
  // others:
  constexpr void rangecheck(size_type i) const {
    return i >= size() ? throw std::out_of_range("array<>: index out of range")
                       : (void)0
                         ;
  }

#if SPROUT_USE_INDEX_ITERATOR_IMPLEMENTATION
  constexpr iterator nth(size_type i) {
    return i < size() ? iterator(*this, i)
                      : (throw std::out_of_range("array<>: index out of range"), iterator())
                        ;
  }
  constexpr const_iterator nth(size_type i) const {
    return i < size() ? const_iterator(*this, i)
                      : (throw std::out_of_range("array<>: index out of range"), const_iterator())
                        ;
  }
  constexpr size_type index_of(iterator p) noexcept {
    return p.index();
  }
  constexpr size_type index_of(const_iterator p) const noexcept {
    return p.index();
  }
#else
  constexpr iterator nth(size_type i) {
    return i < size() ? iterator(elems) + i
                      : (throw std::out_of_range("array<>: index out of range"), iterator())
                        ;
  }
  constexpr const_iterator nth(size_type i) const {
    return i < size() ? const_iterator(elems) + i
                      : (throw std::out_of_range("array<>: index out of range"), const_iterator())
                        ;
  }
  constexpr size_type index_of(iterator p) noexcept {
    return std::distance(begin(), p);
  }
  constexpr size_type index_of(const_iterator p) const noexcept {
    return std::distance(begin(), p);
  }
#endif

  constexpr operator std::array<T, N>() const
  noexcept
  {
    return to_std_array(make_index_tuple<N>::make());
  }
};
template<typename T, std::size_t N>
constexpr typename array<T, N>::size_type array<T, N>::static_size;

//
// swap
//
template<typename T, std::size_t N>
inline constexpr void
swap(array<T, N>& lhs, array<T, N>& rhs)
noexcept
{
  lhs.swap(rhs);
}

//
// to_array
//
template<typename T, std::size_t N>
inline constexpr array<T, N>
to_array(array<T, N> const& arr)
noexcept
{
  return arr;
}
namespace detail {
template<typename T, std::size_t N, index_t... Indexes>
inline constexpr array<std::remove_cv_t<T>, N>
to_array_impl(T (& arr)[N], index_tuple<Indexes...>)
noexcept
{
  return array<std::remove_cv_t<T>, N>{{arr[Indexes]...}};
}
}	// namespace detail
template<typename T, std::size_t N>
inline constexpr array<std::remove_cv_t<T>, N>
to_array(T (& arr)[N])
noexcept
{
  return detail::to_array_impl(arr, make_index_tuple<N>::make());
}
namespace detail {
template<typename T, std::size_t N, index_t... Indexes>
inline constexpr array<T, N>
to_array_impl(std::array<T, N> const& arr, index_tuple<Indexes...>)
noexcept
{
  return array<T, N>{{arr[Indexes]...}};
}
}	// namespace detail
template<typename T, std::size_t N>
inline constexpr array<T, N>
to_array(std::array<T, N> const& arr)
noexcept
{
  return detail::to_array_impl(arr, make_index_tuple<N>::make());
}
