/*
 * wrap_traits.hpp
 */

#ifndef WRAP_TRAITS_HPP_
#define WRAP_TRAITS_HPP_

#include <cstdint>

template <std::size_t S>
struct WrapTraits {
    static uint8_t SHIFT_TO_OFFSET_IDX[][127];
};

#endif /* WRAP_TRAITS_HPP_ */
