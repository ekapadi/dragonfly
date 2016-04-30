/**
 * @file: Factorial_functor.h
 *   Double-precision *only*, lookup-table based factorial function.
 */
#include <assert.h>
#include <complex>
#include <array>
#include <mutex>
#include <functional>

#include "numberTraits.h"
#include "factorial.h"

namespace{


// -------------- Declarations: ----------------------------------------------
const size_t _N_MAX_double = 170;
std::array<double, _N_MAX_double + 1> _table_double = {};
std::once_flag _table_once_double;

const size_t _N_MAX_long = 20;
std::array<long, _N_MAX_long + 1> _table_long = {};
std::once_flag _table_once_long;

template <class T, size_t N_MAX>
void _init_table(std::array<T, N_MAX + 1> &table);

// "call_once" wrappers:
void _init_table_double(std::array<double, _N_MAX_double + 1> &table);

void _init_table_long(std::array<long, _N_MAX_long + 1> &table);

// ---------------------------------------------------------------------------


template <class T, size_t N_MAX>
void _init_table(std::array<T, N_MAX + 1> &table)
{
    table[0] = number::one<T>();
    for(size_t n = 1; n <= N_MAX; ++n)
      table[n] = number::integer<T>(n) * table[n-1];
}

void _init_table_double(std::array<double, _N_MAX_double + 1> &table)
{ _init_table<double, _N_MAX_double>(table); }

void _init_table_long(std::array<long, _N_MAX_long + 1> &table)
{ _init_table<long, _N_MAX_long>(table); }


} // namespace

 
namespace linalg{


// Double-precision, lookup-table based Factorial functor.
template <>
double factorial<double>(size_t n)
{
  std::call_once(_table_once_double, _init_table_double, std::ref(_table_double));

  assert(n <= _N_MAX_double);
  return _table_double[n];
}


// 64-bit, lookup-table based Factorial functor.
template <>
long factorial<long>(size_t n)
{
  std::call_once(_table_once_long, _init_table_long, std::ref(_table_long));

  assert(n <= _N_MAX_long);
  return _table_long[n];
}


} // namespace linalg
