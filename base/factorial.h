#if !defined(__FACTORIAL__H)
#define __FACTORIAL__H
/**
 * @file: factorial.h
 *   Templated, lookup-table based factorial functions.
 */
 
namespace linalg{


template <class T>
T factorial(size_t n);

// ------------------ Specializations: --------------------
template <>
double factorial<double>(size_t n);

template <>
long factorial<long>(size_t n);
// --------------------------------------------------------


} // namespace linalg
 
#endif // __FACTORIAL__H
