#if !defined(__HASH_COMBINE__H)
#define __HASH_COMBINE__H

namespace linalg{


// Combine hash values for any type with a defined "std::hash".
// (This definition aligns with that of the "boost" implementation, and will probably be *eventually* added to the C++ standard.)
template <class T, class HASH>
inline void hash_combine(std::size_t & seed, const T & v)
{
  HASH hasher;
  seed ^= hasher(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
}


} // namespace linalg

#endif // __HASH_COMBINE__H
