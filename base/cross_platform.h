/**
 * @file: cross_platform.h
 *    Macro definitions for cross platform builds.
 */
 
#if !defined(__cross_platform__h)
#define __cross_platform__h

// Formerly namespace name for TR1 features, now in C++11:
#define _STL_EXT_NAMESPACE_ std

#if !defined(_WIN32) && 0
#define __attribute_may_alias__ __attribute__((may_alias))
#else
#define __attribute_may_alias__
#endif

#if !defined(_WIN32) && 0
#define __attribute_unused__ __attribute__((unused))
#else
#define __attribute_unused__
#endif

#endif // __cross_platform__h
