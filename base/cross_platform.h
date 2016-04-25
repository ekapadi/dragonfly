/**
 * @file: cross_platform.h
 *    Macro definitions for cross platform builds.
 */
 
#if !defined(__cross_platform__h)
#define __cross_platform__h

#if !defined(_WIN32)
#define __attribute_may_alias__ __attribute__((may_alias))
#else
#define __attribute_may_alias__
#endif

#endif // __cross_platform__h
