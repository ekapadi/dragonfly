#if !defined(__parallelUtil__h)
#define __parallelUtil__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/parallelUtil.h,v $

/* *********************************************************************** */
/*                                                                         */
/* Copyright (C) 2010  Kort Travis                                         */
/*                                                                         */
/*                                                                         */
/* This program is free software; you can redistribute it and/or modify    */
/* it under the terms of the GNU Lesser General Public License as          */
/* published by the Free Software Foundation; version 2.1 of the License.  */
/*                                                                         */
/* This program is distributed in the hope that it will be useful,         */
/* but WITHOUT ANY WARRANTY; without even the implied warranty of          */
/* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           */
/* GNU Lesser General Public License for more details.                     */
/*                                                                         */
/* You should have received a copy of the GNU Lesser General Public        */
/* License along with this program; if not, write to the Free Software     */
/* Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307,  */
/* USA.                                                                    */
/*                                                                         */
/* *********************************************************************** */


// requires:
// <pthread.h>
// <assert.h>
// <vector>
// <list>
// <algorithm>

/** @file parallelUtil.h
@brief Support classes and methods for object-oriented multi-thread and multi-process coding.
*/

// doxygen comment duplicated here, just in case this file is scanned first:
/**
 * @addtogroup parallel parallel programming
 * @brief Support classes and methods for object-oriented multi-thread and multi-process coding.
 */

/**
 * @namespace parallelUtil @ingroup parallel
 * @brief Support classes and methods for object-oriented multi-thread and multi-process coding.
 
 *  Emphasis is on providing a minimal set of efficient methods, with an interface as unified as possible.
 *  Methods are provided to supercede standard OMP pragmas, and to assist in making MPI utilization completely transparent.
 */
  
namespace parallelUtil{

// OpenMP loops require integer loop variable.
// The following facilitates the replacement of loops with multiple iterator indices.

// example of type multiple iterator loop (from "tmatrix::initLegendreTransformTables"):
#if 0
  // Sub-interval accumulation
  for( typename std::vector<T>::const_iterator itIntervalStart = vSubintervals.begin(), itIntervalEnd = ++vSubintervals.begin(),
         itSubintervalsEnd = vSubintervals.end();

       itIntervalEnd != itSubintervalsEnd;

       ++itIntervalStart,
       ++itIntervalEnd){ 
#endif

// Initial implementation: copy to array of iterators 
//  (works best for iterators on classes like linked lists and hash tables, which do not necessarily have
//   an efficient iterator difference operation.)
// This does not work efficiently when there is an extremely large number of indices.
// Note: iterator may be an integral type, such as "long".

// Allow through initialization, but do not _require_ the two lists to be synchronized, with the same number of elements.
// (For example, if iterators are removed from list1 using "erase", but not from list2, using "erase2").

template <class IT1, class IT2, class U1, class U2>
class loopIteratorList{
  protected:
		
	  U1 uIT1_;
		U2 uIT2_;
    
		void copy(const loopIteratorList& other);

  public:

  typedef typename U1::iterator iterator;
	typedef typename U1::const_iterator const_iterator;
  typedef typename U1::reverse_iterator reverse_iterator;
	typedef typename U1::const_reverse_iterator const_reverse_iterator;	
  typedef typename U2::iterator iterator2;
	typedef typename U2::const_iterator const_iterator2;
  typedef typename U2::reverse_iterator reverse_iterator2;
	typedef typename U2::const_reverse_iterator const_reverse_iterator2;	

  inline iterator begin(void) { return uIT1_.begin(); }
	inline const_iterator begin(void)const { return uIT1_.begin(); }
  inline reverse_iterator rbegin(void) { return uIT1_.rbegin(); }
	inline const_reverse_iterator rbegin(void)const { return uIT1_.rbegin(); }
  inline iterator end(void) { return uIT1_.end(); }
	inline const_iterator end(void)const { return uIT1_.end(); }
  inline reverse_iterator rend(void) { return uIT1_.rend(); }
	inline const_reverse_iterator rend(void)const { return uIT1_.rend(); }

  inline iterator erase(iterator pos) { return uIT1_.erase(pos); }
  inline bool empty(void) const { return uIT1_.empty(); }
  inline size_t size(void) const { return uIT1_.size(); }

	inline void clear(void) 
	{ 
	  uIT1_.clear(); 
		uIT2_.clear(); 
	}
	
  inline iterator2 begin2(void) { return uIT2_.begin(); }
	inline const_iterator2 begin2(void)const { return uIT2_.begin(); }
  inline reverse_iterator2 rbegin2(void) { return uIT2_.rbegin(); }
	inline const_reverse_iterator2 rbegin2(void)const { return uIT2_.rbegin(); }
  inline iterator2 end2(void) { return uIT2_.end(); }
	inline const_iterator2 end2(void)const { return uIT2_.end(); }
  inline reverse_iterator2 rend2(void) { return uIT2_.rend(); }
	inline const_reverse_iterator2 rend2(void)const { return uIT2_.rend(); }
	
  inline iterator erase2(iterator2 pos) { return uIT2_.erase(pos); }
  inline bool empty2(void) const { return uIT2_.empty(); }
  inline size_t size2(void) const { return uIT2_.size(); }

  inline const IT1& operator[](size_t n)const
  { return *(begin() + n); }

  inline loopIteratorList& operator=(const loopIteratorList& other)
	{
	 copy(other);
	 return *this;
	}
	
	// Create iterator range subregion; makeup (if any) occurs in last split.
	//   subregion can be either banded (default) or distributed modulo N_splits (comb==true).
	loopIteratorList split(size_t n_split, size_t N_splits, bool comb=false)const;

  void split(size_t n_split, size_t N_splits, bool comb, loopIteratorList& val)const;

#if 1	
  // retain derived class type:
  inline loopIteratorList* cloneSplit(size_t n_split, size_t N_splits, bool comb=false)const;
#endif
	
	virtual loopIteratorList* clone(void)const;

  void write(std::ostream& os)const;
		
	virtual ~loopIteratorList( void );

  loopIteratorList( const IT1& it1Start, const IT1& it1End, const IT2& it2Start );
	
	// one iterator set only	
  loopIteratorList( const IT1& it1Start, const IT1& it1End );

	// copy constructor:
	loopIteratorList( const loopIteratorList& other );

  // default constructor:		
	loopIteratorList( void );

}; // class loopIteratorList


template <class IT1, class IT2=IT1>
class loopIteratorArray: public loopIteratorList<IT1,IT2,std::vector<IT1>,std::vector<IT2> >{
  protected:
	    
	void copy(const loopIteratorArray& other);

  public:

  typedef std::vector<IT1> U1;
	typedef std::vector<IT2> U2;
  typedef loopIteratorList<IT1,IT2,U1,U2> base_class;

	inline const IT1& it(int n)const
	{ return *(base_class::uIT1_.begin() + n); }
	
	inline const IT2& it2(int n)const
	{ return *(base_class::uIT2_.begin() + n); }

  inline const IT1& operator[](size_t n)const
  { return it(n); }
  
		
  inline loopIteratorArray& operator=(const loopIteratorArray& other)
	{
	 copy(other);
	 return *this;
	}

	
	virtual base_class* clone(void)const;
	
	virtual ~loopIteratorArray( void );	

  loopIteratorArray( const IT1& it1Start, const IT1& it1End, const IT2& it2Start );
	
	// one iterator set only	
  loopIteratorArray( const IT1& it1Start, const IT1& it1End );

	// copy constructor:
	loopIteratorArray( const loopIteratorArray& other );

  // default constructor:		
	loopIteratorArray( void );

}; // class loopIteratorArray



// The following specialized mutex class is applicable to protecting
//   large, complicated structures from race-conditions.  
// It is assumed that threads are unlikely to access the same part of the structure at the same time, 
//   and so each thread is allowed to lock only the portion of the structure that it is working on.
// "Portions" are compared using their pointer addresses.
// Example: each thread is allowed to _write_ to a distinct quadrature cache node at the same time.
// Usage: the blocking "size_t nLock = lock(void* key)" method returns a size_t index which must be used in the corresponding "void unlock(nLock)" method.

class multiMutex{

  pthread_mutex_t controlMutex_;

	struct namedLock{
	  namedLock(void);
			
    // unsigned char __PAD0__[sizeof(pthread_mutex_t)]; // *** DEBUG *** destruction order problem. 
		pthread_mutex_t mutex_;
    pthread_t tid_;
    unsigned long key_;
		size_t refCount_;
    // unsigned char __PAD1__[sizeof(pthread_mutex_t)]; // *** DEBUG *** destruction order problem.
#if defined(_GLIBCXX_IOSTREAM) || (defined(__PGI) && defined(_STLP_IOSTREAM))
    void write(std::ostream& os)const;  
#endif
	};

	size_t N_locks_;
  namedLock *vLock_;
		
public:

  void reset(void);

		
	size_t lock(unsigned long key);

  /*
   * returns true if lock has been acquired
	 */
  bool trylock(unsigned long key, size_t& nLock);


  template <class K>
  size_t lock(const K& k);


  /*
   * returns true if lock has been acquired
	 */
  template <class K>
  bool trylock(const K& k, size_t& nLock);

  
	void unlock(size_t nLock);

  void allocLocks(size_t N_locks);
	
  size_t N_locks(void) { return N_locks_; }

#if defined(_GLIBCXX_IOSTREAM) || (defined(__PGI) && defined(_STLP_IOSTREAM))
    void write(std::ostream& os)const;
#endif	

  multiMutex(size_t N_locks=1);

	~multiMutex(void);
			
};

// some constants to initialize mutexes:
#if defined(__PGI)
  // TACC port: (these are GNU specific, but generally useful, see note in "portability.h"):
  # define PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP \
    {0, 0, 0, PTHREAD_MUTEX_RECURSIVE_NP, __LOCK_INITIALIZER}
  # define PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP \
    {0, 0, 0, PTHREAD_MUTEX_ERRORCHECK_NP, __LOCK_INITIALIZER}
  # define PTHREAD_ADAPTIVE_MUTEX_INITIALIZER_NP \
    {0, 0, 0, PTHREAD_MUTEX_ADAPTIVE_NP, __LOCK_INITIALIZER}
#endif

//! @cond full_docs
const pthread_mutex_t CONST_PTHREAD_MUTEX_INITIALIZER = PTHREAD_MUTEX_INITIALIZER;
const pthread_mutex_t CONST_PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP = PTHREAD_RECURSIVE_MUTEX_INITIALIZER_NP;
const pthread_mutex_t CONST_PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP = PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP;
// and conditionals:
const pthread_cond_t CONST_PTHREAD_COND_INITIALIZER = PTHREAD_COND_INITIALIZER;
///< @endcond


/**
 * @brief A mapping of "#pragma omp parallel for schedule(static,1) default(shared)" onto the PTHREADS API.
 *
 * @param[in] threadFunc  pointer to global or static member "thread function": 
 *    signature: void (*threadFunc)(WS&, const loopIteratorArray<IT1,IT2>& indices)
 * @param[in,out]  workspace  workspace structure containing any shared variables
 * @param[in]      indices  list of index iterators, one iterator for each thread tasks 
 */
template <class WS, class IT1, class IT2, class U1, class U2>
bool OMP_parallel_for( void (*threadFunc)(WS&, const loopIteratorList<IT1,IT2,U1,U2>&),
                       WS& workspace,                        
											 const loopIteratorList<IT1,IT2,U1,U2>& indices );

#if 0
// this version helps prevent compiler-confusion:
template <class WS, class IT1, class IT2>
inline bool OMP_parallel_for( void (*threadFunc)(WS&, const loopIteratorArray<IT1,IT2>&),
                              WS& workspace,                        
											        const loopIteratorArray<IT1,IT2>& indices );
#endif
											 
template <class WS, class IT1, class IT2, class U1, class U2>
void* OMP_parallel_for_aux_const(void *pvArgs);

// for convenience, allow a [non-OMP corresponding] variant that can modify its iterator array:

/**
 * @brief A mapping of "#pragma omp parallel for schedule(static,1) default(shared)" onto the PTHREADS API: special [non-OMP corresponding] variant allows mutable iterator list.
 *
 * @param[in] threadFunc  pointer to global or static member "thread function": 
 *    signature: void (*threadFunc)(WS&, const loopIteratorArray<IT1,IT2>& indices)
 * @param[in,out]  workspace  workspace structure containing any shared variables
 * @param[in,out]      indices  list of index iterators, one iterator for each thread tasks 
 */
template <class WS, class IT1, class IT2, class U1, class U2>
bool OMP_parallel_for( void (*threadFunc)(WS&, loopIteratorList<IT1,IT2,U1,U2>&),
                       WS& workspace,                        
											 loopIteratorList<IT1,IT2,U1,U2>& indices );

template <class WS, class IT1, class IT2, class U1, class U2>
void* OMP_parallel_for_aux_non_const(void *pvArgs);


// get the environment variable: OMP_NUM_THREADS (if it doesn't exist, return 1):
size_t OMP_NUM_THREADS(void);

// ----------------------------------------------------------------------------------------------------------------------------------------

// analogous parallel loop implementation using member-functions:
//   MEMBER_FUNC is of type
//     RVAL (CLASS_::*MEMBER_FUNC)(<typedef iterator_traits<IT>::reference OR " "::value_type>)[<optional "const" spec>][<optional "throw" spec>]
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT>
std::vector<RVAL> parallel_for(
  CLASS_ *instance, 
  MEMBER_FUNC threadFunc, 
  IT begin_, IT end_, size_t N_THREAD = OMP_NUM_THREADS());
  
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT>
void* parallel_for_aux(void *pvArgs);

// +1 ARG version:
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT, class ARG>
std::vector<RVAL> parallel_for(
  CLASS_ *instance, 
  MEMBER_FUNC threadFunc, 
  IT begin_, IT end_, ARG arg, size_t N_THREAD = OMP_NUM_THREADS());
  
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT, class ARG>
void* parallel_for_aux(void *pvArgs);

// ----------------------------------------------------------------------------

//! @cond full_docs
// utility classes for pthread API methods "void*" parameter passing and "void*" return-value implementation under strict-aliasing rules:
template <class S>
union void_ptr_union{
  S *ptr;
  void *void_ptr;
};
///< @endcond

// return true if the calling thread is the process "root" thread:
bool is_process_root(void);

// ----------------------------------------------------------------------------------------------------------------------------------------

/**
 @cond full_docs 
 @brief manipulators to lock an iostream associated mutex.
 */
 
struct _iostream_thread_lock_base{
#if defined(__USE_PTHREAD)  
  static pthread_mutex_t mutex_;
#endif
};

struct _iostream_thread_lock: private _iostream_thread_lock_base {
  static inline void lock(void)
  { 
#if defined(__USE_PTHREAD)   
    if (pthread_mutex_lock(&mutex_))
      throw std::runtime_error("_iostream_thread_lock::lock: pthread_mutex_lock error return");
#endif
  }
};
extern _iostream_thread_lock log_lock; 
 
struct _iostream_thread_unlock: private _iostream_thread_lock_base {
  static inline void unlock(void)
  { 
#if defined(__USE_PTHREAD) 
    if (pthread_mutex_unlock(&mutex_))
      throw std::runtime_error("_iostream_thread_unlock::unlock: pthread_mutex_unlock error return");
#endif
  }
};
extern _iostream_thread_unlock log_unlock;

//! @endcond full_docs


} // namespace parallelUtil

//! @cond full_docs

namespace std{

template<typename _CharT, typename _Traits>
  inline basic_istream<_CharT,_Traits>& 
  operator>>(basic_istream<_CharT,_Traits>& __is, parallelUtil::_iostream_thread_lock __f)
  {
    __f.lock();
    return __is; 
  }

template<typename _CharT, typename _Traits>
  inline basic_istream<_CharT,_Traits>& 
  operator>>(basic_istream<_CharT,_Traits>& __is, parallelUtil::_iostream_thread_unlock __f)
  {
    __f.unlock();
    return __is; 
  }

template<typename _CharT, typename _Traits>
  inline basic_ostream<_CharT,_Traits>& 
  operator<<(basic_ostream<_CharT,_Traits>& __os, parallelUtil::_iostream_thread_lock __f)
  {
    __f.lock();
    return __os; 
  }  
  
template<typename _CharT, typename _Traits>
  inline basic_ostream<_CharT,_Traits>& 
  operator<<(basic_ostream<_CharT,_Traits>& __os, parallelUtil::_iostream_thread_unlock __f)
  {
    __f.unlock();
    return __os; 
  }
} // namespace std   

//! @endcond full_docs

#if defined(__USE_MPI)
#include "mpiDispatcher.h"
#endif

#if 0
  // requires _explicit_inclusion, uses TMatrix::read/writeBinary...
  #if defined(__USE_MPI)
  #include "parallelFunctor.h"
  #endif
#endif

#include "parallelUtil_inline.h"

#if !defined(EXCLUDE_TEMPLATE_BODIES)
#include "parallelUtil_template.h"
#endif

#endif // __parallelUtil__h
