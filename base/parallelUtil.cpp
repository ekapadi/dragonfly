// $Source: /usr/data0/leipzig_work/tmat_cvs/src/parallelUtil.cpp,v $

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


#include <portability.h>

#if defined(__ICC)
  #pragma warning(disable:981) // operands are evaluated in unspecified order
	#if defined(NDEBUG)
	  #pragma warning(disable:593) // variable "mutexErr" was set but never used
	#endif
#endif

#if defined(__USE_PTHREAD)
  #include <pthread.h>
#endif

#if defined(__USE_MPI)
  #include <mpi.h>
#endif

#include <assert.h>
#include <errno.h>
#include <asm/unistd.h> // for "by-number" syscalls to obtain system: "gettid"
#include <stdexcept>
#include <typeinfo>
#include <cstdio>

#include <iostream>
#include <sstream>
#include <iomanip>

#include <cmath>
#include <complex>
#include <vector>
#include <list>
#include <algorithm>

#include <typeinfo>
#include <tr1/type_traits>

#include "commUtil.h"
using commUtil::abstractCommHandle;
using commUtil::fileHandle;
#if defined(__USE_MPI)
using commUtil::processCommHandle;
#endif
using commUtil::open;
using commUtil::read;
using commUtil::write;
using commUtil::close;
using commUtil::readBinary;
using commUtil::writeBinary;

// MPI debugging output messages:
#if !defined(__NO_DEBUG_OUTPUT__)
#define __MPI_DEBUG__
#endif

#include "parallelUtil.h"

#include "statusUtil.h"

using std::cout;
using std::endl;

namespace parallelUtil{

#if defined(_GLIBCXX_IOSTREAM) || (defined(__PGI) && defined(_STLP_IOSTREAM))
void multiMutex::namedLock::write(std::ostream& os)const
{
   os<<"key: "<<key_;
   if (refCount_ > 0)
     os<<"(ref-count: "<<refCount_<<")";
     
   int status __attribute_unused__ (0);
   if ((status = pthread_mutex_trylock(const_cast<pthread_mutex_t*>(&mutex_))))
     os<<"--LOCKED--("<<tid_<<")";
   else
   if (pthread_mutex_unlock(const_cast<pthread_mutex_t*>(&mutex_)))
     throw std::runtime_error("multiMutex::namedLock::write: pthread_mutex_unlock error return");       
}
#endif

multiMutex::namedLock::namedLock(void)
	: mutex_(CONST_PTHREAD_ERRORCHECK_MUTEX_INITIALIZER_NP), // (CONST_PTHREAD_MUTEX_INITIALIZER),
    tid_(0), 
		key_(static_cast<unsigned long>(NULL)), refCount_(0)		
{}

			
size_t multiMutex::lock(unsigned long key)
{
 #if 0
 // *** DEBUG ***
 size_t state(0);
 std::cout<<this<<": lock, at-entry (acquiring lock for "<<key<<"): ";
 write(std::cout);
 std::cout<<endl;
 #endif	
  
 size_t nLock(0);
 namedLock *pLock(NULL);

 if (pthread_mutex_lock(&controlMutex_))
   throw std::runtime_error("multiMutex::lock: pthread_mutex_lock error return");

 for (size_t n=0; n < N_locks_; ++n){
	 namedLock& cLock(vLock_[n]);

	 // Either find an already assigned named lock, or
	 //   return the first available named lock.
	 if (key == cLock.key_){
     #if 0
     // *** DEBUG ***
     state = 1;
     cout<<"matching lock exists"<<endl;
     #endif
		 // there already is a matching lock:
		 pLock = &cLock;
		 nLock = n;
		 break;
	 }	 
	 else	 
	 if ((NULL == pLock) && (static_cast<unsigned long>(NULL) == cLock.key_)){
     #if 0
     // *** DEBUG ***
     state = 2;
     cout<<"lock available"<<endl;
     #endif		 
     // first available lock:
		 pLock = &cLock;
		 nLock = n;
	 }
 }

 if (NULL == pLock)
   throw std::runtime_error("multiMutex::lock: no new lock available");
 // This case should probably be a throw -- it will segfault next...
 // Either there are fewer locks than threads, or a thread took a lock
 //   for another key and didn't give it back.
 
 if (static_cast<unsigned long>(NULL) == pLock->key_)
	 pLock->key_ = key; // name the lock.

 // note: this also will work with recursive mutexes:		 
 ++pLock->refCount_;

 if (pthread_mutex_unlock(&controlMutex_))
   throw std::runtime_error("multiMutex::lock: pthread_mutex_unlock error return");

 #if 0
 // *** DEBUG ***
 std::cout<<"Blocking for named-lock: ";
 pLock->write(std::cout);
 std::cout<<std::endl;
 #endif
 #if 1
 // block until available:
 if (pthread_mutex_lock(&pLock->mutex_))
   throw std::runtime_error("multiMutex::lock: pthread_mutex_lock error return");
 pLock->tid_ = pthread_self();
 assert(pLock->tid_ != 0);
 #else
 // *** DEBUG ***
 // is the mutex just locked??!!
 int status(0);
 if (status = pthread_mutex_trylock(&pLock->mutex_)){
   std::ostringstream oss;
   oss<<"multiMutex::lock: pthread_mutex_trylock returns "<<status<<": "<<strerror(status)<<"\n"
     <<" full mutex: ";
   write(oss);
   throw std::runtime_error(oss.str()); 
 }
 pLock->tid_ = pthread_self();
 assert(pLock->tid_ != 0); 
 #endif
   
 #if 0
 // *** DEBUG ***
 std::cout<<this<<": lock, at-exit (state: "<<state<<"): ";
 write(std::cout);
 std::cout<<endl;
 #endif	 
   
 return nLock;
}


bool multiMutex::trylock(unsigned long key, size_t& nLock)
{
  bool status(true);
  namedLock *pLock(NULL);

  if (pthread_mutex_lock(&controlMutex_))
    throw std::runtime_error("multiMutex::lock: pthread_mutex_lock error return");

  for (size_t n=0; n < N_locks_; ++n){
	  namedLock& cLock(vLock_[n]);

	  // Either find an already assigned named lock, or
	  //   return the first available named lock.
	  if (key == cLock.key_){
		  // there already is a matching lock:
		  pLock = &cLock;
		  nLock = n;
		  break;
	  }	 
	  else	 
	  if ((NULL == pLock) && (static_cast<unsigned long>(NULL) == cLock.key_)){
      // first available lock:
		  pLock = &cLock;
		  nLock = n;
	  }
  }

  if (NULL == pLock)
    throw std::runtime_error("multiMutex::lock: no new lock available");
  // This case should probably be a throw -- it will segfault next...
  // Either there are fewer locks than threads, or a thread took a lock
  //   for another key and didn't give it back.

  if (static_cast<unsigned long>(NULL) == pLock->key_)
	  pLock->key_ = key; // name the lock.

  // note: this also will work with recursive mutexes:		 
  ++pLock->refCount_;

  if (pthread_mutex_unlock(&controlMutex_))
    throw std::runtime_error("multiMutex::lock: pthread_mutex_unlock error return");

  // attempt to lock the mutex:
  int N_STAT(0);
  if ((N_STAT = pthread_mutex_trylock(&pLock->mutex_))){
    if (N_STAT == EBUSY)
      status = false; // mutex is locked
    else
      throw std::runtime_error("multiMutex::trylock: pthread_mutex_trylock error return");
  }
  else{
    pLock->tid_ = pthread_self();
    assert(pLock->tid_ != 0);
  }

  return status;
}



void multiMutex::unlock(size_t nLock)
{
 #if 0
 // *** DEBUG ***
 std::cout<<this<<": unlock, at-entry: ";
 write(std::cout);
 std::cout<<endl;
 #endif	  
 
 // Note: the lockObject _shall_ belong to this multimutex!
 if( nLock >= N_locks_ )
   throw std::runtime_error("multiMutex::unlock: key out-of-range");
	 
 namedLock &cLock(vLock_[nLock]);
 // <number of unlocks> == <number of locks>
 if (vLock_[nLock].refCount_ == 0)
   throw std::runtime_error("multiMutex::unlock: too many unlocks/locks");
	 
 if (pthread_mutex_unlock(&cLock.mutex_))
   throw std::runtime_error("multiMutex::unlock: pthread_mutex_unlock error return");	
 
 // maintain the mutex list:
 if (pthread_mutex_lock(&controlMutex_))
   throw std::runtime_error("multiMutex::unlock: pthread_mutex_lock error return");	

 if (0 == --cLock.refCount_){
	 cLock.key_ = static_cast<unsigned long>(NULL); // remove the name.
   cLock.tid_ = 0;
 }
 
 if (pthread_mutex_unlock(&controlMutex_))
   throw std::runtime_error("multiMutex::unlock: pthread_mutex_unlock error return");		  
   
 #if 0
 // *** DEBUG ***
 std::cout<<this<<": unlock, at-exit: ";
 write(std::cout);
 std::cout<<endl;
 #endif	    
}

void multiMutex::allocLocks(size_t N_locks)
{
  if (vLock_ != NULL)
    delete[] vLock_;
	N_locks_ = N_locks;
	vLock_ = new namedLock[N_locks];
  #if 0
  // *** DEBUG *** spurious memory access trashes mutex 0:
  std::cout
    <<"multiMutex::mutiMutex: &(vLock[0]) start pad0: "<<reinterpret_cast<void*>(&(vLock_[0].__PAD0__[0]))<<"\n"
    <<"                       &(vLock[0]) end   pad0: "<<reinterpret_cast<void*>(&(vLock_[0].__PAD0__[sizeof(pthread_mutex_t)]))<<"\n"
    <<"multiMutex::mutiMutex: &(vLock[0]) start pad1: "<<reinterpret_cast<void*>(&(vLock_[0].__PAD1__[0]))<<"\n"
    <<"                       &(vLock[0]) end   pad1: "<<reinterpret_cast<void*>(&(vLock_[0].__PAD1__[sizeof(pthread_mutex_t)]))<<std::endl;   
  #endif
}

#if defined(_GLIBCXX_IOSTREAM) || (defined(__PGI) && defined(_STLP_IOSTREAM))
void multiMutex::write(std::ostream& os)const
{
  os<<"[";
  for(size_t n = 0; n < N_locks_; ++n){
    if (n > 0)
      os<<", ";
    vLock_[n].write(os);
  }
  os<<"]";  
}
#endif

multiMutex::multiMutex(size_t N_locks)
	:controlMutex_(CONST_PTHREAD_MUTEX_INITIALIZER),
	 vLock_(NULL)
{
  allocLocks(N_locks);
}

multiMutex::~multiMutex(void)
{
 delete[] vLock_;
 // mark so that it is _obvious_ the multiMutex has been destroyed:
 if (pthread_mutex_destroy(&controlMutex_))
   throw std::runtime_error("multiMutex::~multiMutex: pthread_mutex_destroy error return");
}		

// get the environment variable: OMP_NUM_THREADS (if it doesn't exist, return 1):
size_t OMP_NUM_THREADS(void)
{
  size_t val(1);
	
	char *szOMP_NUM_THREADS = getenv("OMP_NUM_THREADS");
	if (NULL != szOMP_NUM_THREADS){
	  char *szTest;
		long nThreads = std::strtol(szOMP_NUM_THREADS, &szTest, 10);
		if ( nThreads > 0 && szOMP_NUM_THREADS[0] != '\0' && szTest[0] == '\0' )
		  val = nThreads;
	}

  return val;
}


// return true if the calling thread is the process "root" thread:
// (this could be an inline, but the requirement for "asm/unistd.h" would seem to cause undesirable global namespace pollution.)
bool is_process_root(void)
{
  return (syscall(__NR_getpid) == syscall(__NR_gettid));
}


#if defined(__USE_PTHREAD)
pthread_mutex_t _iostream_thread_lock_base::mutex_ = PTHREAD_MUTEX_INITIALIZER;
#endif
_iostream_thread_lock log_lock;
_iostream_thread_unlock log_unlock;

} // namespace parallelUtil
