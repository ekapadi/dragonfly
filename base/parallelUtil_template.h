#if !defined(__parallelUtil_template__h)
#define __parallelUtil_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/parallelUtil_template.h,v $

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

namespace parallelUtil{

// implicit assumption here is that "K" is a pointer to some object,
//   and that reinterpret_cast returns the unsigned long equivalent to the pointer:
template <class K>
inline size_t multiMutex::lock(const K& k)
{ return multiMutex::lock(reinterpret_cast<unsigned long>(k)); }

template <class K>
inline bool multiMutex::trylock(const K& k, size_t& nLock)
{ return multiMutex::trylock(reinterpret_cast<unsigned long>(k), nLock); }

  
  
template <class IT1, class IT2, class U1, class U2>    
	void loopIteratorList<IT1,IT2,U1,U2>::copy(const loopIteratorList& other)
	{
	 uIT1_ = other.uIT1_;
	 uIT2_ = other.uIT2_; 
	}

// Create iterator range subregion; makeup (if any) occurs in last split.
template <class IT1, class IT2, class U1, class U2>  
loopIteratorList<IT1,IT2,U1,U2> loopIteratorList<IT1,IT2,U1,U2>::split(size_t n_split, size_t N_splits, bool comb)const
{
 loopIteratorList<IT1,IT2,U1,U2> dest(*this);
 split(n_split, N_splits, comb, dest);
 return dest;
}

// Create iterator range subregion; makeup (if any) occurs in last split.
template <class IT1, class IT2, class U1, class U2>  
void loopIteratorList<IT1,IT2,U1,U2>::split(size_t n_split, size_t N_splits, bool comb, loopIteratorList& dest)const
{
 dest.clear(); // if dest is re-used, or created via "clone", need to eliminate existing contents.

 // number of iterators in a single split:
 const double splitQuotient( double(size())/double(N_splits) ); 
 const size_t splitSize = static_cast<size_t>( std::floor( splitQuotient ) );
 const bool makeUp = ( splitQuotient - double(splitSize) ) > std::numeric_limits<double>::epsilon();
	 
 if (!comb){ 
	 // Extract banded subrange:
	 // Note that split index origin is one.

	 // Calculate start and end iterators for range:
   // Note: this must work with non-random-access iterators:
	 typename U1::const_iterator it1Start(uIT1_.begin());
	 for(size_t i=0; i < (n_split-1)*splitSize; ++i) ++it1Start;

	 typename U1::const_iterator it1End(it1Start);
   for(size_t i=0; i < splitSize; ++i) ++it1End;		 

   size_t destSize = splitSize;		 

	 if ((n_split == N_splits) && makeUp){
	   // add one for cumulative "make-up" fraction:
	   ++it1End;
     ++destSize;			 
	 }

	 // copy the range to the destination:
	 dest.uIT1_.resize(destSize);
	 std::copy( it1Start, it1End, dest.uIT1_.begin() );

	 if (!uIT2_.empty()){
		 // Calculate start and end iterators for range:
		 typename U2::const_iterator it2Start(uIT2_.begin());
		 for(size_t i=0; i < (n_split-1)*splitSize; ++i) ++it2Start;

		 typename U2::const_iterator it2End(it2Start);
     for(size_t i=0; i < splitSize; ++i) ++it2End;		 

     // add one for cumulative "make-up" fraction:
		 if ((n_split == N_splits) && makeUp)
	  	 ++it2End;

     // note: destination uIT1_.size() == uIT2_.size() (destSize from above).
	   dest.uIT2_.resize(destSize);
	   std::copy( it2Start, it2End, dest.uIT2_.begin() );	 
	 }

 }
 else{
	 // Subset is spread-out, modulo N_splits:
	 // We assume iterators have "==" and "!=", but not necessarily ">".
	 // For this reason we span the entire range, rather than incrementing by N_splits.
	 // Note that split index origin is one.

   size_t destSize = 0;

	 // copy the subset to the destination:
   size_t n(0);
	 for(typename U1::const_iterator it1 = uIT1_.begin(), it1End = uIT1_.end(); 
		   it1 != it1End;
			 ++it1,
			 ++n)
	 if ((n_split-1) == n%N_splits){
		 dest.uIT1_.push_back(*it1);
     ++destSize;			 
	 }	 

	 if (!uIT2_.empty()){

	   dest.uIT2_.resize(destSize);
     size_t n(0);
		 for(typename U2::const_iterator it2 = uIT2_.begin(), it2End = uIT2_.end(); 
		     it2 != it2End;
				 ++it2,
				 ++n)
		 if ((n_split-1) == n%N_splits)
			 dest.uIT2_.push_back(*it2); 
	 }

 }

}

#if 1
// retain derived class type:
template <class IT1, class IT2, class U1, class U2>  
inline loopIteratorList<IT1,IT2,U1,U2>* loopIteratorList<IT1,IT2,U1,U2>::cloneSplit(size_t n_split, size_t N_splits, bool comb)const
{
 loopIteratorList* val(clone());
 split(n_split, N_splits, comb, *val);
 return val;
}
#endif

template <class IT1, class IT2, class U1, class U2>
loopIteratorList<IT1,IT2,U1,U2>* loopIteratorList<IT1,IT2,U1,U2>::clone(void)const
{
 return new loopIteratorList(*this);
} 

template <class IT1, class IT2, class U1, class U2>
void loopIteratorList<IT1,IT2,U1,U2>::write(std::ostream& os)const
{
 typename U1::const_iterator itU1=uIT1_.begin(), itU1End=uIT1_.end();
 if (!uIT2_.empty()){
   typename U2::const_iterator itU2=uIT2_.begin();
	 for(;
	     itU1 != itU1End;
			 ++itU1, ++itU2)
		 os<<"("<<*itU1<<","<<*itU2<<") ";	 
 }
 else{
	 for(;
	     itU1 != itU1End;
			 ++itU1)
		 os<<*itU1<<" ";	  
 } 
}
	
template <class IT1, class IT2, class U1, class U2>  	
loopIteratorList<IT1,IT2,U1,U2>::~loopIteratorList(void)
{}

template <class IT1, class IT2, class U1, class U2>  	
loopIteratorList<IT1,IT2,U1,U2>::loopIteratorList( const IT1& it1Start, const IT1& it1End, const IT2& it2Start )
{
#if 0
 // std::list<>::iterator does _not_ have a difference operator.
 size_t N_iterators_ = (it1End-it1Start); // it1End and it2End themselves are excluded: <iterator range> ::= [itStart, itEnd)
 uIT1_.resize(N_iterators_);
 uIT2_.resize(N_iterators_);

 IT1 it1(it1Start);
 IT2 it2(it2Start);
 typename U2::iterator itU2=uIT2_.begin();
 for(typename U1::iterator itU1=uIT1_.begin(), itU1End=uIT1_.end();
	   itU1 != itU1End;
		 ++itU1,
		 ++itU2,
		 ++it1,
		 ++it2){
   *itU1 = it1;
	 *itU2 = it2;
 }
#else
 // This implementation has potential inefficiencies due to possible re-allocation of uIT1_ and uIT2_.
 IT2 it2(it2Start);
 for(IT1 it1 = it1Start;
	   it1 != it1End;
		 ++it1,
		 ++it2){
	 uIT1_.push_back(it1);
	 uIT2_.push_back(it2);	 
 }
#endif 			 	
}
	
	// one iterator set only	
template <class IT1, class IT2, class U1, class U2>  	
loopIteratorList<IT1,IT2,U1,U2>::loopIteratorList( const IT1& it1Start, const IT1& it1End )
{
#if 0
 // std::list<>::iterator does _not_ have a difference operator.
 size_t N_iterators_ = (it1End-it1Start); // it1End itself is excluded: <iterator range> ::= [itStart, itEnd)
 uIT1_.resize(N_iterators_);

 IT1 it1(it1Start);
 for(typename U1::iterator itU1=uIT1_.begin(), itU1End=uIT1_.end();
	   itU1 != itU1End;
		 ++itU1,
		 ++it1)
   *itU1 = it1;
#else
// This implementation has potential inefficiencies due to possible re-allocation of uIT1_.
 for(IT1 it1 = it1Start;
	   it1 != it1End;
		 ++it1)
	 uIT1_.push_back(it1);
#endif 
}

	// copy constructor:
template <class IT1, class IT2, class U1, class U2>  	
  loopIteratorList<IT1,IT2,U1,U2>::loopIteratorList( const loopIteratorList& other )
		{ copy(other); }
		

  // default constructor:		
template <class IT1, class IT2, class U1, class U2>  	
loopIteratorList<IT1,IT2,U1,U2>::loopIteratorList( void )		
{}


template <class IT1, class IT2>	    
void loopIteratorArray<IT1,IT2>::copy(const loopIteratorArray& other)
{
 base_class::copy(*static_cast<const base_class*>(this));
}

template <class IT1, class IT2>  
typename loopIteratorArray<IT1,IT2>::base_class* loopIteratorArray<IT1,IT2>::clone(void)const
{
 return new loopIteratorArray<IT1,IT2>(*this);
}

template <class IT1, class IT2>	
loopIteratorArray<IT1,IT2>::~loopIteratorArray( void )
{}	

template <class IT1, class IT2>  
loopIteratorArray<IT1,IT2>::loopIteratorArray( const IT1& it1Start, const IT1& it1End, const IT2& it2Start )
  : loopIteratorList<IT1,IT2,std::vector<IT1>,std::vector<IT2> >(it1Start, it1End, it2Start)
{}

// one iterator set only	
template <class IT1, class IT2>  
loopIteratorArray<IT1,IT2>::loopIteratorArray( const IT1& it1Start, const IT1& it1End )
  : loopIteratorList<IT1,IT2,std::vector<IT1>,std::vector<IT2> >(it1Start, it1End)
{}

// copy constructor:
template <class IT1, class IT2>  
loopIteratorArray<IT1,IT2>::loopIteratorArray( const loopIteratorArray& other )
{ copy(other); }

// default constructor:		
template <class IT1, class IT2>  
loopIteratorArray<IT1,IT2>::loopIteratorArray( void )
{}



template <class WS, class IT1, class IT2, class U1, class U2>
bool OMP_parallel_for( void (*threadFunc)(WS&, const loopIteratorList<IT1,IT2,U1,U2>&),
                       WS& workspace,                        
											 const loopIteratorList<IT1,IT2,U1,U2>& indices )
{
 bool rVal(true);
 
 size_t OMP_NUM_THREADS_( parallelUtil::OMP_NUM_THREADS() );
 
 typedef  void (*thread_func_t)(WS&, const loopIteratorList<IT1,IT2,U1,U2>&);
 typedef  std::pair<WS*, const loopIteratorList<IT1,IT2,U1,U2>* > thread_parm_pair_t;
 typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;
 
 thread_functor_t vFunctor[OMP_NUM_THREADS_];
 pthread_t vThread[OMP_NUM_THREADS_]; 
  
 // create the threads:
 for(size_t nT = 0; nT < OMP_NUM_THREADS_; ++nT){
   vFunctor[nT].first = threadFunc;
	 vFunctor[nT].second.first = &workspace;
	 vFunctor[nT].second.second = indices.cloneSplit( nT+1, OMP_NUM_THREADS_, true );
 
	 // create the thread  
	 if (pthread_create( &(vThread[nT]), NULL, OMP_parallel_for_aux_const<WS,IT1,IT2,U1,U2>, &(vFunctor[nT])) != 0 ){
	   rVal = false;
		 break;
	 }	 
 }
 
 // join the threads:
 for(size_t nT = 0; nT < OMP_NUM_THREADS_; ++nT){
   // join the thread
	 if (pthread_join( vThread[nT], NULL ) != 0 ){
	   rVal =  false;	 
	 }	 
	 // cleanup
	 delete vFunctor[nT].second.second;
 } 

 return rVal;  
}

template <class WS, class IT1, class IT2, class U1, class U2>
void* OMP_parallel_for_aux_const(void *pvArgs)
{
 typedef  void (*thread_func_t)(WS&, const loopIteratorList<IT1,IT2,U1,U2>&);
 typedef  std::pair<WS*, const loopIteratorList<IT1,IT2,U1,U2>* > thread_parm_pair_t;
 typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;
 
 thread_functor_t &functor(*reinterpret_cast<thread_functor_t*>(pvArgs));
 (functor.first)( *functor.second.first, *functor.second.second );
 return NULL; 
}


template <class WS, class IT1, class IT2, class U1, class U2>
bool OMP_parallel_for( void (*threadFunc)(WS&, loopIteratorList<IT1,IT2,U1,U2>&),
                       WS& workspace,                        
											 loopIteratorList<IT1,IT2,U1,U2>& indices )
{
 bool rVal(true);
 
 size_t OMP_NUM_THREADS_( parallelUtil::OMP_NUM_THREADS() );
 
 typedef  void (*thread_func_t)(WS&, loopIteratorList<IT1,IT2,U1,U2>&);
 typedef  std::pair<WS*, loopIteratorList<IT1,IT2,U1,U2>* > thread_parm_pair_t;
 typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;
 
 thread_functor_t vFunctor[OMP_NUM_THREADS_];
 pthread_t vThread[OMP_NUM_THREADS_]; 
	
 // create the threads:
 for(size_t nT = 0; nT < OMP_NUM_THREADS_; ++nT){
   vFunctor[nT].first = threadFunc;
	 vFunctor[nT].second.first = &workspace;
	 vFunctor[nT].second.second = indices.cloneSplit( nT+1, OMP_NUM_THREADS_, true );
	 
	 // create the thread  
	 if (pthread_create( &(vThread[nT]), NULL, OMP_parallel_for_aux_non_const<WS,IT1,IT2,U1,U2>, &(vFunctor[nT])) != 0 ){
	   rVal = false;
		 break;
	 }	 
 }
 
 // join the threads:
 for(size_t nT = 0; nT < OMP_NUM_THREADS_; ++nT){
   // join the thread
	 if (pthread_join( vThread[nT], NULL ) != 0 ){
	   rVal =  false;	 
	 }	 
	 // cleanup
	 delete vFunctor[nT].second.second;
 } 

 return rVal;  
}

template <class WS, class IT1, class IT2, class U1, class U2>
void* OMP_parallel_for_aux_non_const(void *pvArgs)
{
 typedef  void (*thread_func_t)(WS&, loopIteratorList<IT1,IT2,U1,U2>&);
 typedef  std::pair<WS*, loopIteratorList<IT1,IT2,U1,U2>* > thread_parm_pair_t;
 typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;
 
 thread_functor_t &functor(*reinterpret_cast<thread_functor_t*>(pvArgs));
 (functor.first)( *functor.second.first, *functor.second.second );
 return NULL; 
}



// ----------------------------------------------------------------------------------------------------------------------------------------

// analogous parallel loop implementation using member-functions:
//   MEMBER_FUNC is of type
//     RVAL (CLASS_::*MEMBER_FUNC)(<typedef iterator_traits<IT>::reference OR " "::value_type>)[<optional "const" spec>][<optional "throw" spec>]
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT>
std::vector<RVAL> parallel_for(
  CLASS_ *instance, 
  MEMBER_FUNC threadFunc, 
  IT begin_, IT end_, size_t N_THREAD)
{
  bool status(true);

  // information required to call member-method from each thread:
  typedef loopIteratorArray<IT,IT> IT_list;
  typedef  MEMBER_FUNC thread_func_t;
  typedef  std::pair<CLASS_*, typename IT_list::base_class* > thread_parm_pair_t;
  typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;

  // return-value information:
  //   bool => value is std::vector<RVAL>*
  //     otherwise it is a std::runtime_error* (e.g. a thrown exception)
  typedef std::pair<bool, void*> thread_result_t;
  
  IT_list indices(begin_, end_);
  thread_functor_t vFunctor[N_THREAD];
  pthread_t vThread[N_THREAD]; 

  #if 0
  // *** DEBUG ***
  std::cout<<"parallelUtil::parallel_for: creating "<<N_THREAD<<" threads..."<<std::endl;
  #endif
  
  // create the threads:
  for(size_t nt = 0; status && (nt < N_THREAD); ++nt){
    vFunctor[nt].first = threadFunc;
	  vFunctor[nt].second.first = instance;
	  vFunctor[nt].second.second = indices.cloneSplit( nt+1, N_THREAD, true );

	  // create the thread  
	  if (pthread_create(&(vThread[nt]), NULL, parallel_for_aux<CLASS_,RVAL,MEMBER_FUNC,IT>, &(vFunctor[nt])) != 0)
      status = false;  
  }

  std::vector<RVAL> result;
  result.reserve(indices.size());
  
  // join the threads:
  for(size_t nt = 0; status && (nt < N_THREAD); ++nt){

    // gcc-4.1.2 port: strict-aliasing rules:
    // join the thread
	  void_ptr_union<thread_result_t> p_result_;
    if (pthread_join(vThread[nt], &(p_result_.void_ptr)) != 0){ 
      status = false;
      break;
    }
    // append to overall result:
    if (p_result_.ptr->first){
      std::vector<RVAL>* pv_result_(reinterpret_cast<std::vector<RVAL>*>(p_result_.ptr->second));
      result.insert(result.end(), pv_result_->begin(), pv_result_->end());     
      delete pv_result_;
      delete p_result_.ptr;
    }
    else{
      // thread function threw exception:
      std::runtime_error *px_(reinterpret_cast<std::runtime_error*>(p_result_.ptr->second));
      std::runtime_error x(*px_);
      delete px_;
      delete p_result_.ptr;
      throw std::runtime_error(x);
    }
	  delete vFunctor[nt].second.second;    

  } 

  #if 0
  // *** DEBUG ***
  std::cout<<"parallelUtil::parallel_for: "<<N_THREAD<<" threads joined."<<std::endl;
  #endif

  if (!status)
    throw std::runtime_error("parallel_for: pthread method error return");

  return result;    
}  

template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT>
void* parallel_for_aux(void *pvArgs)
{
  // information required to call member-method from each thread(from above):
  typedef loopIteratorArray<IT,IT> IT_list;
  typedef  MEMBER_FUNC thread_func_t;
  typedef  std::pair<CLASS_*, typename IT_list::base_class* > thread_parm_pair_t;
  typedef  std::pair< thread_func_t, thread_parm_pair_t > thread_functor_t;

  // return-value information:
  //   bool => value is std::vector<RVAL>*
  //     otherwise it is a std::runtime_error* (e.g. a thrown exception)
  typedef std::pair<bool, void*> thread_result_t;
  thread_result_t *presult(new thread_result_t(true, NULL));
  
  try{
    thread_functor_t &functor(*reinterpret_cast<thread_functor_t*>(pvArgs));
    typename IT_list::base_class &indices(*functor.second.second);

    std::vector<RVAL> *pv_result(new std::vector<RVAL>(indices.size()));
    typename IT_list::base_class::iterator itN = indices.begin();
    for(typename std::vector<RVAL>::iterator itR = (*pv_result).begin(), itREnd = (*pv_result).end();
        itR != itREnd;
        ++itR, ++itN)
      *itR = ((functor.second.first)->*(functor.first))(*(*itN));

    presult->second = reinterpret_cast<void*>(pv_result); 
  }
  catch (const std::runtime_error &x){
    presult->first = false; // indicate ERROR return
    assert(presult->second == NULL);
    presult->second = reinterpret_cast<void*>(new std::runtime_error(x));
  }
  
  return presult;
}


// +1 ARG version:
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT, class ARG>
std::vector<RVAL> parallel_for(
  CLASS_ *instance, 
  MEMBER_FUNC threadFunc, 
  IT begin_, IT end_, ARG arg, size_t N_THREAD)
{
  bool status(true);

  // information required to call member-method from each thread:
  typedef loopIteratorArray<IT,IT> IT_list;
  typedef  MEMBER_FUNC thread_func_t;
  typedef  std::pair <typename IT_list::base_class*, ARG*> thread_parm_args_t; 
  typedef  std::pair<CLASS_*, thread_parm_args_t > thread_parms_t;
  typedef  std::pair< thread_func_t, thread_parms_t > thread_functor_t;

  // return-value information:
  //   bool => value is std::vector<RVAL>*
  //     otherwise it is a std::runtime_error* (e.g. a thrown exception)
  typedef std::pair<bool, void*> thread_result_t;
  
  IT_list indices(begin_, end_);
  thread_functor_t vFunctor[N_THREAD];
  pthread_t vThread[N_THREAD]; 

  // create the threads:
  for(size_t nt = 0; status && (nt < N_THREAD); ++nt){
    vFunctor[nt].first = threadFunc;
	  vFunctor[nt].second.first = instance;
	  vFunctor[nt].second.second.first = indices.cloneSplit( nt+1, N_THREAD, true );
	  vFunctor[nt].second.second.second = &arg;
    
	  // create the thread  
	  if (pthread_create(&(vThread[nt]), NULL, parallel_for_aux<CLASS_,RVAL,MEMBER_FUNC,IT,ARG>, &(vFunctor[nt])) != 0)
      status = false;  
  }

  std::vector<RVAL> result;
  result.reserve(indices.size());
  
  // join the threads:
  for(size_t nt = 0; status && (nt < N_THREAD); ++nt){

    // gcc-4.1.2 port: strict-aliasing rules:
    // join the thread
	  void_ptr_union<thread_result_t> p_result_;
    if (pthread_join(vThread[nt], &(p_result_.void_ptr)) != 0){ 
      status = false;
      break;
    }
    // append to overall result:
    if (p_result_.ptr->first){
      std::vector<RVAL>* pv_result_(reinterpret_cast<std::vector<RVAL>*>(p_result_.ptr->second));
      result.insert(result.end(), pv_result_->begin(), pv_result_->end());     
      delete pv_result_;
      delete p_result_.ptr;
    }
    else{
      // thread function threw exception:
      std::runtime_error *px_(reinterpret_cast<std::runtime_error*>(p_result_.ptr->second));
      std::runtime_error x(*px_);
      delete px_;
      delete p_result_.ptr;
      throw x;
    }
	  delete vFunctor[nt].second.second.first;    

  } 

  if (!status)
    throw std::runtime_error("parallel_for: pthread method error return");

  return result;    
}  
  
template <class CLASS_, class RVAL, class MEMBER_FUNC, class IT, class ARG>
void* parallel_for_aux(void *pvArgs)
{
  // information required to call member-method from each thread:
  typedef loopIteratorArray<IT,IT> IT_list;
  typedef  MEMBER_FUNC thread_func_t;
  typedef  std::pair <typename IT_list::base_class*, ARG*> thread_parm_args_t; 
  typedef  std::pair<CLASS_*, thread_parm_args_t > thread_parms_t;
  typedef  std::pair< thread_func_t, thread_parms_t > thread_functor_t;


  // return-value information:
  //   bool => value is std::vector<RVAL>*
  //     otherwise it is a std::runtime_error* (e.g. a thrown exception)
  typedef std::pair<bool, void*> thread_result_t;
  thread_result_t *presult(new thread_result_t(true, NULL));
  
  try{
    thread_functor_t &functor(*reinterpret_cast<thread_functor_t*>(pvArgs));
    typename IT_list::base_class &indices(*functor.second.second.first);
    ARG &arg(*functor.second.second.second);
    
    std::vector<RVAL> *pv_result(new std::vector<RVAL>(indices.size()));
    typename IT_list::base_class::iterator itN = indices.begin();
    for(typename std::vector<RVAL>::iterator itR = (*pv_result).begin(), itREnd = (*pv_result).end();
        itR != itREnd;
        ++itR, ++itN)
      *itR = ((functor.second.first)->*(functor.first))(*(*itN), arg);
      
    presult->second = reinterpret_cast<void*>(pv_result); 
  }
  catch (const std::runtime_error &x){
    presult->first = false; // indicate ERROR return
    assert(presult->second == NULL);
    presult->second = reinterpret_cast<void*>(new std::runtime_error(x));
  }
  
  return presult;
}
// ----------------------------------------------------------------------------------------------------------------------------------------



} // namespace parallelUtil

#endif // __parallelUtil_template__h
