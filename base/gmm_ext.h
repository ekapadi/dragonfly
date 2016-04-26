#if !defined(__gmm_ext__h)
#define __gmm_ext__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/gmm_ext.h,v $

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


#define __specialize_POD_binary__

/** \file gmm_ext.h
  * \brief Extension classes and methods for namespace gmm.
  *
  * The namespace "gmm_ext" is envisioned as a small group of extension methods to "gmm".
  * These methods are implemented using "gmm" style, rather than "numericalFunctor" style of interface.
  * At the point this becomes a _large_ group of extensions, we should consider switching to another
  *  matrix/vector template library.
  *
  *
  * Class list:
  *   "basic_index_ref": an unsorted sub-index that stores only the iterators.
  *   "gmm::linalg_traits<std::dequeue<T>>" : only declared when #include \<deque\> _prior_
  *
  * Current list of methods: 
  *   "QRP_factor": QR reduction with complete pivoting 
  *   "SVD_jacobi": Jacobi method for singular-value-decomposition(SVD)
  *   "SVD_QR": High relative accuracy QR method for SVD
  *
  */

/*
 *
 * Unfortunately, gmm::unsorted_sub_index copies its index-container.  Ref-counting is implemented, but _not_ enforced for
 *   any of the using classes like gmm::sub_col_matrix, which stores copies of the sub-index (rather than pointers).
 * Even worse is the fact that  gmm::tab_ref_index_ref, which is the default row or column sub-vector return type 
 *   for unsorted_sub_index derived sub-matrices simply takes the "std::vector<size_t>::const_iterator" begin and end
 *   of the unsorted_sub_index's index-vector.  The end effect is that this does not correctly update the reference-counting
 *   through a returned reference (e.g. gmm::mat_col( <sub_matrix>, col_num )), and these iterators end up pointing to
 *   an index-vector that has already been destroyed.
 *
 * For these reasons I implement the necessary _minimal_ classes to have an "unsorted_sub_index_ref" which simply retains
 *   appropriate iterators to any suitable container (that _shall_ not be destroyed while the sub_index_ref is valid).
 * To work correctly including required template specializations of "gmm" templates, these must be defined in the "gmm" namespace.
 *
 */

// miscellaneous formatted I/O methods:
namespace gmm{

// place these std::pair methods here, for the moment (output of vector goes through "gmm"):
template <class T1, class T2>
inline std::ostream& operator<<(std::ostream& os, const std::pair<T1,T2>& p_);

}

//! @cond full_docs
// selector for POD value_type:
namespace gmm{

  #if 0 
  // ------------------------- non-TR1 version: ----------------------------------
  /**
    @brief selector for POD value_type:
   */
  
  template <class T>
  struct is_POD{
    typedef linalg_false POD_type;
  }; 

  template <>
  struct is_POD<bool>{
    typedef linalg_true POD_type;
  };  

  template <>
  struct is_POD<long>{
    typedef linalg_true POD_type;
  };

  template <>
  struct is_POD<int>{
    typedef linalg_true POD_type;
  };
  
  template <>
  struct is_POD<size_t>{
    typedef linalg_true POD_type;
  };
  
  template <>
  struct is_POD<double>{
    typedef linalg_true POD_type;
  };
  
  template <>
  struct is_POD<std::complex<double> >{
    typedef linalg_true POD_type;
  };
  #else
  // -------------------------- version using TR1: -------------------------------
  
  // map the TR1 abstract true / false type to the gmm types:
  template <class T>
  struct linalg_type_mapper{
    typedef abstract_null_type type;
  };
  
  template <>
  struct linalg_type_mapper<std::false_type>{
    typedef linalg_false type;
  };  
  
  template <>
  struct linalg_type_mapper<std::true_type>{
    typedef linalg_true type;
  };
    
  template <class T>
  struct is_POD{
    typedef typename linalg_type_mapper<typename std::is_pod<T>::type>::type POD_type;
  };   

  #endif
  
} // namespace gmm

//! @endcond


#if 1 || defined(_GLIBCXX_DEQUE)
// only declare if std::deque is declared prior:
namespace gmm{
  template <typename T, typename alloc>
  struct linalg_traits<std::deque<T, alloc> > {
    typedef std::deque<T, alloc> this_type;
    typedef this_type origin_type;
    typedef linalg_false is_reference;
    typedef abstract_vector linalg_type;
    typedef T value_type;
    typedef T& reference;
    typedef typename this_type::iterator iterator;
    typedef typename this_type::const_iterator const_iterator;
    typedef abstract_dense storage_type;
    typedef linalg_true index_sorted;
    static size_type size(const this_type &v) { return v.size(); }
    static iterator begin(this_type &v) { return v.begin(); }
    static const_iterator begin(const this_type &v) { return v.begin(); }
    static iterator end(this_type &v) { return v.end(); }
    static const_iterator end(const this_type &v) { return v.end(); }
    static origin_type* origin(this_type &v) { return &v; }
    static const origin_type* origin(const this_type &v) { return &v; }
    static void clear(origin_type*, const iterator &it, const iterator &ite)
    { std::fill(it, ite, number_traits<value_type>::zero()); }
    static void do_clear(this_type &v) { std::fill(v.begin(), v.end(), number_traits<T>::zero() ); }
    static value_type access(const origin_type *, const const_iterator &it,
			     const const_iterator &, size_type i)
    { return it[i]; }
    static reference access(origin_type *, const iterator &it,
			    const iterator &, size_type i)
    { return it[i]; }
    static void resize(this_type &v, size_type n) { v.resize(n); }
  };
} // namespace gmm

namespace std {
  template <typename T> ostream &operator <<
  (std::ostream &o, const deque<T>& m) { gmm::write(o,m); return o; }
} // namespace std

namespace gmm {
  template <typename T>
  inline size_type nnz(const std::deque<T>& l) { return l.size(); }
}

#endif // defined(_GLIBCXX_DEQUE)


namespace gmm{
  
  class basic_index_ref{
		public:
		
		typedef size_t size_type;
		
    class generic_iterator_base{
      protected:
			  // allow strict-aliasing: move actual pointer to derived class generic_iterator_wrapper<IT>:
        #if 0
        void *pvIT_;
        #endif
		  public:
				typedef basic_index_ref::size_type size_type;
    		typedef size_type value_type;
    		typedef size_type *pointer;
    		typedef size_type &reference;
    		typedef ptrdiff_t difference_type;
    		typedef std::random_access_iterator_tag  iterator_category;

				// implementation note:
        //   virtuals to be used in "++", "--", "+=", "-="
        //   must have _void_ return type, as the derived-class iterator operators
        //   do not normally return the base-class pointer
        //   (which is the only possibility to return an abstract base class).
        virtual void increment(void)=0;
    		virtual void decrement(void)=0;
    		virtual void add(difference_type n)=0;
    		virtual void subtract(difference_type n)=0;
        virtual difference_type operator-(const generic_iterator_base& other)const=0;

    		// virtual value_type& operator *() const=0;
				virtual const value_type& operator *() const=0;
    		// virtual value_type& operator [](int i) const=0;
        virtual const value_type& operator [](int n) const=0;
				
    		virtual bool operator ==(const generic_iterator_base &other) const=0;
    		virtual bool operator !=(const generic_iterator_base &other) const=0;
    		virtual bool operator < (const generic_iterator_base &other) const=0;

        #if 0
        virtual void copy(const generic_iterator_base& other)=0;
				#endif
        
				virtual generic_iterator_base* clone(void)const=0;
				
				#if 0
				template <class IT>
				virtual void init(const IT& it)=0;
				#endif
				
				#if 0
				generic_iterator_base& operator=(const generic_iterator_base& other)
				{ copy(other); return *this; }
				#endif
				
				generic_iterator_base(void)
				{ }

        #if 0
				generic_iterator_base(const generic_iterator_base& other)
				{ copy(other); }
				#endif
								
        virtual ~generic_iterator_base(void)
				{ }
		};
				
    template <class IT>
		class generic_iterator_wrapper: public generic_iterator_base{
      public:
			  typedef generic_iterator_base base_class; 
						
      private:
        // allow strict-aliasing:
        IT *pit_;

		    IT &it_(void) 
				{return *pit_; }
				
        const IT &it_(void)const 
				{return *pit_; }	
							
			  const IT &it_(const base_class& other)const
				{ return *(static_cast<const generic_iterator_wrapper<IT>&>(other).pit_); }
				
			public:
        typedef base_class::size_type size_type;
    		typedef size_type value_type;
    		typedef size_type *pointer;
    		typedef size_type &reference;
    		typedef ptrdiff_t difference_type;
    		typedef std::random_access_iterator_tag  iterator_category;
 
    		virtual void increment(void)  { it_()++; }
    		virtual void decrement(void)  { it_()--; }
    		virtual void add(difference_type n) { it_() += n; }
    		virtual void subtract(difference_type n) { it_() -= n; }
    		virtual difference_type operator -(const base_class &other) const 
				{ return it_() - it_(other); }

    		// virtual value_type& operator *() const { return *it_(); }
    		virtual const value_type& operator *() const { return *it_(); }
				
    		// virtual value_type& operator [](int i) const { return *(it_()+i); }
    		virtual const value_type& operator [](int n) const { return *(it_()+n); }

    		virtual bool operator ==(const base_class &other) const { return (it_() == it_(other)); }
    		virtual bool operator !=(const base_class &other) const { return (it_() != it_(other)); }
    		virtual bool operator < (const base_class &other) const { return (it_() < it_(other)); }


        virtual ~generic_iterator_wrapper(void)
				{ 
				  if (NULL != pit_)
					  delete pit_;
				  pit_ = NULL;
				}

        void copy(const generic_iterator_wrapper& other)
	      { 
				  if (NULL != pit_){
					  delete pit_;
						pit_ = NULL;
					}
					if (NULL != other.pit_)
				    pit_ = new IT(it_(other));
				}
			
			  virtual base_class* clone(void)const
				{ return new generic_iterator_wrapper<IT>(*this); }
			
							
				generic_iterator_wrapper(void)
          : pit_(NULL)
				{ }
				
				generic_iterator_wrapper(const IT& it)
				  : pit_(new IT(it))
        { }	

		    generic_iterator_wrapper(const generic_iterator_wrapper& other)
				  : pit_(NULL)
        { 
				  copy(other); 
				}
		};				
    
		class generic_iterator{
		  private:
			  generic_iterator_base *piterator_;
				
			  generic_iterator_base& it_(void)
				{ return *piterator_; }
				const generic_iterator_base& it_(void)const
				{ return *piterator_; }
			public:
 
        typedef generic_iterator iterator;
				
				typedef generic_iterator_base::size_type size_type;
    		typedef size_type value_type;
    		typedef size_type *pointer;
    		typedef size_type &reference;
    		typedef ptrdiff_t difference_type;
    		typedef std::random_access_iterator_tag  iterator_category;
 
    		iterator operator ++(int) { iterator tmp = *this; it_().increment(); return tmp; }
    		iterator operator --(int) { iterator tmp = *this; it_().decrement(); return tmp; }
    		iterator &operator ++()   { it_().increment(); return *this; }
    		iterator &operator --()   { it_().decrement(); return *this; }
    		iterator &operator +=(difference_type n) { it_().add(n); return *this; }
    		iterator &operator -=(difference_type n) { it_().subtract(n); return *this; }
    		iterator operator +(difference_type n) const 
    		{ iterator itt = *this; itt.it_().add(n); return itt; }
    		iterator operator -(difference_type n) const
    		{ iterator itt = *this; itt.it_().subtract(n); return itt; }
    		difference_type operator -(const iterator &other) const 
				{ return it_() - other.it_(); }

    		// value_type& operator *() const { return *it_(); }
     		const value_type& operator *() const { return *it_(); }
    		
				// value_type& operator [](int n) const { return it_()[n]; }
				const value_type& operator [](int n) const { return it_()[n]; }

    		bool operator ==(const iterator &other) const { return (it_() == other.it_()); }
    		bool operator !=(const iterator &other) const { return (it_() != other.it_()); }
    		bool operator < (const iterator &other) const { return (it_() < other.it_()); }

			
			  void copy(const generic_iterator& other)
				{
				 if (NULL != piterator_){
				   delete piterator_;
					 piterator_ = NULL;
				 }
				 if (NULL != other.piterator_)
				   piterator_ = other.piterator_->clone();
				}
				
				generic_iterator* clone(void)const
				{ return new generic_iterator(*this); }
				
				generic_iterator& operator=(const generic_iterator& other)
				{ copy(other); return *this; }
				
				~generic_iterator(void)
				{
				 if (NULL != piterator_)
				   delete piterator_;
				 piterator_ = NULL;
				}
				
			  generic_iterator(void)
				  :piterator_(NULL)
				{ }
				
				template <class IT>
				explicit
				generic_iterator(const IT& it)
				  :piterator_(new generic_iterator_wrapper<IT>(it))
				{ }
				
				generic_iterator(const generic_iterator& other)
				  :piterator_(NULL)
				{ copy(other); }
		};
		
		typedef generic_iterator const_iterator;
		
		private:

		const_iterator *pbegin_, *pend_;

    size_type first_, last_; // range of indices covered.
		
		void compute_extremes(void)
		{
		 using std::min;
		 using std::max;
		 
		 first_ = std::numeric_limits<size_type>::max(); 
		 last_ = 0;
		 for(const_iterator it = *pbegin_; it != *pend_; ++it){
		   first_ = min(first_, *it);
			 last_ = max(last_, *it);
		 }		 
		}
				
		public:
				
		size_type size(void)const
		{ return  (*pend_ - *pbegin_); }
		
		size_type first(void)const { return first_; }
		
		size_type last(void)const { return last_; }
		
    size_type operator[](size_type i) const 
		{
		  assert(i < size());
      return (*pbegin_)[i];
    }

    size_type index(size_type i)const { return  operator[](i); }
		
    const_iterator& begin(void)const
		{ return *pbegin_; }
		
    const_iterator& end(void)const
		{ return *pend_; }		
		  
		void copy(const basic_index_ref& other)
		{
		 if (NULL != pbegin_){
		   delete pbegin_;
			 pbegin_ = NULL;
		 }
		 if (NULL != other.pbegin_)
		   pbegin_ = other.pbegin_->clone();
		 if (NULL != pend_){
		   delete pend_;
			 pend_ = NULL;
		 }
		 if (NULL != other.pend_)
		   pend_ = other.pend_->clone();		 
		 first_ = other.first_;
		 last_ = other.last_;
		}	 
		 
		basic_index_ref& operator=(const basic_index_ref& other)
		{ copy(other); return *this; }
		
    basic_index_ref()
		  : pbegin_(NULL), pend_(NULL), first_(0), last_(std::numeric_limits<size_type>::max())
		 { }

    basic_index_ref(const basic_index_ref& other)
		  : pbegin_(NULL), pend_(NULL), first_(0), last_(std::numeric_limits<size_type>::max())
    { copy(other); } 

    template <class IT>
		explicit
    basic_index_ref(const IT &itb, const IT &ite)
		  :pbegin_(new generic_iterator(itb)), pend_(new generic_iterator(ite))
	   { compute_extremes(); }

    template <class U>
		explicit
    basic_index_ref(const U &u)
		  :pbegin_(new generic_iterator(u.begin())), 
			 pend_(new generic_iterator(u.end()))
	   { compute_extremes(); }
		 
		~basic_index_ref(void)
		{ 
		 if (NULL != pbegin_)
		   delete pbegin_;
		 pbegin_ = NULL;
		 if (NULL != pend_)
		   delete pend_;
		 pend_ = NULL;	 
		}
		
  };

  #if 0
  class sub_index_ref {

    public:
		
    typedef basic_index_ref base_type;
    typedef typename base_type::const_iterator const_iterator;

    mutable base_type ind;

    size_type size(void) const { return ind.size(); }

    size_type index(size_type i) const { return ind[i]; }
  
    const_iterator  begin(void) const { return  ind.begin(); }
    const_iterator    end(void) const { return  ind.end();   }

		void copy(const sub_index_ref& other)
		{	ind = other.ind; }
		
    ~sub_index_ref()
    { }	
			
    sub_index_ref()
		{ }
		
		template <class IT>
    sub_index_ref(const IT& itb, const IT& ite)
      : ind(itb,ite)
		{ }
		
		template <class U>
    sub_index_ref(const U &u)
      : ind(u)
    { }
		
    sub_index_ref(const sub_index_ref &si) 
    { copy(si); }
		
    sub_index_ref &operator =(const sub_index_ref &other) 
		{
      copy(other);
      return *this;
    }
		
  };

	
  class unsorted_sub_index_ref : public basic_index_ref {
	
	  public:
		
		typedef typename basic_index_ref::const_iterator const_iterator;
		
		template <class IT>
    unsorted_sub_index_ref(const IT& itb, const IT& ite)
      : basic_index_ref(itb, ite) 
		{ }
		
		template <class U>
    unsorted_sub_index_ref(const U &u)
      : basic_index_ref(u) 
		{ }
		
    unsorted_sub_index_ref(void) 
		{ }
		
    unsorted_sub_index_ref(const unsorted_sub_index_ref &si) 
		  : basic_index_ref(si) 
		{ }
		
    unsorted_sub_index &operator =(const unsorted_sub_index &si)
    { basic_index_ref::operator =(si); return *this; }
    
  };
	#endif
	
  inline std::ostream &operator << (std::ostream &o, const basic_index_ref &si) { 
    o << "basic_index_ref(";
    if (si.size() != 0) o << si.index(0);
    for (size_type i = 1; i < si.size(); ++i) o << ", " << si.index(i);
    o << ")";
    return o;
  }

#if 0 // signature for appropriate tab_ref_index_ref_with_origin constructor:
vect_begin(const_cast<V&>(v)),
					    si.begin(), si.end()),
      origin(linalg_origin(const_cast<V&>(v)))
#endif
  // instantiations for "sub_vector_type" and "return_type" with sub-vector:
	template <typename PT>
  struct svrt_ir<PT, basic_index_ref, abstract_dense> {
    typedef typename std::iterator_traits<PT>::value_type V;
    typedef typename select_ref<typename linalg_traits<V>::const_iterator,
                                typename linalg_traits<V>::iterator, PT>::ref_type iterator;
    typedef tab_ref_index_ref_with_origin<iterator,
      typename basic_index_ref::const_iterator, V> vector_type;
  }; 

  template <typename PT>
  struct svrt_ir<PT, basic_index_ref, abstract_skyline> {
    typedef sparse_sub_vector<PT, basic_index_ref > vector_type;
  };

  template<> struct index_is_sorted<basic_index_ref>
  {  typedef linalg_false bool_type; };

#if 0
	// _try_ moving from gmm/gmm_interface.h:
template <typename IT, typename ITINDEX, typename V>	
	tab_ref_index_ref_with_origin<IT, ITINDEX, V>::	
	tab_ref_index_ref_with_origin(const V &v, const basic_index_ref &si)
    : dal::tab_ref_index_ref<IT, ITINDEX>(vect_begin(const_cast<V&>(v)),
					  si.begin(), si.end()),
            origin(linalg_origin(const_cast<V&>(v))) {}

template <typename IT, typename ITINDEX, typename V>
	tab_ref_index_ref_with_origin<IT, ITINDEX, V>::	
  tab_ref_index_ref_with_origin(V &v, const basic_index_ref &si)
    : dal::tab_ref_index_ref<IT, ITINDEX>(vect_begin(const_cast<V&>(v)),
					  si.begin(), si.end()),
	          origin(linalg_origin(const_cast<V&>(v))) {}		
#endif
		
} // namespace gmm 
 
 
 
/**
 * @ingroup arbitrary_precision
 * @brief Extension classes and methods for namespace gmm.
 *
 * The namespace "gmm_ext" is envisioned as a small group of extension methods to "gmm".
 * These methods are implemented using "gmm" style, rather than "numericalFunctor" style of interface.
 * At the point this becomes a _large_ group of extensions, we should consider switching to another
 *  matrix/vector template library.
 *
 *
 * Class list:
 *   "basic_index_ref": an unsorted sub-index that stores only the iterators.
 *   "gmm::linalg_traits<std::dequeue<T>>" : only declared when #include \<deque\> _prior_
 *
 * Current list of methods: 
 *   "QRP_factor": QR reduction with complete pivoting 
 *   "SVD_jacobi": Jacobi method for singular-value-decomposition(SVD)
 *   "SVD_QR": High relative accuracy QR method for SVD
 *
 */
 
namespace gmm_ext{

template <typename U>
  inline void fill_random(U& u);
  
template <typename U>
  void fill_random(U& u, gmm::abstract_matrix); 
  
template <typename U>
  void fill_random(U& u, gmm::abstract_vector);    


template <typename U, typename P>
  inline void permute_rows(const P& p, U& u);
	
template <typename U, typename P>
  void permute_rows(const P& p, U& u, gmm::abstract_matrix);	
	
template <typename U, typename P>
  void permute_rows(const P& p, U& u, gmm::abstract_vector);		
	
template <typename U, typename P>
  inline void permute_cols(const P& p, U& u);	

template <typename U, typename P>
  void permute_cols(const P& p, U& u, gmm::abstract_matrix);	
	
template <typename U, typename P>
  void permute_cols(const P& p, U& u, gmm::abstract_vector);	
	
		
// compute the QR factorization of A in Q and R, P: column permutations:
//   A [[P]] = Q R, where 
//      [[P]]_{i j} = \left{ 
//            \begin{array}{l}
//               1: j = P[i] //
//               0: j \ne P[i] 
//             \end{array} \right.
// 

// _all_ of these methods do the weird "gmm" const_cast thing on the output arguments:

// A = Q R P^T: QR with _complete_ pivoting, initial row-sorting absorbed into "Q", column-pivoting as "P"               
template <typename MAT1, typename MAT2, typename MAT3>
  void qrp_factor(const MAT1 &A, const MAT2 &Q_, const MAT3 &R_, std::vector<size_t>& P);

// A = U S V^{\daggar}
template <typename MAT1, typename MAT2, typename VEC3, typename MAT4>
  void svd_qr(const MAT1 &A, const MAT2 &U_, const VEC3 &S_, const MAT4 &V_);

// supporting utility methods for svd_QR:

// bi-diagonal decomposition:
//  A = U [D, F] V^{\daggar}, where D and F are diagonal and super-diagonal components of the decomposition:
template <typename MAT1, typename MAT2, typename VEC3, typename VEC4, typename MAT5>
  void bidiag_decomp( const MAT1 &A, const MAT2 &U_, const VEC3 &D_, const VEC4 &F_, const MAT5& V_);

template <typename VEC1, typename MAT2>
  void diag(const VEC1& D, const MAT2& A_);

template <typename VEC1, typename VEC2, typename MAT3>
  void diag(const VEC1& D, const VEC2& F, const MAT3& A_);  

/**
 * @ingroup arbitrary_precision
 * @brief Utility methods to support SVD and QR algorithms in namespace gmm_ext.
 *
 * For the moment, due to what would be the required complexity of the template class arguments, do _not_ put these in a static class,
 *   but temporarily just put them in their own namespace:  gmm_ext::svd_qr_util
 */
 
namespace svd_qr_util{

  // |f[i]| < |d[i] * d[i+1]| *eps -> 0 
  template <typename VEC1, typename VEC2>
    void chop_small_elements( const VEC1 &d, const VEC2 &f );

  // 2 x 2 SVD, by-hand:  
  template <typename VEC1, typename VEC2, typename MAT3, typename MAT4>
    void svd_2x2(const VEC1 &d, const VEC2 &f, const MAT3 &U, const MAT4 &V);

  // chase out zero _on_ d[k_0]:
  template <typename VEC1, typename VEC2, typename MAT3>
    void chase_out_intermediate_zero(const VEC1 &d, const VEC2 &f, const MAT3 &U, size_t k0);

  // chase out zero at end of diagonal d:
  template <typename VEC1, typename VEC2, typename MAT3>
    void chase_out_trailing_zero(const VEC1 &d, const VEC2 &f, const MAT3 &V);

  // trailing eigenvalue of a bidiagonal matrix: D: diagonal, F: super-diagonal.
  //   (note: for the moment, this routine implicitely assumes _real_ D and F; _probably_ OK for complex).
  template <class VEC1, class VEC2>
    typename gmm::number_traits<typename gmm::linalg_traits<VEC1>::value_type>::magnitude_type
  trailing_eigenvalue (const VEC1 &D, const VEC2 &F);

  // create Schur rotation to orthogonalize the columns of the _real_ bidiagonal 2 x 2 matrix [ D; F ]:
  // (note: this is not used for _complex_ matrices above, but eventually it should be extended so it works correctly for complex matrices).
  template <class T1, class T2>
    void schur_rotation (T1 d0, T1 f0, T1 d1, T2 &c, T2 &s);

  // one step of the QR SVD of the bidiagonal matrix:
  // (implicitely this assumes that D and F are _real_)
  template <typename VEC1, typename VEC2, typename MAT3, typename MAT4>
    void qr_step(const VEC1 &D_, const VEC2 &F_, const MAT3 &U_, const MAT4 &V_);

} // namespace svd_qr_util


// null_space of a matrix: return matrix V of column-vectors such that A V(:,k) = 0

// parameter NV for the null space column vectors is either an allocated vector, or an unallocated matrix.

template <typename MAT1, typename MAT2>
inline void null_space( const MAT1& A, const MAT2& NV_,
                        typename gmm::number_traits<typename gmm::linalg_traits<MAT1>::value_type>::magnitude_type nullRelativeTol 
                        = number::epsilon<typename gmm::number_traits<typename gmm::linalg_traits<MAT1>::value_type>::magnitude_type>() );

// NV matrix version:
// For this version the dimension of the null space is unknown => NV cannot be pre-sized (convention: NV shall be 0 x 0):
template <typename MAT1, typename MAT2>
  void null_space( const MAT1& A, const MAT2& NV_,
                   typename gmm::number_traits<typename gmm::linalg_traits<MAT1>::value_type>::magnitude_type nullRelativeTol,
									 gmm::abstract_matrix );

// NV vector version:
// For this version the dimension of the null space is 1, A is M x N matrix => NV shall be allocated with length min(M,N).
//   If more (or less) than one null vector is found, an exception is thrown.
template <typename MAT1, typename VEC2>
  void null_space( const MAT1& A, const VEC2& NV_,
                   typename gmm::number_traits<typename gmm::linalg_traits<MAT1>::value_type>::magnitude_type nullRelativeTol,
									 gmm::abstract_vector );
  
#if 0
// NOT yet quite correct for arbitrary complex matrices.

template <typename MAT1, typename MAT2, typename VEC3, typename MAT4>
  void svd_jacobi(const MAT1 &A, const MAT2 &U_, const VEC3 &S_, const MAT4 &V_);
	
// Complex Jacobi rotation utilities:
template <class T>
void jacobi_rotation( const typename gmm::number_traits<T>::magnitude_type& a_jj,
                      const typename gmm::number_traits<T>::magnitude_type& a_kk,
                      const T& a_jk,
                      typename gmm::number_traits<T>::magnitude_type& c,
                      typename gmm::number_traits<T>::magnitude_type& s,
                      T& e );
                      
// Apply Q* v
template <typename T> inline
void apply_jacobi_rotation_left(const typename gmm::number_traits<T>::magnitude_type &c,
                                const typename gmm::number_traits<T>::magnitude_type &s,
                                const T &e,
                                T& x, T& y);

// Apply v^T Q
template <typename T> inline
void apply_jacobi_rotation_right(const typename gmm::number_traits<T>::magnitude_type &c,
                                 const typename gmm::number_traits<T>::magnitude_type &s,
                                 const T &e,
                                 T& x, T& y);

template <typename MAT, typename T>
void apply_jacobi_rotation_left(const typename gmm::number_traits<T>::magnitude_type &c,
                                const typename gmm::number_traits<T>::magnitude_type &s,
                                const T &e,
                                size_t i, size_t k, const MAT &AA); 

template <typename MAT, typename T>
void apply_jacobi_rotation_right (const typename gmm::number_traits<T>::magnitude_type &c,
                                  const typename gmm::number_traits<T>::magnitude_type &s,
                                  const T &e,
                                  size_t i, size_t k, const MAT &AA); 
#endif
 

// modified LAPACK style Householder reflections:
// The following are variants from the "gmm" methods; they should _not_ conflict because they have differing signatures from "gmm".
// These methods use the LAPACK style complex scale-factor tau, and do not use the _explicit_ temporary "W".
// Differing from LAPACK, the reflection vector is _not_ in packed-format, and it is expected that V(0) = 1.0.
// (at present, the explicit temporary is considered to be an over-optimization)

template <typename VEC1> void house_vector(const VEC1 &V_,
    const typename gmm::linalg_traits<VEC1>::value_type &tau_, 
    const typename gmm::number_traits<typename gmm::linalg_traits<VEC1>::value_type>::magnitude_type &beta_ );
    
// multiply A to the left by the reflector stored in V, tau.
template <typename MAT1, typename VEC2> inline
  void row_house_update(const MAT1 &A_, 
    const VEC2 &V, const typename gmm::linalg_traits<MAT1>::value_type &tau);
    
// multiply A to the right by the reflector stored in V, tau. 
template <typename MAT1, typename VEC2> inline
  void col_house_update(const MAT1 &A_, const VEC2 &V, const typename gmm::linalg_traits<MAT1>::value_type &tau);          
    
} // namespace gmm_ext

// -------------------- omit, for the moment: ------------------------
#if 0
#include "gmm_ext_householder.h"
#include "gmm_ext_template.h"
#endif
// -------------------------------------------------------------------

#endif // __gmm_ext__h 
