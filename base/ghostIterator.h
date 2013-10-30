#if !defined(__ghostIterator__h)
#define __ghostIterator__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/ghostIterator.h,v $

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


// A helper class to allow such things as sorting two vectors, driven by the sort of the first, with the STL sort routine.

namespace linalg{


template <class IT_A, class IT_B>
class ghostIterator{
  private:
    IT_A itA_;
    IT_B itB_; 
  public:
  typedef ghostIterator<IT_A, IT_B> this_type;

  class A_ref0; // forward
  class A_ref{
    public:
      // mutable references to the actual values:
      typename std::iterator_traits<IT_A>::reference A(void);
      typename std::iterator_traits<IT_B>::reference B(void);
      // and non-mutable references:
      /* const */ typename std::iterator_traits<IT_A>::reference A(void)const;
      /* const */ typename std::iterator_traits<IT_B>::reference B(void)const;
 
      void swap(A_ref other);
       
      A_ref& operator=(const A_ref& other); // assign from reference
      A_ref& operator=(const A_ref0& other); 
                  
      bool operator==(const A_ref& other)const;
      
      bool operator!=(const A_ref& other)const;
           
      bool operator < (const A_ref& other)const;
      
      bool operator > (const A_ref& other)const;

      bool operator==(const A_ref0& other)const;
      
      bool operator!=(const A_ref0& other)const;
           
      bool operator < (const A_ref0& other)const;
      
      bool operator > (const A_ref0& other)const;
                  
      A_ref(ghostIterator<IT_A, IT_B>& parent); // attached to iterators

    protected:
      ghostIterator<IT_A, IT_B>* pParent;

  };  
  friend class A_ref;

  class A_ref0{
    public:
      // mutable references to the actual values:
      typename std::iterator_traits<IT_A>::reference A(void);
      typename std::iterator_traits<IT_B>::reference B(void);
#if 0
      // and non-mutable references:
      const typename std::iterator_traits<IT_A>::reference A(void)const;
      const typename std::iterator_traits<IT_B>::reference B(void)const;
#endif
      const typename std::iterator_traits<IT_A>::value_type A(void)const;
      const typename std::iterator_traits<IT_B>::value_type B(void)const;

      void swap(A_ref0 other);

      // assign from A_ref reference       
      A_ref0& operator=(const A_ref& other);
      A_ref0& operator=(const A_ref0& other);
            
      bool operator==(const A_ref0& other)const;
      
      bool operator!=(const A_ref0& other)const;
           
      bool operator < (const A_ref0& other)const;
      
      bool operator > (const A_ref0& other)const;

      bool operator==(const A_ref& other)const;
      
      bool operator!=(const A_ref& other)const;
           
      bool operator < (const A_ref& other)const;
      
      bool operator > (const A_ref& other)const;
           
      A_ref0(const A_ref& other);                // local case: functions as (IT_A::value_type, IT_B::value_type)
      A_ref0(const A_ref0& other);               // local case: functions as (IT_A::value_type, IT_B::value_type)
      A_ref0(void);                              // local case: functions as (IT_A::value_type, IT_B::value_type)   

    protected:
      typename std::iterator_traits<IT_A>::value_type A_;
      typename std::iterator_traits<IT_B>::value_type B_;            
  };
  friend class A_ref0;
  
  // In order to allow std::iterator_traits to function:
  typedef typename std::iterator_traits< IT_A >::iterator_category iterator_category;
  typedef A_ref0        value_type;
  typedef typename std::iterator_traits< IT_A >::difference_type difference_type;
  typedef A_ref*           pointer;
  typedef A_ref         reference;
            
  A_ref operator*(void);

  
  bool operator==(const ghostIterator<IT_A, IT_B>& other)const;

  bool operator!=(const ghostIterator<IT_A, IT_B>& other)const;

  bool operator < (const ghostIterator<IT_A, IT_B>& other)const;
  
  bool operator > (const ghostIterator<IT_A, IT_B>& other)const;  
  
  ghostIterator<IT_A, IT_B>& operator+=(size_t n);
  
  ghostIterator<IT_A, IT_B> operator+(size_t n)const; 

  ghostIterator<IT_A, IT_B>& operator++(void);

  ghostIterator<IT_A, IT_B> operator++(int);

  ghostIterator<IT_A, IT_B>& operator-=(size_t n);

  ghostIterator<IT_A, IT_B> operator-(size_t n)const; 

  ghostIterator<IT_A, IT_B>& operator--(void);

  ghostIterator<IT_A, IT_B> operator--(int);

  difference_type operator-(const ghostIterator<IT_A, IT_B>& other)const;

  ghostIterator(IT_A itA, IT_B itB);

  ghostIterator(const ghostIterator<IT_A, IT_B>& other);

};

} // namespace linalg


#if 0 // STL headers can't _see_ this definition, for some reason:
namespace std{
template <class IT_A, class IT_B>
void swap( typename linalg::ghostIterator<IT_A, IT_B>::A_ref t1, typename linalg::ghostIterator<IT_A, IT_B>::A_ref t2 )
{
  t1.swap(t2);
}
} // namespace std
#endif

#if defined(__PGI) // STLport standard template library specific:
namespace std{
template <class IT_A, class IT_B>
void iter_swap( linalg::ghostIterator<IT_A, IT_B>& it1, linalg::ghostIterator<IT_A, IT_B>& it2 )
{
  (*it1).swap((*it2));
}
} // namespace std
#endif

#include "ghostIterator_template.h"

#endif // __ghostIterator__h
