#if !defined(__ghostIterator_template__h)
#define __ghostIterator_template__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/ghostIterator_template.h,v $

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
  inline typename std::iterator_traits<IT_A>::reference 
    ghostIterator<IT_A, IT_B>::A_ref::A(void)
  {
    return *(pParent->itA_);
  }  

  template <class IT_A, class IT_B>
  inline typename std::iterator_traits<IT_B>::reference 
    ghostIterator<IT_A, IT_B>::A_ref::B(void)
  {
    return *(pParent->itB_);
  }     

  template <class IT_A, class IT_B>
  inline /* const */ typename std::iterator_traits<IT_A>::reference 
    ghostIterator<IT_A, IT_B>::A_ref::A(void)const
  {
    return *(pParent->itA_);
  }  

  template <class IT_A, class IT_B>
  inline /* const */ typename std::iterator_traits<IT_B>::reference 
    ghostIterator<IT_A, IT_B>::A_ref::B(void)const
  {
    return *(pParent->itB_);
  } 

  template <class IT_A, class IT_B>
  inline void ghostIterator<IT_A, IT_B>::A_ref::swap(A_ref other)
  {
   std::swap(A(), other.A());
   std::swap(B(), other.B());
  }

   
  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::A_ref& ghostIterator<IT_A, IT_B>::A_ref::operator=(const A_ref& other)
  { 
    A() = other.A();
    B() = other.B();
   return *this;
  } 

  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::A_ref& ghostIterator<IT_A, IT_B>::A_ref::operator=(const A_ref0& other)
  { 
    A() = other.A();
    B() = other.B();
   return *this;
  } 
    
  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator==(const A_ref& other)const
  { 
    return A() == other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator!=(const A_ref& other)const
  { return !operator==(other); }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator<(const A_ref& other)const
  { 
    return A() < other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator>(const A_ref& other)const
  { 
    return A() > other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator==(const A_ref0& other)const
  { 
    return A() == other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator!=(const A_ref0& other)const
  { return !operator==(other); }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator<(const A_ref0& other)const
  { 
    return A() < other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref::operator>(const A_ref0& other)const
  { 
    return A() > other.A(); 
  }


  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::A_ref::A_ref(ghostIterator<IT_A, IT_B>& parent)
    : pParent(&parent)
  {}


  template <class IT_A, class IT_B>
  inline typename std::iterator_traits<IT_A>::reference 
    ghostIterator<IT_A, IT_B>::A_ref0::A(void)
  {
    return A_;
  }  

  template <class IT_A, class IT_B>
  inline typename std::iterator_traits<IT_B>::reference 
    ghostIterator<IT_A, IT_B>::A_ref0::B(void)
  {
    return B_;
  }     

#if 0
  template <class IT_A, class IT_B>
  inline const typename std::iterator_traits<IT_A>::reference 
    ghostIterator<IT_A, IT_B>::A_ref0::A(void)const
  {
    return A_;
  }  

  template <class IT_A, class IT_B>
  inline const typename std::iterator_traits<IT_B>::reference 
    ghostIterator<IT_A, IT_B>::A_ref0::B(void)const
  {
    return B_;
  } 
#endif

  template <class IT_A, class IT_B>
  inline const typename std::iterator_traits<IT_A>::value_type 
    ghostIterator<IT_A, IT_B>::A_ref0::A(void)const
  {
    return A_;
  }  

  template <class IT_A, class IT_B>
  inline const typename std::iterator_traits<IT_B>::value_type 
    ghostIterator<IT_A, IT_B>::A_ref0::B(void)const
  {
    return B_;
  } 
  
  template <class IT_A, class IT_B>
  inline void ghostIterator<IT_A, IT_B>::A_ref0::swap(A_ref0 other)
  {
   std::swap(A(), other.A());
   std::swap(B(), other.B());
  }

   
  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::A_ref0& ghostIterator<IT_A, IT_B>::A_ref0::operator=(const A_ref& other)
  { 
    A() = other.A();
    B() = other.B();
   return *this;
  } 

  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::A_ref0& ghostIterator<IT_A, IT_B>::A_ref0::operator=(const A_ref0& other)
  { 
    A() = other.A();
    B() = other.B();
   return *this;
  } 
    
  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator==(const A_ref0& other)const
  { 
    return A() == other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator!=(const A_ref0& other)const
  { return !operator==(other); }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator<(const A_ref0& other)const
  { 
    return A() < other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator>(const A_ref0& other)const
  { 
    return A() > other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator==(const A_ref& other)const
  { 
    return A() == other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator!=(const A_ref& other)const
  { return !operator==(other); }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator<(const A_ref& other)const
  { 
    return A() < other.A(); 
  }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::A_ref0::operator>(const A_ref& other)const
  { 
    return A() > other.A(); 
  }




  // local case: functions as (IT_A::value_type, IT_B::value_type)
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::A_ref0::A_ref0(const A_ref0& other)
     : A_(other.A()),
       B_(other.B()) 
  {}

  // local case: functions as (IT_A::value_type, IT_B::value_type)
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::A_ref0::A_ref0(const A_ref& other)
     : A_(other.A()),
       B_(other.B()) 
  {}

  // local case: functions as (IT_A::value_type, IT_B::value_type)      
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::A_ref0::A_ref0(void)
     : A_(number::zero<typename std::iterator_traits<IT_A>::value_type>()),
       B_(number::zero<typename std::iterator_traits<IT_B>::value_type>()) 
  {}  


  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::A_ref ghostIterator<IT_A, IT_B>::operator*(void)
  { return A_ref(*this); }
  
  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::operator==(const ghostIterator<IT_A, IT_B>& other)const
  { return itA_==other.itA_; }
  
  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::operator!=(const ghostIterator<IT_A, IT_B>& other)const
  { return !operator==(other); }

  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::operator < (const ghostIterator<IT_A, IT_B>& other)const
  { return itA_ < other.itA_; }


  template <class IT_A, class IT_B>
  inline bool ghostIterator<IT_A, IT_B>::operator > (const ghostIterator<IT_A, IT_B>& other)const
  { return itA_ > other.itA_; }
      
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>& ghostIterator<IT_A, IT_B>::operator+=(size_t n)
  {
   itA_ += n;
   itB_ += n;
   return *this;
  }
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B> ghostIterator<IT_A, IT_B>::operator+(size_t n)const 
  {
   ghostIterator<IT_A, IT_B> val(*this);
   val.itA_ += n;
   val.itB_ += n;
   return val;
  }
  
  template <class IT_A, class IT_B> 
  inline ghostIterator<IT_A, IT_B>& ghostIterator<IT_A, IT_B>::operator++(void)
  {
   ++itA_;
   ++itB_;
   return *this;
  }
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B> ghostIterator<IT_A, IT_B>::operator++(int)
  {
   ghostIterator<IT_A, IT_B> val(*this);
   operator++();
   return val;
  }
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>& ghostIterator<IT_A, IT_B>::operator-=(size_t n)
  {
   itA_ -= n;
   itB_ -= n;
   return *this;
  }
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B> ghostIterator<IT_A, IT_B>::operator-(size_t n)const 
  {
   ghostIterator<IT_A, IT_B> val(*this);
   val.itA_ -= n;
   val.itB_ -= n;
   return val;
  }
  
  template <class IT_A, class IT_B> 
  inline ghostIterator<IT_A, IT_B>& ghostIterator<IT_A, IT_B>::operator--(void)
  {
   --itA_;
   --itB_;
   return *this;
  }
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B> ghostIterator<IT_A, IT_B>::operator--(int)
  {
   ghostIterator<IT_A, IT_B> val(*this);
   operator--();
   return val;
  }
  
  template <class IT_A, class IT_B>
  inline typename ghostIterator<IT_A, IT_B>::difference_type ghostIterator<IT_A, IT_B>::operator-(const ghostIterator<IT_A, IT_B>& other)const
  { return itA_ - other.itA_; } 
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::ghostIterator(IT_A itA, IT_B itB)
    : itA_(itA), itB_(itB)
  {}
  
  template <class IT_A, class IT_B>
  inline ghostIterator<IT_A, IT_B>::ghostIterator(const ghostIterator<IT_A, IT_B>& other)
    : itA_(other.itA_), itB_(other.itB_)
  {}
 

} // namespace linalg

#endif // __ghostIterator_template__h
