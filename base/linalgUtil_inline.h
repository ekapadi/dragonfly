#if !defined(__linalgUtil_inline__h)
#define __linalgUtil_inline__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/linalgUtil_inline.h,v $

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


namespace linalg{

inline std::vector<size_t> inverse_row_major_index(size_t n, const std::vector<size_t>& shape)
{
  std::vector<size_t> dest;
  inverse_row_major_index(n, shape, dest);
  return dest;
}

} // namespace linalg


namespace commUtil{

inline bool writeBinary(commUtil::abstractCommHandle* fp, const linalg::CSYS_KIND& e) 
{ return writeBinary(fp,static_cast<long>(e)); }

inline bool readBinary(commUtil::abstractCommHandle* fp, linalg::CSYS_KIND& e)
{
 bool val(true);

 long n;
 val = (val && readBinary(fp,n));
 e = static_cast<linalg::CSYS_KIND>(n);
 return val;
}

inline size_t binarySize(const linalg::CSYS_KIND& e) { return binarySize(static_cast<long>(e)); }


inline void write(std::ostream& os, const linalg::CSYS_KIND& e)
{
 std::string s;
 switch (e){
   case linalg::NO_CSYS:
	   s = "NO_CSYS";
	 break;
	 case linalg::CARTESIAN:
	   s = "CARTESIAN";
	 break;
	 case linalg::CYLINDRICAL:
	   s = "CYLINDRICAL";
	 break;
	 case linalg::SPHERICAL:
	   s = "SPHERICAL";
	 break;
	 default:
	   throw std::runtime_error("write(std::ostream& os, const CSYS_KIND& e): invalid CSYS_KIND");
//	 break;
 };
 os<<s;
}

inline void read(std::istream& is, linalg::CSYS_KIND& e)
{
 std::string s;
 is >> s;
 if (is){
   if (s == "NO_CSYS")
     e = linalg::NO_CSYS;
	 else
	 if (s == "CARTESIAN")
	   e = linalg::CARTESIAN;
	 else
	 if (s == "CYLINDRICAL")
	   e = linalg::CYLINDRICAL;
	 else
	 if (s == "SPHERICAL")
	   e = linalg::SPHERICAL;
	 else
	   throw std::runtime_error("read(std::istream& is, CSYS_KIND& e): invalid CSYS_KIND");	    	 
 }
}

inline std::ostream& operator<<(std::ostream& os, const linalg::CSYS_KIND& e) { write(os,e); return os; }

inline std::istream& operator>>(std::istream& is, linalg::CSYS_KIND& e) { read(is,e); return is; }

#if defined(__specialize_POD_binary__)	
// specializations for contiguous POD types:
template<>
inline bool writeBinary<double>(abstractCommHandle *fp, const std::vector<double>& V)
{
   bool rVal = true;
   
   // write size:
   size_t nSize = V.size();
   rVal = (rVal && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     rVal = (rVal && (1 == write(&(*V.begin()), sizeof(double)*nSize, 1, fp)));
        
   return rVal;
}

template<>
inline bool readBinary<double>(abstractCommHandle *fp, std::vector<double>& V )
{
   bool rVal = true;
   
   // read size
   size_t nSize;
   rVal = (rVal && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   V.resize( nSize ); 

   if (0 < nSize)
     rVal = (rVal && (1 == read(&(*V.begin()), sizeof(double)*nSize, 1, fp)));

   return rVal;
}

template<>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const std::vector<std::complex<double> >& V)
{
   bool rVal = true;
   
   // write size:
   size_t nSize = V.size();
   rVal = (rVal && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     rVal = (rVal && (1 == write(&(*V.begin()), sizeof(std::complex<double>)*nSize, 1, fp)));
        
   return rVal;
}

template<>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, std::vector<std::complex<double> >& V )
{
   bool rVal = true;
   
   // read size
   size_t nSize;
   rVal = (rVal && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   V.resize( nSize ); 

   if (0 < nSize)
     rVal = (rVal && (1 == read(&(*V.begin()), sizeof(std::complex<double>)*nSize, 1, fp)));

   return rVal;
}

template<>
inline bool writeBinary<long>(abstractCommHandle *fp, const std::vector<long>& V)
{
   bool rVal = true;
   
   // write size:
   size_t nSize = V.size();
   rVal = (rVal && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     rVal = (rVal && (1 == write(&(*V.begin()), sizeof(long)*nSize, 1, fp)));
        
   return rVal;
}

template<>
inline bool readBinary<long>(abstractCommHandle *fp, std::vector<long>& V )
{
   bool rVal = true;
   
   // read size
   size_t nSize;
   rVal = (rVal && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   V.resize( nSize ); 

   if (0 < nSize)
     rVal = (rVal && (1 == read(&(*V.begin()), sizeof(long)*nSize, 1, fp)));

   return rVal;
}

template<>
inline bool writeBinary<size_t>(abstractCommHandle *fp, const std::vector<size_t>& V)
{
   bool rVal = true;
   
   // write size:
   size_t nSize = V.size();
   rVal = (rVal && (1 == write(&nSize, sizeof(size_t), 1, fp)));

   if (0 < nSize)
     rVal = (rVal && (1 == write(&(*V.begin()), sizeof(size_t)*nSize, 1, fp)));
        
   return rVal;
}

template<>
inline bool readBinary<size_t>(abstractCommHandle *fp, std::vector<size_t>& V )
{
   bool rVal = true;
   
   // read size
   size_t nSize;
   rVal = (rVal && (1 == read(&nSize, sizeof(size_t), 1, fp)));
   V.resize( nSize ); 

   if (0 < nSize)
     rVal = (rVal && (1 == read(&(*V.begin()), sizeof(size_t)*nSize, 1, fp)));

   return rVal;
}
#endif


#if defined(__specialize_POD_binary__)
// specializations for contiguous POD types:
template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::dense_matrix<double>& M)
{
   bool rVal = true;
   
   // write size:
   const size_t 
     nr(gmm::mat_nrows(M)),
     nc(gmm::mat_ncols(M));
   rVal = (rVal && (1 == write(&nr, sizeof(size_t), 1, fp)));
   rVal = (rVal && (1 == write(&nc, sizeof(size_t), 1, fp)));

   if (0 < nr*nc)
     rVal = (rVal && (1 == write(&(M(0,0)), sizeof(double)*nr*nc, 1, fp)));
        
   return rVal;
}

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::dense_matrix<double>& M )
{
   bool rVal = true;
   
   // read size
   size_t nr,nc;
   rVal = (rVal && (1 == read(&nr, sizeof(size_t), 1, fp)));
   rVal = (rVal && (1 == read(&nc, sizeof(size_t), 1, fp)));

   if (rVal)
     gmm::resize(M, nr, nc); 

   if (0 < nr*nc)
     rVal = (rVal && (1 == read(&(M(0,0)), sizeof(double)*nr*nc, 1, fp)));

   return rVal;
}

template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::dense_matrix<std::complex<double> >& M)
{
   bool rVal = true;
   
   // write size:
   const size_t 
     nr(gmm::mat_nrows(M)),
     nc(gmm::mat_ncols(M));
   rVal = (rVal && (1 == write(&nr, sizeof(size_t), 1, fp)));
   rVal = (rVal && (1 == write(&nc, sizeof(size_t), 1, fp)));

   if (0 < nr*nc)
     rVal = (rVal && (1 == write(&(M(0,0)), sizeof(std::complex<double>)*nr*nc, 1, fp)));
        
   return rVal;
}

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::dense_matrix<std::complex<double> >& M )
{
   bool rVal = true;
   
   // read size
   size_t nr,nc;
   rVal = (rVal && (1 == read(&nr, sizeof(size_t), 1, fp)));
   rVal = (rVal && (1 == read(&nc, sizeof(size_t), 1, fp)));

   if (rVal)
     gmm::resize(M, nr, nc); 

   if (0 < nr*nc)
     rVal = (rVal && (1 == read(&(M(0,0)), sizeof(std::complex<double>)*nr*nc, 1, fp)));

   return rVal;
}

template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::slvector<double>& V)
{
 bool status(true);
 const size_t 
   size_(gmm::vect_size(V)),
   data_size_(gmm::nnz(V)),
   shift_(V.first()); 
   
 status = ( status && (1 == write(&size_, sizeof(size_t), 1, fp) ) ); 
 status = ( status && (1 == write(&data_size_, sizeof(size_t), 1, fp) ) );
 status = ( status && (1 == write(&shift_, sizeof(size_t), 1, fp) ) ); 
  
 // write the contiguous section:
 status = (status && (1 == write(&(*V.data_begin()), data_size_*sizeof(double), 1, fp))); 

 return status;    
}

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::slvector<double>& V )
{
 bool status(true);
 size_t 
   size_(0),
   data_size_(0),
   shift_(0); 
   
 status = ( status && (1 == read(&size_, sizeof(size_t), 1, fp) ) ); 
 status = ( status && (1 == read(&data_size_, sizeof(size_t), 1, fp) ) );
 status = ( status && (1 == read(&shift_, sizeof(size_t), 1, fp) ) ); 
  
 if (status){ 
   // This is the only efficient (hopefully) and _legal_ way to set the "shift" (apart from re-writing gmm code).
   // (here it is assumed that "swap" swaps the data pointers).
   gmm::slvector<double> u_(size_, data_size_, shift_);

   // read the contiguous section:
   status = (status && (1 == read(&(*u_.data_begin()), data_size_*sizeof(double), 1, fp))); 
     
   if (status)
     std::swap(V, u_);   	 
 }
 
 return status;    
}

template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::slvector<std::complex<double> >& V)
{
 bool status(true);
 const size_t 
   size_(gmm::vect_size(V)),
   data_size_(gmm::nnz(V)),
   shift_(V.first()); 
   
 status = ( status && (1 == write(&size_, sizeof(size_t), 1, fp) ) ); 
 status = ( status && (1 == write(&data_size_, sizeof(size_t), 1, fp) ) );
 status = ( status && (1 == write(&shift_, sizeof(size_t), 1, fp) ) ); 
  
 // write the contiguous section:
 status = (status && (1 == write(&(*V.data_begin()), data_size_*sizeof(std::complex<double>), 1, fp))); 

 return status;    
}

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::slvector<std::complex<double> >& V )
{
 bool status(true);
 size_t 
   size_(0),
   data_size_(0),
   shift_(0); 
   
 status = ( status && (1 == read(&size_, sizeof(size_t), 1, fp) ) ); 
 status = ( status && (1 == read(&data_size_, sizeof(size_t), 1, fp) ) );
 status = ( status && (1 == read(&shift_, sizeof(size_t), 1, fp) ) ); 
  
 if (status){ 
   // This is the only efficient (hopefully) and _legal_ way to set the "shift" (apart from re-writing gmm code).
   // (here it is assumed that "swap" swaps the data pointers).
   gmm::slvector<std::complex<double> > u_(size_, data_size_, shift_);

   // read the contiguous section:
   status = (status && (1 == read(&(*u_.data_begin()), data_size_*sizeof(std::complex<double>), 1, fp))); 
     
   if (status)
     std::swap(V, u_);   	 
 }
 
 return status;    
}

template <>
inline bool writeBinary<double>(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<double> >& M)
{
 bool status(true);
 
 typedef gmm::linalg_traits<gmm::row_matrix< gmm::slvector<double> > > traits;
   
 // write size:
 const size_t 
   nr(gmm::mat_nrows(M)),
   nc(gmm::mat_ncols(M));
 status = (status && (1 == write(&nr, sizeof(size_t), 1, fp)));
 status = (status && (1 == write(&nc, sizeof(size_t), 1, fp)));

 // iterate and write by row vector:
 for(traits::const_row_iterator 
       itRow = gmm::mat_row_const_begin(M), itRowEnd = gmm::mat_row_const_end(M);
     status && (itRow != itRowEnd);
     ++itRow)  // write the row slvector<double> itself:
   status = (status && writeBinary(fp,*gmm::linalg_traits<traits::const_sub_row_type>::origin(traits::row(itRow)) ));   

 return status;   
}

template <>
inline bool readBinary<double>(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<double> >& M )
{
 bool status(true);
 
 typedef gmm::linalg_traits<gmm::row_matrix< gmm::slvector<double> > > traits;
   
 // read size
 size_t nr,nc;
 status = (status && (1 == read(&nr, sizeof(size_t), 1, fp)));
 status = (status && (1 == read(&nc, sizeof(size_t), 1, fp)));

 if (status)
   gmm::resize(M, nr, nc); 

 // iterate and read by row vector:
 for(traits::row_iterator 
       itRow = gmm::mat_row_begin(M), itRowEnd = gmm::mat_row_end(M);
     status && (itRow != itRowEnd);
     ++itRow)  // read the row slvector<double> itself:
   status = (status && readBinary(fp,*const_cast<gmm::slvector<double>* >(gmm::linalg_traits<traits::sub_row_type>::origin(traits::row(itRow))) ));   

 return status;   
}

template <>
inline bool writeBinary<std::complex<double> >(abstractCommHandle *fp, const gmm::row_matrix< gmm::slvector<std::complex<double> > >& M)
{
 bool status(true);
 
 typedef gmm::linalg_traits<gmm::row_matrix< gmm::slvector<std::complex<double> > > > traits;

 // write size:
 const size_t 
   nr(gmm::mat_nrows(M)),
   nc(gmm::mat_ncols(M));
 status = (status && (1 == write(&nr, sizeof(size_t), 1, fp)));
 status = (status && (1 == write(&nc, sizeof(size_t), 1, fp)));
 
 // iterate and write by row vector:
 for(traits::const_row_iterator 
       itRow = gmm::mat_row_const_begin(M), itRowEnd = gmm::mat_row_const_end(M);
     status && (itRow != itRowEnd);
     ++itRow)  // write the row slvector<double> itself:
   status = (status && writeBinary(fp,*gmm::linalg_traits<traits::const_sub_row_type>::origin(traits::row(itRow)) ));   

 return status;   
}

template <>
inline bool readBinary<std::complex<double> >(abstractCommHandle *fp, gmm::row_matrix< gmm::slvector<std::complex<double> > >& M )
{
 bool status(true);
 
 typedef gmm::linalg_traits<gmm::row_matrix< gmm::slvector<std::complex<double> > > > traits;
   
 // read size
 size_t nr,nc;
 status = (status && (1 == read(&nr, sizeof(size_t), 1, fp)));
 status = (status && (1 == read(&nc, sizeof(size_t), 1, fp)));

 if (status)
   gmm::resize(M, nr, nc); 
 
 // iterate and read by row vector:
 for(traits::row_iterator 
       itRow = gmm::mat_row_begin(M), itRowEnd = gmm::mat_row_end(M);
     status && (itRow != itRowEnd);
     ++itRow)  // read the row slvector<double> itself:
   status = (status && readBinary(fp,*const_cast<gmm::slvector<std::complex<double> >* >(gmm::linalg_traits<traits::sub_row_type>::origin(traits::row(itRow))) ));   

 return status;   
}

#endif

								 
} // namespace commUtil


#endif // __linalgUtil_inline__h
