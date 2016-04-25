#if !defined(__commUtil_inline__h)
#define __commUtil_inline__h

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/commUtil_inline.h,v $

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


namespace commUtil{
      
inline void abstractCommHandle::buffer::alloc(size_t size)
{
 if (allocated())
   free();
 begin_ = new unsigned char[size + sizeof(header)];
 if (NULL == begin_)
   throw std::runtime_error("abstractCommHandle::buffer::alloc: allocation error");
 end_ = begin_ + (size + sizeof(header));
 set_data_size(end_ - begin_ - header_size()); // actually write the data_size into the buffer.
 ownsData_ = true;
 init();
}

inline void abstractCommHandle::buffer::attach(unsigned char *begin,
                                               unsigned char *end)
{
 if (allocated())
   free();
 if ((NULL == begin) || (NULL == end) || (end < begin) || (static_cast<size_t>(end - begin) < header_size()))
   throw std::runtime_error("abstractCommHandle::buffer::attach: invalid pointers");
 begin_ = begin;
 end_ = end;
 set_data_size(end_ - begin_ - header_size()); // actually write the data_size into the buffer.
 ownsData_ = false;
 init();
}

inline void abstractCommHandle::buffer::free(void)
{
 if (allocated()){
   if (ownsData_){
     delete[] begin_;
     ownsData_ = false;
   }
   begin_ = NULL;
   end_ = NULL;
   data_begin_ = NULL;
   data_end_ = NULL;
 }
}

inline void abstractCommHandle::buffer::init(size_t data_size_)
{
  set_data_size(data_size_);
  init();
}

inline void abstractCommHandle::buffer::init(void)
{
 if (data_size() > (end_ - begin_ - sizeof(header)))
   throw std::runtime_error("abstractCommHandle::buffer::init: invalid data size");
 data_begin_ = begin_ + sizeof(header);
 data_end_ = data_begin_ + data_size();
 current_ = data_begin_;
}

inline size_t abstractCommHandle::buffer::data_size(void)const  // doesn't include header size
{
 size_t val(0);
 if (allocated())
   val = reinterpret_cast<const header*>(begin_)->size_;
 return val;  
}

inline void abstractCommHandle::buffer::set_data_size(const size_t size)
{
 if (!allocated())
   throw std::runtime_error("abstractCommHandle::buffer::set_data_size: can't set data_size of unallocated buffer");         
 header &data_size_(*reinterpret_cast<header*>(begin_));
 data_size_.size_ = size;
} 

inline size_t abstractCommHandle::buffer::size(void)const // includes header size
{
 size_t val(0);
 if (allocated())
   val =  end_ - begin_;
 return val;
}

inline size_t abstractCommHandle::buffer::used(void)const // does not include header size
{
 size_t val(0);
 if (allocated())
   val =  current_ - data_begin_;
 return val;
}

inline size_t abstractCommHandle::buffer::header_size(void)const
{
 return sizeof(header);
}

inline bool abstractCommHandle::buffer::allocated(void)const
{
 return (NULL != begin_);
}      

inline abstractCommHandle::buffer::buffer(void)
  : ownsData_(false),
    begin_(NULL), data_begin_(NULL), current_(NULL), data_end_(NULL), end_(NULL)
{ }
    
inline size_t abstractCommHandle::readBufAvail(void) const
{
 size_t val(0);
 if (readBuf_.allocated())
   val = readBuf_.data_end_ - readBuf_.current_;
 return val;  
}

inline size_t abstractCommHandle::writeBufAvail(void) const
{
 size_t val(0);
 if (writeBuf_.allocated())
   val = writeBuf_.data_end_ - writeBuf_.current_;
 return val;  
}

// public method:
inline bool abstractCommHandle::buffered(void) const
{
 return 
   (reading() && readBuf_.allocated())
   || (writing() && writeBuf_.allocated());
}	

// protected methods: ------------------
inline bool abstractCommHandle::readBuffered(bool ignore_sync) const
{
 return 
   readBuf_.allocated() && (ignore_sync || !readBufIO_);
}

inline bool abstractCommHandle::writeBuffered(bool ignore_flush) const
{
 return 
   writeBuf_.allocated() && (ignore_flush || !writeBufIO_);
}	
// --------------------------------------  
  
inline bool abstractCommHandle::valid(void) const { return valid_; }

inline bool abstractCommHandle::valid(bool state)
{
 bool val(valid_);
 valid_ = state;
 return val;
}

inline unsigned short abstractCommHandle::flags(void) const { return flags_; }

inline unsigned short abstractCommHandle::flags(unsigned short newFlags)
{
 unsigned short val(flags_);
 flags_ = newFlags;
 return val;
}	
	
inline bool abstractCommHandle::reading(void) const { return flags_&READING; }

inline bool abstractCommHandle::writing(void) const { return flags_&WRITING; }
	
// assume "h" on heap:
inline void abstractCommHandle::close(abstractCommHandle*& h)
{
 if (h != NULL){
   h->close();
   delete h;
 }  
 h = NULL;	 
}
		
inline size_t read(void *ptr, size_t size, size_t nitems, abstractCommHandle* h) { return h->read(ptr, size, nitems); }

inline size_t write(const void *ptr, size_t size, size_t nitems, abstractCommHandle* h) { return h->write(ptr, size, nitems); }

inline void close(abstractCommHandle*& h) { abstractCommHandle::close(h); }


		
// from an existing (open) stream:
inline fileHandle* fileHandle::clone(FILE *fp)
{
 return fileHandle::clone(fileno(fp));
}


inline abstractCommHandle* open(const char *filename, const char *mode) { return fileHandle::open(filename, mode); }

inline abstractCommHandle* open(const char *pathname, int flags, mode_t mode) { return fileHandle::open(pathname, flags, mode); }

inline void seek(abstractCommHandle* fp, long offset, int whence)
{
 fileHandle *fp_(dynamic_cast<fileHandle*>(fp));
 if (NULL == fp_)
   throw std::runtime_error("commUtil::seek: cast to derived fileHandle* fails");
 fp_->seek(offset,whence); 
}

inline long tell(abstractCommHandle* fp)
{
 fileHandle *fp_(dynamic_cast<fileHandle*>(fp));
 if (NULL == fp_)
   throw std::runtime_error("commUtil::tell: cast to derived fileHandle* fails");
 return fp_->tell(); 
}

inline size_t memoryCommHandle::buf_remaining(void)const
{
 size_t N(0);
 if (buf_ != NULL){
   const unsigned char *buf_end_(reinterpret_cast<const unsigned char*>(buf_) + buf_size_);
   N = buf_end_ - reinterpret_cast<const unsigned char*>(cp_);
 }  
 return N;
}

// gcc-4.1.2 work-around:  use of std::copy with uchar_alias* causes "symbol already defined" assembler error:
inline void memoryCommHandle::copy_(const uchar_alias* src_b, const uchar_alias* src_e, uchar_alias* dest_b)
{
  #if ((__GNUC__ == 4) && (__GNUC_MINOR__ == 1) && (__GNUC_PATCHLEVEL__ == 2))  
  // dissallow overlap:
  assert( (dest_b + (src_e - src_b) <= src_b) || (dest_b >= src_e) );
  uchar_alias *dest(dest_b);
  for(const uchar_alias *src = src_b; src != src_e; ++src, ++dest)
    *dest = *src;
  #else
  std::copy(src_b, src_e, dest_b); // I assume that this is faster...
  #endif  
}


#if defined(__USE_MPI) // ----------------------- MPI ONLY SECTION ------------------------------------------

inline int processCommHandle::rank(void)const
{ return rank_; }

inline const MPI::Comm& processCommHandle::comm(void)const
{ return *pComm_; }
   		
inline int processCommHandle::tag(void)const { return tag_; }

inline int processCommHandle::tag(int newTag)
{ 
 int val(tag_);
 tag_ = newTag;
 return val;
}

inline processCommHandle::TRANSFER_KIND processCommHandle::transfer_type(void)const { return transfer_type_; }

inline processCommHandle::TRANSFER_KIND  processCommHandle::transfer_type(TRANSFER_KIND newTransferType)
{ 
 TRANSFER_KIND val(transfer_type_);
 transfer_type_ = newTransferType;
 return val;
}

#if defined(__USE_PTHREAD) && 0
inline bool processCommHandle::signalThreads(void)
{ return signalThreads_; }

inline bool processCommHandle::signalThreads(bool flag)
{
 bool val = signalThreads_;
 signalThreads_ = flag;
 return val;
}

inline bool processCommHandle::threadDataReady(void)
{ return threadDataReady_; }
#endif

#endif // ----------------------- end, MPI ONLY SECTION ------------------------------------------


// **********************************************************************
//
// _generic_ linalg writeBinary and readBinary.
//

#if 0
// --------------- non-TR1 versions: ------------------------------------------------------
// note: these _usually_ confuse the compiler, because of ambiguous resolution: -----------
template<class U>
inline bool writeBinary(abstractCommHandle *fp, const U& u)
{ return u.writeBinary(fp); }

template<class U>
inline bool readBinary(abstractCommHandle *fp, U& u )
{ return u.readBinary(fp); }
#else
// ---------------- versions using TR1: ---------------------------------------------------
template <class U>
inline bool writeBinary_POD_dispatch_(abstractCommHandle *fp, const U& u, std::tr1::false_type)
{ return u.writeBinary(fp); }
template <class U>
inline bool writeBinary_POD_dispatch_(abstractCommHandle *fp, const U& u, std::tr1::true_type)
{ 
  bool status(true);
  status = (status && (1 == write(&u, sizeof(u), 1, fp) ) );	
  return status; 
}

template<class U>
inline bool writeBinary(abstractCommHandle *fp, const U& u)
{ return writeBinary_POD_dispatch_(fp, u, typename std::tr1::is_pod<U>::type()); }


template <class U>
inline bool readBinary_POD_dispatch_(abstractCommHandle *fp, U& u, std::tr1::false_type)
{ return u.readBinary(fp); }
template <class U>
inline bool readBinary_POD_dispatch_(abstractCommHandle *fp, U& u, std::tr1::true_type)
{ 
  bool status(true);
  status = (status && (1 == read(&u, sizeof(u), 1, fp) ) );	
  return status; 
}

template<class U>
inline bool readBinary(abstractCommHandle *fp, U& u)
{ return readBinary_POD_dispatch_(fp, u, typename std::tr1::is_pod<U>::type()); }


template <class U>
inline size_t binarySize_POD_dispatch_(U& u, std::tr1::false_type)
{ return u.binarySize(); }
template <class U>
inline size_t binarySize_POD_dispatch_(U& u, std::tr1::true_type)
{ return sizeof(u); }

template <class U>
inline size_t binarySize(const U& u)
{ return binarySize_POD_dispatch_(u, typename std::tr1::is_pod<U>::type()); }

#endif

// **********************************************************************

				
} // namespace commUtil

#endif // __commUtil_inline__h
