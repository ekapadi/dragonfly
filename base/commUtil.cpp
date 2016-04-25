// $Source: /usr/data0/leipzig_work/tmat_cvs/src/commUtil.cpp,v $

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

#if defined(__USE_PTHREAD)
  #include <pthread.h>
#endif

#if defined(__USE_MPI)
  #include <mpi.h>
#endif

//  includes to check file v. directory, and status:
#include <sys/types.h> 
#include <sys/stat.h> 
#include <errno.h>
// 

#include <assert.h>
#include <unistd.h>
#include <fcntl.h>

#include <cstdlib>
#include <cstdio>
#include <stdexcept>
#include <iostream>
#include <iomanip>
#include <string>
#include <sstream>
#include <vector>

#include <typeinfo>
#include <tr1/type_traits>

#include "commUtil.h"

using std::cout;
using std::endl;

namespace commUtil{


comm_error& comm_error::operator=(const comm_error& other)
{ 
  base_class::operator=(static_cast<const base_class&>(other));     
  return *this;
}

comm_error::~comm_error(void)
{ }

comm_error::comm_error(void)
  : base_class("unspecified comm_error")
{ }

comm_error::comm_error(const comm_error& other)
  : base_class(static_cast<const std::runtime_error&>(other))
{ }

comm_error::comm_error(const std::string& msg)
  : base_class(msg)
{ }


size_t abstractCommHandle::bufferedRead(void *ptr, size_t size, size_t nitems)
{
  assert(buffered()); // this is an internal method  

  size_t transferSize(size*nitems);
  unsigned char *ptr_(reinterpret_cast<unsigned char*>(ptr));
  while (transferSize != 0){
    if (!readBufAvail())
      sync(); 
    if (!readBufAvail())
      break; 
    size_t increment(transferSize);
    if (increment > readBufAvail())
      increment = readBufAvail();
    std::copy(readBuf_.current_, readBuf_.current_ + increment, ptr_);
    readBuf_.current_ += increment;
    transferSize -= increment;
    ptr_ += increment;
  } 

  const size_t 
    bytesTransferred(size*nitems - transferSize),
    itemsTransferred(bytesTransferred/size);
    
  // don't allow partial transfer of items:  
  assert(bytesTransferred%size==0);
    
  return itemsTransferred;
}

    
size_t abstractCommHandle::bufferedWrite(const void *ptr, size_t size, size_t nitems)
{
  assert(buffered()); // this is an internal method
  
  size_t transferSize(size*nitems);
  const unsigned char *ptr_(reinterpret_cast<const unsigned char*>(ptr));
  while (transferSize != 0){
    if (!writeBufAvail())
      flush(); 
    if (!writeBufAvail())
      break;
    size_t increment(transferSize);
    if (increment > writeBufAvail())
      increment = writeBufAvail();
    std::copy(ptr_, ptr_+increment, writeBuf_.current_);
    writeBuf_.current_ += increment;
    transferSize -= increment;
    ptr_ += increment;
  }
  
  const size_t 
    bytesTransferred(size*nitems - transferSize),
    itemsTransferred(bytesTransferred/size);
    
  // don't allow partial transfer of items:  
  assert(bytesTransferred%size==0);
  return itemsTransferred;  
}

   
void abstractCommHandle::sync(void)
{
 if (!valid() || !reading())
   throw std::runtime_error("abstractCommHandle::sync: invalid comm handle or handle not open for read");
 if (readBuffered()){
   readBufIO_ = true; // internal transfer flag. 
   
   unsigned char *transfer_buf(NULL);
   size_t transfer_size(0);
   if (transferBufferHeader()){
     transfer_buf = readBuf_.begin_;
     transfer_size = readBuf_.size();
   }
   else{
     transfer_buf = readBuf_.data_begin_;
     transfer_size = readBuf_.size() - readBuf_.header_size();  // note: data_size() would return encoded size of last transfer (not what we want...)
   }
   
   transfer_size = read(transfer_buf,1,transfer_size);  // max: transfer_size bytes (i.e. current value of transfer_size), actual number depends on MODE.   
   readBufIO_ = false;
   if (transferBufferHeader()){
     if (transfer_size == 0) // transfer_size == 0 is not necessarily an error, read errors throw exception,
       throw std::runtime_error("abstractCommHandle::sync: zero sized read transfer -- unable to encode data size");
     readBuf_.init(); // this obtains the size of the transfer from the buffer's own header (compatible with both MPI::Recv, _and MPI::Bcast).
   }
   else
     readBuf_.init(transfer_size);
   
   #if 0
   // *** DEBUG ***
   cout<<"sync, transferred: "<<transfer_size<<" bytes"
       <<", available data to read: "<<readBufAvail()<<endl;
   #endif
 }   
}
      
void abstractCommHandle::flush(void)
{
 if (!valid() || !writing())
   throw std::runtime_error("abstractCommHandle::flush: invalid comm handle or handle not open for write");
 if (writeBuffered() && writeBuf_.used()){
 
   const size_t 
     data_size(partialBufferTransfer()? writeBuf_.used(): writeBuf_.data_size()),
     transfer_size(transferBufferHeader()? data_size + writeBuf_.header_size(): data_size);
   writeBuf_.set_data_size(data_size); // mark the actual data-size in the buffer header.
   
   unsigned char *transfer_buf(transferBufferHeader()? writeBuf_.begin_: writeBuf_.data_begin_);
   writeBufIO_ = true; // internal transfer flag
   if (transfer_size != write(transfer_buf,1,transfer_size))
     throw std::runtime_error("abstractCommHandle::flush: I/O error");
   writeBufIO_ = false;
   
   // reset the data size as marked in the buffer:
   writeBuf_.set_data_size(writeBuf_.size() - writeBuf_.header_size());
   writeBuf_.init(); 
   #if 0
   // *** DEBUG ***
   cout<<"flush, transferred "<<transfer_size<<" bytes"
       <<", available space to write: "<<writeBufAvail()<<endl;
   #endif    
 }   
}    
    
bool abstractCommHandle::buffered(bool flag, size_t bytes)
{
 if (!valid())
   throw std::runtime_error("abstractCommHandle::buffered: invalid comm handle");
   
 bool val(buffered()); // return previous state

 if (reading()){
   if (flag){
     readBuf_.alloc(bytes);
     readBuf_.current_ = readBuf_.data_end_; // set to empty at start
   }
   else
     readBuf_.free();
 }     
 if (writing()){
   if (flag)
     writeBuf_.alloc(bytes);
   else
     writeBuf_.free();
 } 
 return val;   
}   

// attach to external buffer, don't transfer ownership:
bool abstractCommHandle::buffered(unsigned char *begin, unsigned char *end)
{
 if (!valid())
   throw std::runtime_error("abstractCommHandle::buffered: invalid comm handle");
 if (reading() && writing())
   throw std::runtime_error("abstractCommHandle::buffered: read/write mode requires two buffers");

 bool val(buffered()); // return previous state

 if (reading()){
   readBuf_.attach(begin,end);
   readBuf_.current_ = readBuf_.data_end_; // set to empty at start
 }
      
 if (writing())
   writeBuf_.attach(begin,end);

 return val;   
}   

// attach to external buffers, don't transfer ownership:
bool abstractCommHandle::buffered(unsigned char *beginRead, unsigned char *endRead,
              unsigned char *beginWrite, unsigned char *endWrite)
{
 if (!valid())
   throw std::runtime_error("abstractCommHandle::buffered: invalid comm handle");
 if (!(reading() && writing()))
   throw std::runtime_error("abstractCommHandle::buffered: read-only or write-only mode requires one buffer");
 if (!( ((beginRead < endRead) && (endRead <= beginWrite) && (beginWrite < endWrite))
       ||((beginWrite < endWrite) && (endWrite <= beginRead) && (beginRead < endRead))))
   throw std::runtime_error("abstractCommHandle::buffered: read and write buffers must be disjoint");
    
 bool val(buffered()); // return previous state
             
 readBuf_.attach(beginRead,endRead);
 readBuf_.current_ = readBuf_.data_end_; // set to empty at start
 writeBuf_.attach(beginWrite,endWrite);

 return val;   
}        


// encode the data-size of the buffer portion used in the transfer data stream itself:
bool abstractCommHandle::transferBufferHeader(void) const
{
 return false;
}

// allow transfer of fractional buffers during "sync" and "flush" methods:
bool abstractCommHandle::partialBufferTransfer(void) const
{
 return true;
}
    
void abstractCommHandle::close(void)
{
  // in general case: nothing to do; 
  //   however, do not _require_ implementation of this method.
}
  
abstractCommHandle::abstractCommHandle(void)
  : valid_(false), flags_(0),
    readBufIO_(false), writeBufIO_(false)    
{}	

abstractCommHandle::~abstractCommHandle(void)
{ 
 valid(false); 
 flags(0);
 readBuf_.free();
 writeBuf_.free();
}

// encode the data-size of the buffer portion used in the transfer data stream itself:
bool fileHandle::transferBufferHeader(void) const
{
 return false;
}

// allow transfer of fractional buffers during "sync" and "flush" methods:
bool fileHandle::partialBufferTransfer(void) const
{
 return true;
}

fileHandle::fileHandle(void)
  :fd_(-1)
{}

// allocate on heap:
fileHandle* fileHandle::open(const char *filename, const char *mode)
{
 int flags_(O_SYNC); 

 if (NULL != std::strchr(mode,'+'))
   flags_ = flags_ | O_RDWR | O_CREAT; 
 else
 if (NULL != std::strchr(mode,'r'))
   flags_ = flags_ | O_RDONLY;
 else
 if (NULL != std::strchr(mode,'w'))
   flags_ = flags_ | O_WRONLY | O_CREAT;
 else
   throw std::runtime_error("fileHandle::open: one of \"+\", \"r\", or \"w\" must be specified in the mode string");  
	
 if (NULL != std::strchr(mode,'a'))
   flags_ = flags_ | O_APPEND;
 else
 if (NULL != std::strchr(mode,'w'))
   flags_ = flags_ | O_TRUNC;  

 mode_t mode_ = S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH; // "0666" permissions, as modified by process umask value.
 
 return fileHandle::open(filename, flags_, mode_);
}

// allocate on heap:
fileHandle* fileHandle::open(const char* pathname, int flags, mode_t mode)
{
 fileHandle *h(new fileHandle);

 #if 1
 // *** DEBUG ***
 size_t state(0);
 #endif 
 
 if (NULL != h){
	 h->fd_ = ::open(pathname, flags, mode);

   #if 1
   // *** DEBUG ***
   state = 1;
   #endif
   
	 if (-1 != h->fd_){
  	 if ((flags & 0x3) == O_RDWR){
		   h->flags(READING | WRITING);
		   h->valid(true);

   #if 1
   // *** DEBUG ***
   state = 2;
   #endif

		 } 
		 else
		 if ((flags & 0x3) == O_RDONLY){
	  	 h->flags(READING);
		   h->valid(true);

   #if 1
   // *** DEBUG ***
   state = 3;
   #endif

		 }
		 else
		 if ((flags & 0x3) == O_WRONLY){
	  	 h->flags(WRITING);
		   h->valid(true);

   #if 1
   // *** DEBUG ***
   state = 4;
   #endif

		 }
   #if 1
   // *** DEBUG ***
   state = 5;
   #endif     
	 }
 }
 
 if (NULL == h || !h->valid()){
	   if (NULL != h){
		   delete h;
		   h = NULL; // could return NULL if fp invalid, but for now, we will throw, to be consistent with processCommHandle::open
		 }
     #if 0
		 throw std::runtime_error("fileHandle::open: unable to open file ") + pathname;
     #else
     // *** DEBUG ***
     std::ostringstream oss;
     oss<<"fileHandle::open: unable to open file "<<pathname<<"\n"
       <<"  error string: "<<strerror(errno)<<"\n"
       <<"  state: "<<state<<"\n"
       <<"  flags: "<<std::oct<<flags<<"\n"
       <<"  mode:  "<<std::oct<<mode;
     throw std::runtime_error(oss.str());
     #endif
 } 	 

 return h;
}

// file exists:
bool fileHandle::exists(const char *fileName)
{
  struct stat sTestBuf;
  bool val=true;
  
  if ( -1==stat( fileName, &sTestBuf ) && (ENOENT==errno) ) // file not found in any form
     val = false;
     
  return val;
}

// is directory:
bool fileHandle::isDirectory(const char *fileName)
{
  struct stat sTestBuf;
  bool val=true;
  
  if ( -1==stat( fileName, &sTestBuf ) && (ENOENT==errno) ) // file not found in any form
     val = false;
  else if ( !S_ISDIR(sTestBuf.st_mode) )
     val = false;
		 
  return val;  
}

// file size:
size_t fileHandle::sizeOnMedia(const char *fileName)
{
  struct stat sTestBuf;
  
  if ( -1==stat( fileName, &sTestBuf ) && (ENOENT==errno) ) // file not found in any form
    throw std::runtime_error("fileHandle::sizeOnMedia: file not found");  
     
  return static_cast<size_t>(sTestBuf.st_size);
}

// from an existing (open) file-descriptor:
fileHandle* fileHandle::clone(int fd)
{
 fileHandle *h(new fileHandle);
 
 if (NULL != h){
	 h->fd_ = fd;

	 if (-1 != h->fd_){
     const int OPEN_MODE(O_RDONLY|O_WRONLY|O_RDWR); // mask: these are not flag-bits
     int flags = /* std:: */ fcntl(fd, F_GETFL );
		 if (-1 == flags)
		   throw std::runtime_error("fileHandle::clone: fcntl error");
			 
  	 if ((flags & OPEN_MODE) == O_RDONLY)
	  	 h->flags(READING);
	   else
		 if ((flags & OPEN_MODE) == O_WRONLY)
	  	 h->flags(WRITING);
		 else
		 if ((flags & OPEN_MODE) == O_RDWR)
		   h->flags(READING|WRITING); 
		 else
		   throw std::runtime_error("fileHandle::clone: unsupported file open state ( not read-only, write-only, or read/write )"); 
			 
  	 h->valid(true);
	 }
   	
 }	 

 return h;
}

size_t fileHandle::read(void *ptr, size_t size, size_t nitems)
{ 
 if (!valid() || !reading() ) 
   throw std::runtime_error("fileHandle::read: invalid handle or not open for read");
 size_t val(0);
 if (readBuffered())
   val = bufferedRead(ptr, size, nitems);
 else{
   size_t count(size*nitems);
   
   if (readBuffered(true)){
     // for efficiency reasons, this check is _only_ performed in the buffered case (true flag => also check if reading during actual "sync"):
     
     // adjust for actual amount of data available:
     // (Note: this shouldn't really be necessary, but for some reason read(2) hangs if more data is requested than exists in a file?!
     //    It should just return ssize_t(-1) according to the documentation...)
     struct stat finfo;
     off_t current_offset(static_cast<off_t>(-1));
     if ((-1 == ::fstat(fd_, &finfo)) || (static_cast<off_t>(-1) == (current_offset = ::lseek(fd_, 0, SEEK_CUR))))
       throw std::runtime_error("fileHandle::read: I/O error: fstat or lseek error return");
     if (count > static_cast<size_t>(finfo.st_size - current_offset))
       count = static_cast<size_t>(finfo.st_size - current_offset);        
   }
   assert(count <= SSIZE_MAX);
   
   unsigned char *buf(reinterpret_cast<unsigned char*>(ptr));
   // "read" can be interrupted, in which case it returns a partial count:
   // (note: bytes_read == zero is a possibility, which also => a zero-length file => dissallow it (or this is an infinite loop)).
   do{
     ssize_t bytes_read = ::read(fd_, buf, count);
     if ((bytes_read != ssize_t(-1)) && (bytes_read > 0)){
       count -= bytes_read;
       buf += bytes_read;
     } 
     else
       throw std::runtime_error("fileHandle::read: I/O error (or zero-length file)");  
   }
   while (count > 0);
   val = nitems;
 }
 return val;
}

size_t fileHandle::write(const void *ptr, size_t size, size_t nitems)
{ 
 if (!valid() || !writing() ) 
   throw std::runtime_error("fileHandle::write( const void*, size_t, size_t ): invalid handle or not open for write");
 size_t val(0);
 if (writeBuffered())
   val = bufferedWrite(ptr, size, nitems);
 else{  
   size_t count(size*nitems);
   assert(count <= SSIZE_MAX);
   const unsigned char *buf(reinterpret_cast<const unsigned char*>(ptr));
   // "write" can be interrupted, in which case it returns a partial count:
   do{
     ssize_t bytes_written = ::write(fd_, buf, count);
     if (bytes_written != ssize_t(-1)){
       count -= bytes_written;
       buf += bytes_written;
     } 
     else
       throw std::runtime_error("fileHandle::write: I/O error");  
   }
   while (count > 0);
   val = nitems; 
 }   
 return val;
}

void fileHandle::close(void)
{
 if (fd_ != -1){
   if (-1 == ::close(fd_))
     throw std::runtime_error("fileHandle::~fileHandle: error closing file-descriptor");
	 fd_ = -1;
 } 
}
    
fileHandle::~fileHandle(void)		
{
 #if 0
   close();
 #else
   // Invalidate, but require _explicit_ close(handle)
   //   _otherwise_ cannot daisy-chain file operations from *FILE or <file descriptor>.
   fd_ = -1;
 #endif
}


void fileHandle::seek(long offset, int whence)
{
 if (!valid())
   throw std::runtime_error("fileHandle::seek( long, int ): invalid handle");
 if (off_t(-1) == ::lseek(fd_, static_cast<off_t>(offset), whence))
   throw std::runtime_error("fileHandle::seek: I/O error"); 
}
    
long fileHandle::tell(void)
{
 if (!valid())
   throw std::runtime_error("fileHandle::tell(void): invalid handle");
 off_t offset = ::lseek(fd_, 0, SEEK_CUR);
 if (offset == off_t(-1))
   throw std::runtime_error("fileHandle::tell: I/O error");
     
 long pos = static_cast<long>(offset);
 if (buffered()){
   if (reading() && writing())
     throw std::runtime_error("fileHandle::tell: not implemented for simultaneous buffered read/write mode");
   if (reading())
     pos -= readBufAvail();
   else
   if (writing())
     pos += writeBufAvail();  
 }
 return pos;
} 

// ----------------------- memoryCommHandle: -----------------------------------------------

memoryCommHandle::memoryCommHandle(void) // protected constructor
 :buf_(NULL), cp_(NULL), buf_size_(0)
{ }

// encode the data-size of the buffer portion used in the transfer data stream itself:
bool memoryCommHandle::transferBufferHeader(void) const
{
 return false;
} 
   
bool memoryCommHandle::partialBufferTransfer(void) const
{ return true; }
		
size_t memoryCommHandle::read(void *ptr, size_t size, size_t nitems)
{ 
 if (!valid() || !reading() ) 
   throw std::runtime_error("memoryCommHandle::read( void*, size_t, size_t ): invalid handle or not open for read");
 size_t val(0);
 if (readBuffered())
   val = bufferedRead(ptr, size, nitems);
 else{
   if (size*nitems > buf_remaining())
     throw std::runtime_error("memoryCommHandle::read: buffer underflow");   
   uchar_alias *ptr_(reinterpret_cast<uchar_alias*>(ptr));
   copy_(cp_, cp_ + size*nitems, ptr_);
   cp_ += size*nitems;
   val = nitems; 
 }
 return val;
}

size_t memoryCommHandle::write(const void *ptr, size_t size, size_t nitems)
{
 if (!valid() || !writing() ) 
   throw std::runtime_error("memoryCommHandle::write( const void*, size_t, size_t ): invalid handle or not open for write");
 size_t val(0);
 if (writeBuffered())
   val = bufferedWrite(ptr, size, nitems);
 else{
   if (size*nitems > buf_remaining())
     throw std::runtime_error("memoryCommHandle::write: buffer overflow");
   const uchar_alias *ptr_(reinterpret_cast<const uchar_alias*>(ptr));
   copy_(ptr_, ptr_ + size*nitems, cp_);
   cp_ += size*nitems;
   val = nitems; 
 }
 return val;
}

// allocate on heap (local structures, _excepting_ buffer itself):
memoryCommHandle* memoryCommHandle::open(void *buffer, size_t buf_size, const char *mode)
{
 memoryCommHandle *h(new memoryCommHandle);

 if (NULL != h){
   h->buf_ = reinterpret_cast<uchar_alias*>(buffer);
   h->cp_ = reinterpret_cast<uchar_alias*>(buffer);
   h->buf_size_ = buf_size;	 

	 if (NULL != h->buf_){
  	 if (NULL != std::strchr(mode,'+')){
		   h->flags(READING | WRITING);
		   h->valid(true);
		 } 
		 else
		 if (NULL != std::strchr(mode,'r')){
	  	 h->flags(READING);
		   h->valid(true);
		 }
		 else
		 if ((NULL != std::strchr(mode,'w'))
		     || (NULL != std::strchr(mode,'a'))) {
	  	 h->flags(WRITING);
		   h->valid(true);
		 }
	 }
 }
 
 if (NULL == h || !h->valid()){
	   if (NULL != h){
		   delete h;
		   h = NULL; // could return NULL if fp invalid, but for now, we will throw, to be consistent with processCommHandle::open
		 }
		 throw std::runtime_error("memoryCommHandle::open: invalid buffer or mode specified");
 } 	 

 return h;
}

memoryCommHandle::~memoryCommHandle(void)		
{
  // mark to assist "double-free" error detection:
	buf_ = NULL;
  cp_ = NULL;
  buf_size_ = 0;
} 

// -----------------------------------------------------------------------------------------

#if defined (__USE_MPI) // -------------------------- MPI ONLY SECTION -------------------------------------

// encode the data-size of the buffer portion used in the transfer data stream itself:
bool processCommHandle::transferBufferHeader(void) const
{
 return true;
}

// allow transfer of fractional buffers during "sync" and "flush" methods:
bool processCommHandle::partialBufferTransfer(void) const
{
 if (POINT2POINT == transfer_type())
   return true;
 return false;
}


processCommHandle::processCommHandle(int rank, const MPI::Comm& comm, int tag, TRANSFER_KIND transfer_type)
  : rank_(rank), /* pComm_( &comm.Clone() ), */ tag_(tag), transfer_type_(transfer_type)
{
 #if 1
 pComm_ = &comm; 
 #else
 pComm_ = &(MPI::COMM_WORLD); // *** DEBUG ***, comm.Clone() hangs for some reason with comm==MPI_COMM_WORLD
 #endif
}

size_t processCommHandle::read(void *ptr, size_t size, size_t nitems)
{ 
 if (!valid() || !reading()) 
   throw std::runtime_error("processCommHandle::read( void*, size_t, size_t ): invalid handle or not open for read");

 size_t val(0);
 // MPI_TYPE_CONTIGUOUS(size, MPI_CHARACTER, datatype) 
 // MPI_TYPE_COMMIT(datatype) 
 // int MPI_Recv(void* buf, int count, MPI_Datatype datatype, int source, int tag, MPI_Comm comm, MPI_Status *status)  
 // MPI_TYPE_FREE(datatype);

 #if defined(__USE_PTHREAD) && 0
 if ( (rank() == 0) && (comm().Get_rank() == 0) &&  signalThreads() ){
			 while (!threadDataReady_){
		  	 if ( pthread_cond_wait(&threadSignalCond_, &threadSignalMutex_) ){
			  	 pthread_mutex_unlock(&threadSignalMutex_);
			  	 throw std::runtime_error("processCommHandle::read: pthread_cond_wait error return");
				 }
				 pthread_mutex_unlock(&threadSignalMutex_); // do not keep the mutex locked...
			 }
 }			 
 #endif 

 if (readBuffered())
   val = bufferedRead(ptr, size, nitems);
 else {
   switch (transfer_type()){
     case POINT2POINT:
     {
       MPI::Status status;
       size_t bytesRead(0);
       #if defined(__CHUNK_TRANSFER__) // break into small chunks:
       for(unsigned char *buf = reinterpret_cast<unsigned char *>(ptr), 
                         *bufEnd = reinterpret_cast<unsigned char *>(ptr) + nitems*size;
           buf < bufEnd;
           buf += CHUNK_SIZE_){
         size_t readSize(CHUNK_SIZE_);
         if (buf + readSize > bufEnd)
           readSize = (bufEnd - buf);
         pComm_->Recv(buf, readSize, MPI::BYTE, rank_, tag_, status);
         bytesRead += status.Get_count(MPI::BYTE);       
       }
       #else   
       pComm_->Recv(ptr, nitems*size, MPI::BYTE, rank_, tag_, status);
       bytesRead = status.Get_count(MPI::BYTE);
       #endif
       assert(bytesRead%size == 0); // don't allow read of partial items
       val = bytesRead/size;        
     }
     break;

     case BROADCAST:
     {
       // if current rank is same as rank_, we are the root, and for "read" this is an error. 
       // possibly also so for an intraComm: root's group: *** DEBUG *** to be implemented...
       if (pComm_->Get_rank() == rank_)
         throw std::runtime_error("processCommHandle::read: can't BROADCAST mode to self"); 
       
       // only an "MPI::Intracomm" can Bcast:
       const MPI::Intracomm* pIntracomm(dynamic_cast<const MPI::Intracomm*>(pComm_));
       if (NULL == pIntracomm)
         throw std::runtime_error("processCommHandle::read: cast of communicator to derived type \"Intracomm\" fails");
       
       #if defined(__CHUNK_TRANSFER__) // break into small chunks:
       for(unsigned char *buf = reinterpret_cast<unsigned char *>(ptr), 
                               *bufEnd = reinterpret_cast<unsigned char *>(ptr) + nitems*size;
           buf < bufEnd;
           buf += CHUNK_SIZE_){
         size_t readSize(CHUNK_SIZE_);
         if (buf + readSize > bufEnd)
           readSize = (bufEnd - buf);
         pIntracomm->Bcast(buf, readSize, MPI::BYTE, rank_);       
       }     
       #else
       // if current rank is same as rank_, we are the root, and for "read" this is an error. 
       // possibly also so for an intraComm: root's group: *** DEBUG *** to be implemented...
        pIntracomm->Bcast(ptr, nitems*size, MPI::BYTE, rank_);
       #endif
       val = nitems; // for collective communication, this is an MPI restriction: full transfer is mandated.
     }
     break;

     default:
       throw std::runtime_error("processCommHandle::read: unknown transfer_type");
     // break;
   }
 }

 #if defined(__USE_PTHREAD) && 0
 // signal some other (usually rank-0 multithreaded MPI) thread that data has been sent):
 if ( (rank() == 0)  && (comm().Get_rank() == 0) &&  signalThreads()){
   pthread_mutex_lock(&threadSignalMutex_);
   threadDataReady_ = false;
   if (pthread_cond_signal(&threadSignalCond_)){
	   pthread_mutex_unlock(&threadSignalMutex_);
		 throw std::runtime_error("processCommHandle::write: pthread_cond_signal error");
	 }
	 pthread_mutex_unlock(&threadSignalMutex_);
 }
 #endif
 
 return val; // note: non-error return behavior for std::fread; however, errors throw exceptions
}

size_t processCommHandle::write(const void *ptr, size_t size, size_t nitems)
{ 
 if (!valid() || !writing()) 
   throw std::runtime_error("processCommHandle::write( const void*, size_t, size_t ): invalid handle or not open for write");

 size_t val(0);
 
 #if defined(__USE_PTHREAD) && 0
 if ( (rank() == 0) && (comm().Get_rank() == 0) &&  signalThreads() ){
			 while (threadDataReady_){
		  	 if ( pthread_cond_wait(&threadSignalCond_, &threadSignalMutex_) ){
			  	 pthread_mutex_unlock(&threadSignalMutex_);
			  	 throw std::runtime_error("processCommHandle::read: pthread_cond_wait error return");
				 }
				 pthread_mutex_unlock(&threadSignalMutex_); // do not keep the mutex locked...
			 }
 }			 
 #endif 
 
 if (writeBuffered())
   val = bufferedWrite(ptr, size, nitems);
 else{
   switch (transfer_type()){
     case POINT2POINT:
     {
       #if defined(__CHUNK_TRANSFER__) // break into small chunks:
       for(const unsigned char *buf = reinterpret_cast<const unsigned char *>(ptr), 
           *bufEnd = reinterpret_cast<const unsigned char *>(ptr) + nitems*size;
           buf < bufEnd;
           buf += CHUNK_SIZE_){
         size_t writeSize(CHUNK_SIZE_);
         if (buf + writeSize > bufEnd)
           writeSize = (bufEnd - buf);
         pComm_->Ssend(buf, writeSize, MPI::BYTE, rank_, tag_);
       }
       #else    
       pComm_->Ssend(ptr, nitems*size, MPI::BYTE, rank_, tag_);
       #endif
       val = nitems;
     }
     break;

     case BROADCAST:
     {       
       // only an "MPI::Intracomm" can Bcast:
       const MPI::Intracomm* pIntracomm(dynamic_cast<const MPI::Intracomm*>(pComm_));
       if (NULL == pIntracomm)
         throw std::runtime_error("processCommHandle::write: cast of communicator to derived type \"Intracomm\" fails");
   
       #if defined(__CHUNK_TRANSFER__) // break into small chunks:
       for(const unsigned char *buf = reinterpret_cast<const unsigned char *>(ptr), 
           *bufEnd = reinterpret_cast<const unsigned char *>(ptr) + nitems*size;
           buf < bufEnd;
           buf += CHUNK_SIZE_){
         size_t writeSize(CHUNK_SIZE_);
         if (buf + writeSize > bufEnd)
           writeSize = (bufEnd - buf);
         pIntracomm->Bcast(const_cast<void*>(reinterpret_cast<const void*>(buf)), 
                       writeSize, MPI::BYTE, pComm_->Get_rank() /* (rank_==pComm_->Get_rank()? MPI_ROOT: MPI_PROC_NULL) */  );       
       }     
       #else
       // if current rank is same as rank_, we are the root, otherwise we are in the root's group (otherwise "write" in BROADCAST mode is an MPI error):
       pIntracomm->Bcast(const_cast<void*>(ptr), nitems*size, MPI::BYTE, pComm_->Get_rank() /* (rank_==pComm_->Get_rank()? MPI_ROOT: MPI_PROC_NULL) */ );
       #endif
       val = nitems; // for collective communication, this is an MPI restriction: full transfer is mandated.
     }
     break;

     default:
       throw std::runtime_error("processCommHandle::write: unknown transfer_type");
     // break;
   }
 }

 #if defined(__USE_PTHREAD) && 0
 // signal some other (usually rank-0 multithreaded MPI) thread that data has been sent):
 if ( (rank() == 0) && (comm().Get_rank() == 0) && signalThreads()){
   pthread_mutex_lock(&threadSignalMutex_);
   threadDataReady_ = true;
   if (pthread_cond_signal(&threadSignalCond_)){
	   pthread_mutex_unlock(&threadSignalMutex_);
		 throw std::runtime_error("processCommHandle::write: pthread_cond_signal error");
	 }
	 pthread_mutex_unlock(&threadSignalMutex_);
 }
 #endif
 
 return val;  // note: non-error return behavior for std::fwrite; however, errors throw exceptions 
}

bool processCommHandle::readTag(abstractCommHandle *h_, int& tag)
{ 
 processCommHandle *h( dynamic_cast<processCommHandle*>(h_) );
 
 if ((NULL == h) || !h->valid() || !h->reading()) 
   throw std::runtime_error("processCommHandle::readTag(abstractCommHandle*, int&): invalid handle or not open for read");

 #if defined(__USE_PTHREAD) && 0
 if ( (h->rank() == 0) && (h->comm().Get_rank() == 0) && signalThreads() ){
			 while (!threadDataReady_){
		  	 if ( pthread_cond_wait(&threadSignalCond_, &threadSignalMutex_) ){
			  	 pthread_mutex_unlock(&threadSignalMutex_);
			  	 throw std::runtime_error("processCommHandle::readTag: pthread_cond_wait error return");
				 }
				 pthread_mutex_unlock(&threadSignalMutex_); // do not keep the mutex locked...
			 }
 }			 
 #endif
 			 
 MPI::Status status;      
 #if 0
 unsigned char buf; // zero-length messages break some MPI implementations.
 #endif
 switch (h->transfer_type()){
   case POINT2POINT:   
     // Assume that there is a tag waiting, however if incoming tag is MPI::ANY_TAG, need to find out _what _ the tag is,
     //   otherwise, receive a _specific_ requested tag.
     if (tag == MPI::ANY_TAG){
	     h->pComm_->Probe(h->rank_, MPI::ANY_TAG, status);
	     if (status.Get_error() != MPI::SUCCESS)
	  	    throw std::runtime_error("processCommHandle::readTag: MPI::Comm::Probe error return");
	     tag = status.Get_tag();
     }

     h->pComm_->Recv( /* &buf, 1, */ NULL, 0, MPI::BYTE, h->rank_, tag); // actually receive the message
   break;
   
   case BROADCAST:
     // Apparently, "BCast" uses the tag field.  I don't see any way to "Probe" this, so dissallow it:
     throw std::runtime_error("processCommHandle::readTag: can't readTag with BROADCAST mode handle");
   // break;
   
   default:
     throw std::runtime_error("processCommHandle::readTag: unknown transfer_type");
   // break;
 }



 #if defined(__USE_PTHREAD) && 0
 // signal some other (usually rank-0 multithreaded MPI) thread that data has been sent):
 if ( (h->rank() == 0) && (h->comm().Get_rank() == 0) && signalThreads()){
   pthread_mutex_lock(&threadSignalMutex_);
   threadDataReady_ = false;
   if (pthread_cond_signal(&threadSignalCond_)){
	   pthread_mutex_unlock(&threadSignalMutex_);
		 throw std::runtime_error("processCommHandle::readTag: pthread_cond_signal error");
	 }
	 pthread_mutex_unlock(&threadSignalMutex_);
 }
 #endif
   
 return true;
}


bool processCommHandle::writeTag(abstractCommHandle *h_, int tag)
{ 
 processCommHandle *h( dynamic_cast<processCommHandle*>(h_) );

 if ((NULL == h) || !h->valid() || !h->writing()) 
   throw std::runtime_error("processCommHandle::writeTag(abstractCommHandle*, int): invalid handle or not open for write");

 #if defined(__USE_PTHREAD) && 0
 if ( (h->rank() == 0) && (h->comm().Get_rank() == 0) && signalThreads() ){
			 while (threadDataReady_){
		  	 if ( pthread_cond_wait(&threadSignalCond_, &threadSignalMutex_) ){
			  	 pthread_mutex_unlock(&threadSignalMutex_);
			  	 throw std::runtime_error("processCommHandle::readTag: pthread_cond_wait error return");
				 }
				 pthread_mutex_unlock(&threadSignalMutex_); // do not keep the mutex locked...
			 }
 }			 
 #endif
 		
 #if 0
 const unsigned char buf = '\0';  // don't send zero-length messages (breaks some MPI implementations).    		 
 #endif
 
 switch (h->transfer_type()){
   case POINT2POINT:
     h->pComm_->Ssend( /* &buf, 1, */ NULL, 0, MPI::BYTE, h->rank_, tag); 
   break;
   
   case BROADCAST:
     // Apparently, "BCast" uses the tag field.  I don't see any way to "Probe" this, so dissallow it:
     throw std::runtime_error("processCommHandle::writeTag: can't writeTag with BROADCAST mode handle");
   // break;
   
   default:
     throw std::runtime_error("processCommHandle::writeTag: unknown transfer_type");
   // break;
 }

 #if defined(__USE_PTHREAD) && 0
 // signal some other (usually rank-0 multithreaded MPI) thread that data has been sent):
 if ( (h->rank() == 0) && (h->comm().Get_rank() == 0) && signalThreads()){
   pthread_mutex_lock(&threadSignalMutex_);
   threadDataReady_ = true;
   if (pthread_cond_signal(&threadSignalCond_)){
	   pthread_mutex_unlock(&threadSignalMutex_);
		 throw std::runtime_error("processCommHandle::write: pthread_cond_signal error");
	 }
	 pthread_mutex_unlock(&threadSignalMutex_);
 }
 #endif
  
 return true;
}

// allocate on heap, allow only modes that make sense for processes:
processCommHandle* processCommHandle::open(const int rank, const MPI::Comm& comm, const char *mode, int tag, TRANSFER_KIND transfer_type)
{ 
 processCommHandle *h(new processCommHandle(rank, comm, tag, transfer_type));

 if (NULL != h){ 

     if (NULL != std::strchr(mode,'+'))
       // rank is _root_ rank:
		   h->flags(READING | WRITING); 
     else
     if (NULL != std::strchr(mode,'r'))
       // rank is _source_ rank:
       h->flags(READING);
     else
     if (NULL != std::strchr(mode,'w'))
		   // rank is _destination_ rank:
		   h->flags(WRITING);
     else
       throw std::runtime_error("processCommHandle::open: one of \"+\", \"r\", or \"w\" must be specified in the mode string");  

    h->valid(true);
 }

 // invalid handle => throw exception.
 
 return h;
}

processCommHandle::~processCommHandle(void)
{
 // how to close...
 flags(0);
#if 0
 delete (pComm_);
#endif
// abstractCommHandle::~abstractCommHandle();
}	

// static:
void *processCommHandle::MPIsendBuffer_ = NULL;

#if defined(__USE_PTHREAD) && 0		
// thread send-receive handoff (mostly for MPI rank-0 multi-threaded singleton):
bool processCommHandle::signalThreads_ = false; 
bool processCommHandle::threadDataReady_ = false; 
pthread_cond_t processCommHandle::threadSignalCond_ = PTHREAD_COND_INITIALIZER;
pthread_mutex_t processCommHandle::threadSignalMutex_ = PTHREAD_MUTEX_INITIALIZER;
#endif	

void processCommHandle::attachMPISendBuffer(size_t bufSize)
{
 if (NULL != MPIsendBuffer_)
   detachMPISendBuffer();
 MPIsendBuffer_ = std::malloc(bufSize);
 if (NULL == MPIsendBuffer_)
   throw std::runtime_error("processCommHandle::attachMPISendBuffer: can't allocate buffer");
 MPI::Attach_buffer(MPIsendBuffer_, bufSize); 
}

void processCommHandle::detachMPISendBuffer(void)
{
 if (NULL != MPIsendBuffer_){
   MPI::Detach_buffer(MPIsendBuffer_);
	 std::free(MPIsendBuffer_);
	 MPIsendBuffer_ = NULL;
 }
}

bool processCommHandle::MPIbuffered(void)
{ return !(NULL == MPIsendBuffer_); }

#endif // -------------------------- end, MPI ONLY SECTION -------------------------------------

// cannot use _other_ read/write methods which may not have been declared yet:
bool readBinary(abstractCommHandle* h, std::string& s)
{
 bool status(true);
 size_t N(0);
 status = (status && (1 == read(&N, sizeof(size_t), 1, h)));
 s.resize(N);
 for(std::string::iterator itS = s.begin(), itSEnd = s.end();
     status && (itS != itSEnd);
		 ++itS)
	 status = (status && (1 == read(&(*itS), sizeof(char), 1, h) ) );	
 return status;	  
}

bool writeBinary(abstractCommHandle* h, const std::string& s)
{
 bool status(true);
 const size_t N(s.size());
 status = (status && (1 == write(&N, sizeof(size_t), 1, h)));
 for(std::string::const_iterator itS = s.begin(), itSEnd = s.end();
     status && (itS != itSEnd);
		 ++itS)
	 status = (status && (1 == write(&(*itS), sizeof(char), 1, h) ) );	
 return status;	  
}

size_t binarySize(const std::string& s)
{
  size_t val(0);

  val += sizeof(size_t);
  // assume 8-bit characters:
  val += s.size();
	return val;
}  


} // namespace commUtil

// end commUtil.cpp
