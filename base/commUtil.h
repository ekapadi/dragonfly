#if !defined(__commUtil__h)
#define __commUtil__h

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


// #define __CHUNK_TRANSFER__

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/commUtil.h,v $

/**
 * @namespace commUtil
 * @ingroup parallel
 * @brief Unified low-level binary interface to disk-resident, memory-resident, and MPI inter-process data transfer.
 *
 * Using these methods to implement ``readBinary'', ''writeBinary'', and ``binarySize'' member or static methods for an end-user class
 * is sufficient to provide unified disk, memory, and/or inter-process I/O of class instances.
 */

namespace commUtil{


/**
 * @brief Base-class for binary-I/O exceptions.
 */
class comm_error: public std::runtime_error{
  private:

  public:
    typedef std::runtime_error base_class;
    
    comm_error& operator=(const comm_error& other);    
    
    virtual ~comm_error(void);
    
    comm_error(void);
    comm_error(const comm_error& other);
    comm_error(const std::string& msg);
}; // class comm_error 


/**
 * @brief An abstract base class to unify functionality of binary-I/O streams.
 *
 *    These classes include: file-streams, MPI::communicator-streams, and memory I/O-streams.
 */
class abstractCommHandle{
  protected:
    enum COMM_FLAGS {ANY_SET=0x0003, READING=0x0001, WRITING=0x0002};
			
  private:
	  bool valid_; // false if handle has been closed (or never opened).
		
		unsigned short flags_;


    class buffer{
      private:
        bool ownsData_;
        
      public:
        unsigned char *begin_, *data_begin_, *current_, *data_end_, *end_;     
 
        struct header{
          size_t size_;
        };

        inline void alloc(size_t size);
                
        // attach to external buffer,
        //   do not transfer ownership:
        inline void attach(unsigned char *begin, unsigned char *end);

        inline void free(void);

        inline void init(size_t data_size_);
        
        inline void init(void);

        inline size_t data_size(void)const;  // doesn't include header size

        inline void set_data_size(const size_t size);

        inline size_t size(void)const; // overall allocated size including header

        inline size_t header_size(void)const;    

        inline size_t used(void)const; // doesn't include header size

        inline bool allocated(void)const;      

        inline buffer(void);
    };
    
    buffer readBuf_, writeBuf_;
        
    bool readBufIO_, writeBufIO_; // true if I/O on a buffer itself, false otherwise.
      		
	protected:

		abstractCommHandle(void); // hide
	
	  inline bool valid(void) const;
		
		inline bool valid(bool state);
		
		inline unsigned short flags(void) const;
		
		inline unsigned short flags(unsigned short newFlags);
		
		inline bool reading(void)const;
		
		inline bool writing(void)const;

    size_t bufferedRead(void *ptr, size_t size, size_t nitems);
    
    size_t bufferedWrite(const void *ptr, size_t size, size_t nitems);

    // protected method: off during sync:
    inline bool readBuffered(bool ignore_sync=false) const;
    
    // protected method: off during flush:
    inline bool writeBuffered(bool ignore_flush=false) const;

    // protected: to allow use by fileHandle::seek
    inline size_t readBufAvail(void) const; 

    // protected: to allow use by fileHandle::seek    
    inline size_t writeBufAvail(void) const;    
    
    // the following are really "flags" governing the behavior of "sync" and "flush"
    //   during buffered data transfer:
    
    // "transferBufferHeader": transfer the buffer header along with the rest of the data
    //    this indicates that the data-size of the buffer portion used will be encoded into the data stream itself:
    virtual bool transferBufferHeader(void) const;  
    
    // "partialBufferTransfer": transfer only whole buffer blocks (or not): some modes (namely MPI BROADCAST) cannot
    //    encode any data-size information, so for these modes the entire buffer must be transfered, even if only a 
    //    portion of it is used (note that this is _independent_ of whether or not the data-size is _also_ encoded on the stream).
    virtual bool partialBufferTransfer(void) const;
    
  public:
	
		// for the moment, exact match to FILE signatures:
		virtual size_t read(void *ptr, size_t size, size_t nitems)=0;

		virtual size_t write(const void *ptr, size_t size, size_t nitems)=0;

    virtual void close(void);
    
    inline static void close(abstractCommHandle*& h);
		
    void sync(void);
    
    void flush(void);

    // public method: do not use _within_ specializations of read, write:    
    inline bool buffered(void) const;
    
    // alloc or free/release buffer:
    bool buffered(bool flag, size_t bytes=10*1024*1024);    

    // attach to external buffer, don't transfer ownership:
    bool buffered(unsigned char *begin, unsigned char *end);
    
    // attach to external buffers, don't transfer ownership:
    bool buffered(unsigned char *beginRead, unsigned char *endRead,
                  unsigned char *beginWrite, unsigned char *endWrite);
     
    
  	virtual ~abstractCommHandle(void);

};

// global methods:
inline size_t read(void *ptr, size_t size, size_t nitems, abstractCommHandle* h);
inline size_t write(const void *ptr, size_t size, size_t nitems, abstractCommHandle* h);
inline void close(abstractCommHandle*& h);

class fileHandle: public abstractCommHandle{
  private:

		int fd_; // file-descriptor

  protected:
  
	  fileHandle(void); // hide

    // the following are really "flags" governing the behavior of "sync" and "flush"
    //   during buffered data transfer:
    
    // "transferBufferHeader": transfer the buffer header along with the rest of the data
    //    this indicates that the data-size of the buffer portion used will be encoded into the data stream itself:
    virtual bool transferBufferHeader(void) const;  
    
    // "partialBufferTransfer": transfer only whole buffer blocks (or not): some modes (namely MPI BROADCAST) cannot
    //    encode any data-size information, so for these modes the entire buffer must be transfered, even if only a 
    //    portion of it is used (note that this is _independent_ of whether or not the data-size is _also_ encoded on the stream).
    virtual bool partialBufferTransfer(void) const;
		
  public:

    int fd(void)const { return fd_; }

		// for the moment, exact match to FILE signatures:
		virtual size_t read(void *ptr, size_t size, size_t nitems);

		virtual size_t write(const void *ptr, size_t size, size_t nitems);

    virtual void close(void);
    
    // file-position related methods (abstractCommHandle must be a fileHandle):
    void seek(long offset, int whence);
    
    long tell(void);
    
    
    // allocate on heap:
		
    // "open" signature and behavior:
    static fileHandle* open(const char *pathname, int flags, mode_t mode);
    // "fopen" signature and behavior:
    static fileHandle* open(const char *filename, const char *mode);

    // file exists:
		static bool exists(const char *fileName);
		
		// is directory:
		static bool isDirectory(const char *fileName);
		
		// file size:
		static size_t sizeOnMedia(const char *fileName);
				
#if 1
    // from an existing (open) stream:
		static fileHandle* clone(FILE *fp);
     // from an existing (open) file-descriptor:
		static fileHandle* clone(int fd);    
#endif
		
  	virtual ~fileHandle(void);		
};

// global methods:

// "fopen" signature and behavior: 
inline abstractCommHandle* open(const char *filename, const char *mode);
// "open" signature and behavior:
inline abstractCommHandle* open(const char *pathname, int flags, mode_t mode);

inline void seek(abstractCommHandle* fp, long offset, int whence);

inline long tell(abstractCommHandle* fp);



/** The following class "memoryCommHandle" implements a uniform interface to abstractCommHandle to allow reading and writing to in-memory buffers.
 *  The primary purpose of this class is to allow binary format conversion for complicated objects to occur entirely in user or kernel space, 
 *  at the time of _submission_ of an asynchronous-IO request, without the previous requirement for user-kernel transitions occuring at the moment of the
 *  actual I/O action (i.e. when this conversion occured previously in the AIO signal handler).
 */
class memoryCommHandle: public abstractCommHandle{
  private:

    typedef unsigned char __attribute_may_alias__ uchar_alias;
		mutable uchar_alias *buf_, *cp_; // externally managed buffer (i.e. not owned by memoryCommHandle).
    size_t buf_size_;

    // gcc-4.1.2 work-around:  use of std::copy with uchar_alias* causes "symbol already defined" assembler error:
    static void copy_(const uchar_alias* src_b, const uchar_alias* src_e, uchar_alias* dest_b);
    
  protected:
  
	  memoryCommHandle(void); // hide
    
    // the following are really "flags" governing the behavior of "sync" and "flush"
    //   during buffered data transfer:
    
    // "transferBufferHeader": transfer the buffer header along with the rest of the data
    //    this indicates that the data-size of the buffer portion used will be encoded into the data stream itself:
    virtual bool transferBufferHeader(void) const;  
    
    // "partialBufferTransfer": transfer only whole buffer blocks (or not): some modes (namely MPI BROADCAST) cannot
    //    encode any data-size information, so for these modes the entire buffer must be transfered, even if only a 
    //    portion of it is used (note that this is _independent_ of whether or not the data-size is _also_ encoded on the stream).
    virtual bool partialBufferTransfer(void) const;
		
    size_t buf_remaining(void)const;
    
  public:

		// for the moment, exact match to FILE signatures:
		virtual size_t read(void *ptr, size_t size, size_t nitems);

		virtual size_t write(const void *ptr, size_t size, size_t nitems);


    // allocate on heap (local structures, _excepting_ buffer itself):
		static memoryCommHandle* open(void *buffer, size_t buf_size, const char *mode);
		
  	virtual ~memoryCommHandle(void);		
}; 
 


#if defined(__USE_MPI) // ------------------------------------------- MPI ONLY SECTION ------------------------------------------------------------

#if 0
// allow stand-alone header:
namespace MPI;
class MPI::Comm; // forward
#endif

class processCommHandle: public abstractCommHandle{
  public:
  
    enum TRANSFER_KIND {POINT2POINT=0, BROADCAST=1}; 
   
  private:

    #if defined(__CHUNK_TRANSFER__)
    // break transfers into smaller chunks:
    static const size_t CHUNK_SIZE_ = 2*1024;
    #endif
    
		int rank_;
    const MPI::Comm *pComm_;  // pointer: allow stand-alone header
	
	  int tag_;

    TRANSFER_KIND transfer_type_;
    
    static void* MPIsendBuffer_;

#if defined(__USE_PTHREAD) && 0		
		// thread send-receive handoff (mostly for MPI rank-0 multi-threaded singleton):
		static bool signalThreads_, threadDataReady_; 
		static pthread_cond_t threadSignalCond_;
		static pthread_mutex_t threadSignalMutex_;
#endif		
				
	protected:

    // the following are really "flags" governing the behavior of "sync" and "flush"
    //   during buffered data transfer:
    
    // "transferBufferHeader": transfer the buffer header along with the rest of the data
    //    this indicates that the data-size of the buffer portion used will be encoded into the data stream itself:
    virtual bool transferBufferHeader(void) const;  
    
    // "partialBufferTransfer": transfer only whole buffer blocks (or not): some modes (namely MPI BROADCAST) cannot
    //    encode any data-size information, so for these modes the entire buffer must be transfered, even if only a 
    //    portion of it is used (note that this is _independent_ of whether or not the data-size is _also_ encoded on the stream).
    virtual bool partialBufferTransfer(void) const;
		
    processCommHandle(int rank, const MPI::Comm& comm, int tag, TRANSFER_KIND transfer_kind); // hide
		
  public:    
	     
		inline int rank(void)const;
		
		inline const MPI::Comm& comm(void)const;
		
		inline int tag(void)const;
	
    inline TRANSFER_KIND transfer_type(void)const;

		inline TRANSFER_KIND transfer_type(TRANSFER_KIND newTransferType);
    
	  // if incoming tag is MPI::ANY_TAG, receive the tag, 
		//   otherwise _only_ receive a specific tag...
		static bool readTag(abstractCommHandle *h, int& tag);		
    
    // write the specified tag, regardless of whether it is the same as h->tag():
		static bool writeTag(abstractCommHandle *h, int tag);
    
		inline int tag(int newTag);
		
		// for the moment, exact match to FILE signatures:
		virtual size_t read(void *ptr, size_t size, size_t nitems);

		virtual size_t write(const void *ptr, size_t size, size_t nitems);

    // ------ MPI collective communication: --------
		// The following methods implicitely buffer and transfer object U atomically
		//   with respect to the corresponding MPI call:
		// processCommHandle::rank() is target rank for the gather/scatter, _or_
		//   rank of running process is _root_ rank, as determined by signature.
		
		// gather to another rank:
		template <class U>
		void gather(const U& u);
		
		// gather to this rank:
		template <class U>
		void gather(const U& u, std::vector<U>& vU);
		
		// gather to all ranks:
		template <class U>
		void all_gather(const U& u, std::vector<U>& vU);
		
		// scatter from another rank:
		template <class U>
		void scatter(U& u);
		
		// scatter from this rank:
		template <class U>
		void scatter(const std::vector<U>& vU, U& u);
				
		
		// ---------------------------------------------

    // allocate on heap, allow only modes that make sense for processes:
		static processCommHandle* open(const int rank, const MPI::Comm& comm, const char *mode, 
                                    int tag=0, TRANSFER_KIND transfer_type=POINT2POINT);

    static void attachMPISendBuffer(size_t bufSize=32*1024*1024);

		static void detachMPISendBuffer(void);

#if defined(__USE_PTHREAD) && 0		
		inline static bool signalThreads(void);
		
		inline static bool signalThreads(bool flag);
		
		inline static bool threadDataReady(void);
#endif
		
		static bool MPIbuffered(void);
		
  	virtual ~processCommHandle(void);		
};

// ----- global methods: ---------
template <class U>
inline void gather(const U& u, abstractCommHandle *h);

template <class U>
inline void gather(const U& u, std::vector<U>& vU, abstractCommHandle *h);

template <class U>
inline void all_gather(const U& u, std::vector<U>& vU, abstractCommHandle *h);

template <class U>
inline void scatter(U& u, abstractCommHandle *h);

template <class U>
inline void scatter(const std::vector<U>& vU, U& u, abstractCommHandle *h);

#endif // ------------------------------------------- end, MPI ONLY SECTION ------------------------------------------------------------

// **********************************************************************
//
// _generic_ writeBinary, readBinary and binarySize:
//

template<class U>
bool writeBinary(abstractCommHandle *fp, const U& u);

template<class U>
bool readBinary(abstractCommHandle *fp, U& u );

// in consideration of variable-precision number types, the following "binarySize" method requires an instance of class U:
template <class U>
size_t binarySize(const U& u);

// **********************************************************************


// The following string methods are declared _here_ in order to allow this file to be included _first_, prior to other read/write declarations.
// These methods are required _now_ for the _transfer_ of error messages (e.g. between threads or processes):

bool readBinary(abstractCommHandle* h, std::string& s);

bool writeBinary(abstractCommHandle* h, const std::string& s);

size_t binarySize(const std::string& s);

} // namespace commUtil

#include "commUtil_inline.h"

#include "commUtil_template.h"

#endif // __commUtil__h
