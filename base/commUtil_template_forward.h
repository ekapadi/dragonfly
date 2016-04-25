#if !defined(__commUtil_template_forward__h)
#define __commUtil_template_forward__h

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


// this file requires prior inclusion of "gmm_ext.h", and possibly "linalgUtil.h"

// $Source: /usr/data0/leipzig_work/tmat_cvs/src/commUtil_template_forward.h,v $

namespace commUtil{

#if defined(__USE_MPI) // ------------------------------------------- MPI ONLY SECTION ------------------------------------------------------------
// gather to another rank:
template <class U>
void processCommHandle::gather(const U& u)
{	
	if (!valid() || !writing()) 
  	throw std::runtime_error("processCommHandle::gather<U>(const U& u): invalid handle or not open for write");

  // allocate the send-buffer:
	const size_t size_(binarySize(u));
	void *buf(new unsigned char[size_]);
	abstractCommHandle *hmem 
	  = memoryCommHandle::open(buf, size_, "wb");
	if (!writeBinary(hmem, u))
	  throw std::runtime_error("processCommHandle::gather<U>: I/O error writing to memory");
	close(hmem);
		
	dynamic_cast<const MPI::Intracomm*>(pComm_)->Gather(buf, size_, MPI::BYTE, NULL, 0, 0, rank());
	
	delete[] reinterpret_cast<unsigned char*>(buf);
}

// gather to this rank:
template <class U>
void processCommHandle::gather(const U& u, std::vector<U>& vU)
{
	if (!valid() || !writing() || !reading()) 
  	throw std::runtime_error("processCommHandle::gather<U>(const U& u, std::vector<U>& vU): invalid handle or not open for read/write");

  // allocate the send and receive buffers:
	const size_t 
	  size_(binarySize(u)),
		procs_(pComm_->Get_size());
	
	void 
	  *sbuf(new unsigned char[size_]),
	  *rbuf(new unsigned char[size_*procs_]);
	abstractCommHandle 
	  *hsend_mem = memoryCommHandle::open(sbuf, size_, "wb"),
	  *hrecv_mem = memoryCommHandle::open(rbuf, size_*procs_, "rb");
	if (!writeBinary(hsend_mem, u))
	  throw std::runtime_error("processCommHandle::gather<U>: I/O error writing to memory");
	close(hsend_mem);
		
	dynamic_cast<const MPI::Intracomm*>(pComm_)->Gather(sbuf, size_, MPI::BYTE, rbuf, size_/*  *procs_ */, MPI::BYTE, pComm_->Get_rank());

  // transfer to destination vector:
	vU.resize(procs_);
	bool status(true);
	for(typename std::vector<U>::iterator itU = vU.begin(), itUEnd = vU.end();
	    status && (itU != itUEnd);
			++itU){
    status = (status && readBinary(hrecv_mem, *itU));
	}
	if (!status)
	  throw std::runtime_error("processCommHandle::gather<U>: I/O error reading from memory");
	close(hrecv_mem);
		
	delete[] reinterpret_cast<unsigned char*>(sbuf);
	delete[] reinterpret_cast<unsigned char*>(rbuf);	
}

// gather to all ranks:
template <class U>
void processCommHandle::all_gather(const U& u, std::vector<U>& vU)
{
	if (!valid() || !writing() || !reading()) 
  	throw std::runtime_error("processCommHandle::all_gather<U>(const U& u, std::vector<U>& vU): invalid handle or not open for read/write");

  // allocate the send and receive buffers:
	const size_t 
	  size_(binarySize(u)),
		procs_(pComm_->Get_size());
	
	void 
	  *sbuf(new unsigned char[size_]),
	  *rbuf(new unsigned char[size_*procs_]);
	abstractCommHandle 
	  *hsend_mem = memoryCommHandle::open(sbuf, size_, "wb"),
	  *hrecv_mem = memoryCommHandle::open(rbuf, size_*procs_, "rb");
	if (!writeBinary(hsend_mem, u))
	  throw std::runtime_error("processCommHandle::all_gather<U>: I/O error writing to memory");
	close(hsend_mem);
		
	dynamic_cast<const MPI::Intracomm*>(pComm_)->Allgather(sbuf, size_, MPI::BYTE, rbuf, size_ /* *procs_ */, MPI::BYTE);

  // transfer to destination vector:
	vU.resize(procs_);
	bool status(true);
	for(typename std::vector<U>::iterator itU = vU.begin(), itUEnd = vU.end();
	    status && (itU != itUEnd);
			++itU){
    status = (status && readBinary(hrecv_mem, *itU));
	}
	if (!status)
	  throw std::runtime_error("processCommHandle::all_gather<U>: I/O error reading from memory");
	close(hrecv_mem);
		
	delete[] reinterpret_cast<unsigned char*>(sbuf);
	delete[] reinterpret_cast<unsigned char*>(rbuf);	
}

// scatter from another rank:
template <class U>
void processCommHandle::scatter(U& u)
{
	if (!valid() || !reading()) 
  	throw std::runtime_error("processCommHandle::scatter<U>(U& u): invalid handle or not open for read");

  // allocate the receive buffer:
	const size_t 
	  size_(binarySize(u));
	
	void 
	  *rbuf(new unsigned char[size_]);
	abstractCommHandle 
	  *hrecv_mem = memoryCommHandle::open(rbuf, size_, "rb");
		
	dynamic_cast<const MPI::Intracomm*>(pComm_)->Scatter(NULL, 0, 0, rbuf, size_, MPI::BYTE, rank());

  // transfer to destination:
  if (!readBinary(hrecv_mem, u))
	  throw std::runtime_error("processCommHandle::scatter<U>(U& u): I/O error reading from memory");
	close(hrecv_mem);
	
	delete[] reinterpret_cast<unsigned char*>(rbuf);	
}

// scatter from this rank:
template <class U>
void processCommHandle::scatter(const std::vector<U>& vU, U& u)
{
	if (!valid() || !writing() || !reading()) 
  	throw std::runtime_error("processCommHandle::scatter<U>(const std::vector<U>& vU, U& u): invalid handle or not open for read/write");

  // allocate the send and receive buffers:
	const size_t 
	  size_(binarySize(u)),
		procs_(pComm_->Get_size());
	
	void 
	  *sbuf(new unsigned char[size_*procs_]),
	  *rbuf(new unsigned char[size_]);
	abstractCommHandle 
	  *hsend_mem = memoryCommHandle::open(sbuf, size_*procs_, "wb"),
	  *hrecv_mem = memoryCommHandle::open(rbuf, size_, "rb");
	bool status(true);
	for(typename std::vector<U>::const_iterator itU = vU.begin(), itUEnd = vU.end();
	    status && (itU != itUEnd);
			++itU){
    status = (status && writeBinary(hsend_mem, *itU));
	}
	if (!status)
	  throw std::runtime_error("processCommHandle::scatter<U>(const std::vector<U>& vU, U& u): I/O error writing to memory");
	close(hsend_mem);
		
	dynamic_cast<const MPI::Intracomm*>(pComm_)->Scatter(sbuf, size_/* *procs_ */, MPI::BYTE, rbuf, size_, MPI::BYTE, pComm_->Get_rank());

	if (!readBinary(hrecv_mem, u))
	  throw std::runtime_error("processCommHandle::scatter<U>(const std::vector<U>& vU, U& u): I/O error reading from memory");
	close(hrecv_mem);
		
	delete[] reinterpret_cast<unsigned char*>(sbuf);
	delete[] reinterpret_cast<unsigned char*>(rbuf);	
}				
#endif // --------------------------------------end:   MPI ONLY SECTION ------------------------------------------------------------	
	
} // namespace commUtil

#endif // __commUtil_template_forward__
