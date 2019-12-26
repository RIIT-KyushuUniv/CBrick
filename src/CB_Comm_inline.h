#ifndef _CB_COMM_INLINE_H_
#define _CB_COMM_INLINE_H_

/*
###################################################################################
#
# CBrick
#
# Copyright (c) 2017-2019 Research Institute for Information Technology(RIIT),
#                    Kyushu University.  All rights reserved.
#
####################################################################################
*/

/**
 * @file   CB_Comm_inline.h
 * @brief  BrickComm class inline Header
 */

#include <typeinfo>

// #############################################################
// MPI_Datatypeを取得
template<class T> inline
MPI_Datatype BrickComm::GetMPI_Datatype( T *ptr )
{
  if( !ptr )
  {
    return MPI_DATATYPE_NULL;
  }
  
  //  if     ( typeid(ptr) == typeid(char*)  )              return MPI_CHAR;
  //  else if( typeid(ptr) == typeid(short*) )              return MPI_SHORT;
  if     ( typeid(ptr) == typeid(int*) )                return MPI_INT;
  else if( typeid(ptr) == typeid(long*) )               return MPI_LONG;
  else if( typeid(ptr) == typeid(float*) )              return MPI_FLOAT;
  else if( typeid(ptr) == typeid(double*) )             return MPI_DOUBLE;
  //  else if( typeid(ptr) == typeid(long double*) )        return MPI_LONG_DOUBLE;
  //  else if( typeid(ptr) == typeid(unsigned char*) )      return MPI_UNSIGNED_CHAR;
  //  else if( typeid(ptr) == typeid(unsigned short*) )     return MPI_UNSIGNED_SHORT;
  else if( typeid(ptr) == typeid(unsigned*) )           return MPI_UNSIGNED;
  else if( typeid(ptr) == typeid(unsigned int*) )       return MPI_UNSIGNED;
  else if( typeid(ptr) == typeid(unsigned long*) )      return MPI_UNSIGNED_LONG;
  //#ifdef MPI_LONG_LONG_INT
  //  else if( typeid(ptr) == typeid(long long int*) )      return MPI_LONG_LONG_INT;
  //#endif
#ifdef MPI_LONG_LONG
   else if( typeid(ptr) == typeid(long long*) )          return MPI_LONG_LONG;
#endif
  //#ifdef MPI_UNSIGNED_LONG_LONG
  //  else if( typeid(ptr) == typeid(unsigned long long*) ) return MPI_UNSIGNED_LONG_LONG;
  //#endif
  else {
    printf("CBrick error : MPI_DATATYPE_NULL\n");
    exit(-1);
  }
  
  return MPI_DATATYPE_NULL;
}


// #############################################################
// 隣接間通信 IsendIrecv Wrapper Interface 
template <class T> inline
bool BrickComm::IsendIrecv(T* ms,
                           T* mr,
                           T* ps,
                           T* pr,
                           int msz,
                           int nIDm,
                           int nIDp,
                           MPI_Request *req)
{
  if( !ms || !mr || !ps || !pr ) return false;
  
  MPI_Datatype dtype = BrickComm::GetMPI_Datatype(ms);
  if( dtype == MPI_DATATYPE_NULL ) return false;
  
  return IsendIrecv(dtype,
                    (void*)ms, (void*)mr, (void*)ps, (void*)pr,
                    msz, nIDm, nIDp, req);
}


// #############################################################
// Irecv Wrapper Interface
template <class T> inline
bool BrickComm::IrecvData(T* ptr,
                          int sz,
                          int nID,
                          MPI_Request *req)
{
  if( !ptr ) return false;
  
  MPI_Datatype dtype = BrickComm::GetMPI_Datatype(ptr);
  if( dtype == MPI_DATATYPE_NULL ) return false;
  
  return IrecvData(dtype, (void*)ptr, sz, nID, req);
}


// #############################################################
// Isend Wrapper Interface
template <class T> inline
bool BrickComm::IsendData(T* ptr,
                          int sz,
                          int nID,
                          MPI_Request *req)
{
  if( !ptr ) return false;
  
  MPI_Datatype dtype = BrickComm::GetMPI_Datatype(ptr);
  if( dtype == MPI_DATATYPE_NULL ) return false;
  
  return IsendData(dtype, (void*)ptr, sz, nID, req);
}



#endif // _CB_COMM_INLINE_H_
