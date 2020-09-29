/****************************************************************************
* This product contains one or more programs protected under international
* and U.S. copyright laws as unpublished works.  They are confidential and
* proprietary to Dolby Laboratories.  Their reproduction or disclosure, in
* whole or in part, or the production of derivative works therefrom without
* the express permission of Dolby Laboratories is prohibited.
*
*             Copyright 2011 - 2013 by Dolby Laboratories.
*                     All rights reserved.
****************************************************************************/

/*! @file DeApi.h
    @brief API header file

    This header file defines the interface for API.

    @note

*/

#ifndef DE_API_H
#define DE_API_H

#include "CdeType.h"


/*! @brief kernel handle

    It is a pointer to kernel context, a private structure for kernel.

*/
typedef struct DeKsFlt_t_  *HDeKsFlt_t;

#if defined(c_plusplus) || defined(__cplusplus)
extern "C"
{
#endif

    /*! @brief  process one frame

  These function are the interface to data path.


  DolbyDe takes one input frame, one graphic overlay if supported, together with 
  handle to generate one output frame.
  Both input and output frame buffers must be allocated by caller.
  These are synchronized (blocking) call.

    @param[in] hKs handle ti the kernel
    @param[in] pMeanDe place holder to store DE ITP and DE TP mean score
    @param[in] maxDe place holder to store DE ITP and DE TP Max score
    @param[in] pSdDe place holder to store DE ITP and DE TP Standard deviation score
    @param[in] outBuf Buffer to hold intermediate DE ITP scores for each pixel location.
    @param[in] outBuf1 Buffer to hold intermediate DE TP scores for each pixel location
    @return
        @li 0 Succeed.
        @li <0 Error condition (not implement at moment).
*/

    int DolbyDe(const HDeKsFlt_t hKs,
                float pMeanDe[2], float maxDe[2], float pSdDe[2],
                unsigned char *outBuf, unsigned char *outBuf1);

#ifdef __cplusplus
}
#endif

#endif /* DE_API_H */
