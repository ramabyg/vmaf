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
#ifndef K_DE_UP_DOWN_SAMP_FLT_H
#define K_DE_UP_DOWN_SAMP_FLT_H

#include "KdeType.h"

extern int
UYVYToYUV444InU16(const DmKsUdsFlt_t *pKsUds,
                  const DmKsFrmFmt_t *pFrmFmtUYVY, const unsigned short *pUyVy,
                  FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                 );

extern int
YUV420ToYUV444InU16(const DmKsUdsFlt_t *pKsUds,
                    const DmKsFrmFmt_t *pFrmFmt420,
                    const unsigned short *pInBuf0, const unsigned short *pInBuf1, const unsigned short *pInBuf2,
                    FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                   );
extern int
UYVYToYUV444InU8(const DmKsUdsFlt_t *pKsUds,
                 const DmKsFrmFmt_t *pFrmFmtUYVY, const unsigned char *pUyVy,
                 FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                );
#if 0
extern int
YUV420ToYUV444InU8(const DmKsUdsFlt_t *pKsUds,
                   const DmKsFrmFmt_t *pFrmFmt420,
                   const unsigned char *pInBuf0, const unsigned char *pInBuf1, const unsigned char *pInBuf2,
                   FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                  );
#else
extern int
YUV420ToYUV444InU8(FloatComp_t * ref_data_cb, FloatComp_t * ref_data_cr,
                   FloatComp_t * temp_data, int src_width,
                   int src_height, int dst_stride, int bdp);
extern int 
YUV422ToYUV444InU8(FloatComp_t *ref_data_cb, FloatComp_t *temp_data, 
                  int src_width, int src_height, int dst_stride, int bdp);
#endif

#endif // K_DE_UP_DOWN_SAMP_FLT_H
