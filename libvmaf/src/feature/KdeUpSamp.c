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
#include <assert.h>
#include "KdeType.h"
#include "KdeUtil.h"
#include "KdeUpSamp.h"



#define Tap_6_FilterUvRowUsHalfSize    3
static const FloatComp_t Tap_6_FilterUvRowUs0_m[Tap_6_FilterUvRowUsHalfSize<<1] = {
  (FloatComp_t)(2/256.0), -(FloatComp_t)(12/256.0), (FloatComp_t)(65/256.0),
  (FloatComp_t)(222/256.0), -(FloatComp_t)(25/256.0), (FloatComp_t)(4/256.0)
};

static const FloatComp_t Tap_6_FilterUvRowUs1_m[Tap_6_FilterUvRowUsHalfSize<<1] = {
  (FloatComp_t)(4/256.0), -(FloatComp_t)(25/256.0), (FloatComp_t)(222/256.0),
  (FloatComp_t)(65/256.0), -(FloatComp_t)(12/256.0), (FloatComp_t)(2/256.0)
};

#define FilterUvColUsHalfSize    4
static const FloatComp_t FilterUvColUs_m[FilterUvColUsHalfSize<<1] = {
  (FloatComp_t)(22/4096.0), (FloatComp_t)(94/4096.0), -(FloatComp_t)(524/4096.0), (FloatComp_t)(2456/4096.0),
  (FloatComp_t)(2456/4096.0), -(FloatComp_t)(524/4096.0), (FloatComp_t)(94/4096.0), (FloatComp_t)(22/4096.0)
};


//// double the column of a single compnnent
static int
ColDoublePlnr(const FloatComp_t *pIn, int rowNum, int colNum, int rowPitchNumIn,
              FloatComp_t minVal, FloatComp_t maxVal,
              int filterUsHalfSize, const FloatComp_t *filterUs_m,
              FloatComp_t *pOut, int rowPitchNumOut
              )
{
  const FloatComp_t *pInRow = pIn;
  const FloatComp_t *pInRowLast;

  const FloatComp_t *filterUs = filterUs_m + filterUsHalfSize;

  FloatComp_t *pOutRow = pOut;

  int r, c, k;
  FloatComp_t fSum;

  for (r = 0; r < rowNum; ++r) {
    pIn = pInRow;
    pOut = pOutRow;

    pInRowLast = pInRow + colNum - 1;

    //// left margin
    for (c = 0; c < filterUsHalfSize - 1; ++c) {
      // copy the first value
      *pOut++ = (FloatComp_t)(*pIn++);

      fSum = 0;
      for (k = -filterUsHalfSize; k < -c -1; ++k) {
        fSum += filterUs[k];
      }
      fSum *= *pInRow;

      for (; k < filterUsHalfSize; ++k) {
        fSum += filterUs[k]*pIn[k];
      }

      // second y is the upsampled one
      *pOut++ = (FloatComp_t)(CLAMP(fSum, minVal, maxVal));
    }

    //// sweet spot
    for (; c < colNum - filterUsHalfSize; ++c) {
      //// upsample y
      // copy the first value
      *pOut++ = (FloatComp_t)(*pIn++);

      fSum = 0;
      for (k = -filterUsHalfSize; k < 0; ++k) {
        fSum += filterUs[k]*(pIn[k] + pIn[-k -1]);
      }

      // second y is the upsampled one
      *pOut++ = (FloatComp_t)(CLAMP(fSum, minVal, maxVal));
    }

    //// right margin
    for (; c < colNum; ++c) {
      //// upsample y
      // copy the first y
      *pOut++ = (FloatComp_t)(*pIn++);

      fSum = 0;
      for (k = -filterUsHalfSize; k < colNum -1 -c; ++k) {
        fSum += filterUs[k]*pIn[k];
      }

      for (; k < filterUsHalfSize; ++k) {
        fSum += filterUs[k]*(*pInRowLast);
      }

      // second y is the upsampled one
      *pOut++ = (FloatComp_t)(CLAMP(fSum, minVal, maxVal));
    }

    pInRow  += rowPitchNumIn;
    pOutRow += rowPitchNumOut;
  }

  return 0;
}

//// convert uyvy to 444 by double U's and V's column
int
UYVYToYUV444InU16(const DmKsUdsFlt_t *pKsUds,
                    const DmKsFrmFmt_t *pFrmFmtUYVY, const unsigned short *pUyVy,
                    FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                    )
{
  const int filterColUsHalfSize = pKsUds->filterUvColUsHalfSize;
  const FloatComp_t *filterColUs_m =  pKsUds->filterUvColUs;

  const int rowNum = pFrmFmtUYVY->rowNum;
  const int colNum = pFrmFmtUYVY->colNum;
  const int rowPitchUyVy = pFrmFmtUYVY->rowPitch;
  const unsigned short *pUyVyRow = pUyVy;
  const unsigned short *pUyVyRow1st[2], *pUyVyRowLast[2];
  const FloatComp_t minVal = pKsUds->minUs, maxVal = pKsUds->maxUs;

  FloatComp_t *pYRow = pY, *pURow = pU, *pVRow = pV;
  FloatComp_t *pUv[2];
  int uvIdx;

  ////// build a filter to align with U/V to improve efficiency
  FloatComp_t filterColUsQ_m[64];
  const int filterColUsHalfSizeQ = 4*filterColUsHalfSize;
  const FloatComp_t *filterColUsQ = filterColUsQ_m + (filterColUsHalfSizeQ - 4);

  int r, c, k;
  FloatComp_t fSum;

  assert(pFrmFmtUYVY->dtp == KDtpU16);

  // build a filter with idx [-(filterColUsHalfSizeQ - 4), filterColUsHalfSizeQ]
  // filterColUsHalfSizeQ instead of filterColUsHalfSizeQ + 3 since the last three is never reached
  // filterColUsQ[k] at k = 0 is filterColUsQ_m[ilterColUsHalfSizeD - 4]
  assert(64 >= 2*filterColUsHalfSizeQ);
  for (k = 0; k < 2*filterColUsHalfSize; ++k) {
    filterColUsQ_m[4*k] = filterColUs_m[k];
    filterColUsQ_m[4*k+1] = 0;
    filterColUsQ_m[4*k+2] = 0;
    filterColUsQ_m[4*k+3] = 0;
  }

  ////// upsample U, V and format conversion
  for (r = 0; r < rowNum; ++r) {
    pUyVy = pUyVyRow;
    pUyVyRow1st[0] = pUyVyRow;      // first U
    pUyVyRow1st[1] = pUyVyRow + 2;  // first V
    pUyVyRowLast[0] = pUyVyRow + 2*colNum - 4; // last U
    pUyVyRowLast[1] = pUyVyRow + 2*colNum - 2; // last V
    pY = pYRow;
    pUv[0] = pURow;
    pUv[1] = pVRow;
    uvIdx = 0;

    //// left boundary condition
    for (c = 0; c < filterColUsHalfSizeQ-4; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k < -c; k += 4) {
        fSum += filterColUsQ[k];
      }
      fSum *= *pUyVyRow1st[uvIdx]; // first U/V

      for (; k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      //// Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    // sweet spot
    for (; c < 2*colNum - filterColUsHalfSizeQ; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      //// Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    // right margine
    for (; c < 2*colNum; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k < 2*colNum - c; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      for (; k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*(*pUyVyRowLast[uvIdx]); // last U/V
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      // Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    pUyVyRow = (const unsigned short *)( (const char *)pUyVyRow + rowPitchUyVy );
    pYRow += rowPitchNum;
    pURow += rowPitchNum;
    pVRow += rowPitchNum;
  }

  return 0;
}

//// double row
static int
RowDoublePlnrInU16(const unsigned short *pIn, int rowNum, int colNum, int rowPitchIn,
                   FloatComp_t minVal, FloatComp_t maxVal,
                   int filterUsHalfSize, const FloatComp_t *filterUs0_m, const FloatComp_t *filterUs1_m,
                   FloatComp_t *pOut, int rowPitchNumOut
                   )
{
  // value at date point dp_ of type dt_ advance by pitch in byte pb_
  #define VAL_AT_ADDR_ADV_BY_PB(dt_, dp_, pb_)  (*( (const dt_ *)( (const char *)(dp_)  + (pb_) ) ))

  const unsigned short *pInRow0 = pIn;

  const FloatComp_t *filterUs0 = filterUs0_m + filterUsHalfSize;
  const FloatComp_t *filterUs1 = filterUs1_m + filterUsHalfSize;

  FloatComp_t *pOutRow0 = pOut;

  ////// vertical upsampling
  int r, c, k, idx0, idx1;
  FloatComp_t fSum0, fSum1;

  for (c = 0; c < colNum; ++c) {
    pIn = pInRow0;
    pOut = pOutRow0;
    //// take care of top bounary
    for (r = 0; r < filterUsHalfSize; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        idx0 = (r + k >= 0) ? k : -r;
        idx1 = (r + k + 1 >= 0) ? k + 1 : -r;

        fSum0 += filterUs0[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, idx0*rowPitchIn);
        fSum1 += filterUs1[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, idx1*rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn  = (const unsigned short *)( (const char *)pIn + rowPitchIn );
    } // r loop for top boundary

    //// take care of center part: no boundary issue
    for (; r < rowNum - filterUsHalfSize; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        fSum0 += filterUs0[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, k*rowPitchIn);
        fSum1 += filterUs1[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, (k+1)*rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn  = (const unsigned short *)( (const char *)pIn + rowPitchIn );
    } // r loop for center of row

    //// take care of bottom bounary
    for (; r < rowNum; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        idx0 = (r + k < rowNum) ? k : rowNum - 1 - r;
        idx1 = (r + k + 1 < rowNum) ? k + 1 : rowNum - 1 - r;

        fSum0 += filterUs0[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, idx0*rowPitchIn);
        fSum1 += filterUs1[k]*VAL_AT_ADDR_ADV_BY_PB(unsigned short, pIn, idx1*rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn  = (const unsigned short *)( (const char *)pIn + rowPitchIn );
    } // r loop for bottom boundary

    ++pInRow0;
    ++pOutRow0;
  } // c loop

  return 0;
}

//// convert data type for a planar plane
static void
FrameAssignPlnrInU16(const unsigned short *pIn, int rowNum, int colNum, int rowPitchIn, FloatComp_t *pOut, int rowPitchNumOut)
{
  const unsigned short *pInRow = pIn;
  FloatComp_t *pOutRow = pOut;
  int r, c;
  for (r = 0; r < rowNum; ++r) {
    pIn = pInRow;
    pOut = pOutRow;
    for (c = 0; c < colNum; ++c) {
      *pOut++ = (FloatComp_t)(*pIn++);
    }
    pInRow  = (const unsigned short *)( (const char *)pInRow  + rowPitchIn );
    pOutRow += rowPitchNumOut;
  }
}

int
YUV420ToYUV444InU16(const DmKsUdsFlt_t *pKsUds,
                    const DmKsFrmFmt_t *pFrmFmt420,
                    const unsigned short *pInBuf0, const unsigned short *pInBuf1, const unsigned short *pInBuf2,
                    FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                    )
{
  assert(pFrmFmt420->dtp == KDtpU16);
  assert(pFrmFmt420->weav == KWeavPlnr);

  // U
  RowDoublePlnrInU16(pInBuf1,
    pFrmFmt420->rowNum>>1, pFrmFmt420->colNum>>1, pFrmFmt420->rowPitchC, pKsUds->minUs, pKsUds->maxUs,
    pKsUds->filterUvRowUsHalfSize, pKsUds->filterUvRowUs0, pKsUds->filterUvRowUs1,
    pY, rowPitchNum);

  ColDoublePlnr(pY,
    pFrmFmt420->rowNum, pFrmFmt420->colNum>>1, rowPitchNum, pKsUds->minUs, pKsUds->maxUs,
    pKsUds->filterUvColUsHalfSize, pKsUds->filterUvColUs,
    pU, rowPitchNum);

  // V
  RowDoublePlnrInU16(pInBuf2,
    pFrmFmt420->rowNum>>1, pFrmFmt420->colNum>>1, pFrmFmt420->rowPitchC, pKsUds->minUs, pKsUds->maxUs,
    pKsUds->filterUvRowUsHalfSize, pKsUds->filterUvRowUs0, pKsUds->filterUvRowUs1,
    pY, rowPitchNum);

  ColDoublePlnr(pY,
    pFrmFmt420->rowNum, pFrmFmt420->colNum>>1, rowPitchNum, pKsUds->minUs, pKsUds->maxUs,
    pKsUds->filterUvColUsHalfSize, pKsUds->filterUvColUs,
    pV, rowPitchNum);

  // Y
  FrameAssignPlnrInU16(pInBuf0,
    pFrmFmt420->rowNum, pFrmFmt420->colNum, pFrmFmt420->rowPitch,
    pY, rowPitchNum);

  return 0;
}

//// convert uyvy to 444 by double U's and V's column
int
UYVYToYUV444InU8(const DmKsUdsFlt_t *pKsUds,
                    const DmKsFrmFmt_t *pFrmFmtUYVY, const unsigned char *pUyVy,
                    FloatComp_t * pY, FloatComp_t * pU, FloatComp_t * pV, int rowPitchNum
                    )
{
  const int filterColUsHalfSize = pKsUds->filterUvColUsHalfSize;
  const FloatComp_t *filterColUs_m =  pKsUds->filterUvColUs;

  const int rowNum = pFrmFmtUYVY->rowNum;
  const int colNum = pFrmFmtUYVY->colNum;
  const int rowPitchUyVy = pFrmFmtUYVY->rowPitch;
  const unsigned char *pUyVyRow = pUyVy;
  const unsigned char *pUyVyRow1st[2], *pUyVyRowLast[2];
  const FloatComp_t minVal = pKsUds->minUs, maxVal = pKsUds->maxUs;

  FloatComp_t *pYRow = pY, *pURow = pU, *pVRow = pV;
  FloatComp_t *pUv[2];
  int uvIdx;

  ////// build a filter to align with U/V to improve efficiency
  FloatComp_t filterColUsQ_m[64];
  const int filterColUsHalfSizeQ = 4*filterColUsHalfSize;
  const FloatComp_t *filterColUsQ = filterColUsQ_m + (filterColUsHalfSizeQ - 4);

  int r, c, k;
  FloatComp_t fSum;

  assert(pFrmFmtUYVY->dtp == KDtpU8);

  // build a filter with idx [-(filterColUsHalfSizeQ - 4), filterColUsHalfSizeQ]
  // filterColUsHalfSizeQ instead of filterColUsHalfSizeQ + 3 since the last three is never reached
  // filterColUsQ[k] at k = 0 is filterColUsQ_m[ilterColUsHalfSizeD - 4]
  assert(64 >= 2*filterColUsHalfSizeQ);
  for (k = 0; k < 2*filterColUsHalfSize; ++k) {
    filterColUsQ_m[4*k] = filterColUs_m[k];
    filterColUsQ_m[4*k+1] = 0;
    filterColUsQ_m[4*k+2] = 0;
    filterColUsQ_m[4*k+3] = 0;
  }

  ////// upsample U, V and format conversion
  for (r = 0; r < rowNum; ++r) {
    pUyVy = pUyVyRow;
    pUyVyRow1st[0] = pUyVyRow;      // first U
    pUyVyRow1st[1] = pUyVyRow + 2;  // first V
    pUyVyRowLast[0] = pUyVyRow + 2*colNum - 4; // last U
    pUyVyRowLast[1] = pUyVyRow + 2*colNum - 2; // last V
    pY = pYRow;
    pUv[0] = pURow;
    pUv[1] = pVRow;
    uvIdx = 0;

    //// left boundary condition
    for (c = 0; c < filterColUsHalfSizeQ-4; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k < -c; k += 4) {
        fSum += filterColUsQ[k];
      }
      fSum *= *pUyVyRow1st[uvIdx]; // first U/V

      for (; k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      //// Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    // sweet spot
    for (; c < 2*colNum - filterColUsHalfSizeQ; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      //// Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    // right margine
    for (; c < 2*colNum; c += 2) {
      //// U/V
      // copy for the even one at planar
      *pUv[uvIdx]++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSizeQ - 4); k < 2*colNum - c; k += 4) {
          fSum += filterColUsQ[k]*pUyVy[k];
      }
      for (; k <= filterColUsHalfSizeQ; k += 4) {
          fSum += filterColUsQ[k]*(*pUyVyRowLast[uvIdx]); // last U/V
      }
      *pUv[uvIdx]++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      // Y: copy
      ++pUyVy; // pass U/V to Y
      *pY++ = (FloatComp_t)(*pUyVy++); // pass Y to V/U

      uvIdx ^= 1;
    }

    pUyVyRow = (const unsigned char *)( (const char *)pUyVyRow + rowPitchUyVy );
    pYRow += rowPitchNum;
    pURow += rowPitchNum;
    pVRow += rowPitchNum;
  }

  return 0;
}

//// double row
static int
RowDoublePlnrInU8(const FloatComp_t *pIn, int rowNum, int colNum, int rowPitchIn,
                  FloatComp_t minVal, FloatComp_t maxVal,
                  int filterUsHalfSize, const FloatComp_t *filterUs0_m, const FloatComp_t *filterUs1_m,
                  FloatComp_t *pOut, int rowPitchNumOut)
{
  // value at date point dp_ of type dt_ advance by pitch in byte pb_
#define VAL_AT_ADDR_ADV_BY_PB(dt_, dp_, pb_) (*((const dt_ *)((const dt_ *)(dp_) + (pb_))))

  const FloatComp_t *pInRow0 = pIn;

  const FloatComp_t *filterUs0 = filterUs0_m + filterUsHalfSize;
  const FloatComp_t *filterUs1 = filterUs1_m + filterUsHalfSize;

  FloatComp_t *pOutRow0 = pOut;

  ////// vertical upsampling
  int r, c, k, idx0, idx1;
  FloatComp_t fSum0, fSum1;

  for (c = 0; c < colNum; ++c) {
    pIn = pInRow0;
    pOut = pOutRow0;
    //// take care of top bounary
    for (r = 0; r < filterUsHalfSize; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        idx0 = (r + k >= 0) ? k : -r;
        idx1 = (r + k + 1 >= 0) ? k + 1 : -r;

        fSum0 += filterUs0[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, idx0 * rowPitchIn);
        fSum1 += filterUs1[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, idx1 * rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn = (const FloatComp_t *)((FloatComp_t *)pIn + rowPitchIn);
    } // r loop for top boundary

    //// take care of center part: no boundary issue
    for (; r < rowNum - filterUsHalfSize; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        fSum0 += filterUs0[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, k * rowPitchIn);
        fSum1 += filterUs1[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, (k + 1) * rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn = (const FloatComp_t *)((FloatComp_t *)pIn + rowPitchIn);
    } // r loop for center of row

    //// take care of bottom bounary
    for (; r < rowNum; ++r) {
      fSum0 = 0; fSum1 = 0;
      for (k = -filterUsHalfSize; k < filterUsHalfSize; ++k) {
        // row to frame
        idx0 = (r + k < rowNum) ? k : rowNum - 1 - r;
        idx1 = (r + k + 1 < rowNum) ? k + 1 : rowNum - 1 - r;

        fSum0 += filterUs0[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, idx0 * rowPitchIn);
        fSum1 += filterUs1[k] * VAL_AT_ADDR_ADV_BY_PB(FloatComp_t, pIn, idx1 * rowPitchIn);
      }
      *pOut = (FloatComp_t)(CLAMP(fSum0, minVal, maxVal)); // even one
      pOut += rowPitchNumOut;

      *pOut = (FloatComp_t)(CLAMP(fSum1, minVal, maxVal));// old one
      pOut += rowPitchNumOut;

      pIn = (const FloatComp_t *)((FloatComp_t *)pIn + rowPitchIn);
    } // r loop for bottom boundary

    ++pInRow0;
    ++pOutRow0;
  } // c loop

  return 0;
}

//// convert data type for a planar plane
static void
FrameAssignPlnrInU8(const unsigned char *pIn, int rowNum, int colNum, int rowPitchIn, FloatComp_t *pOut, int rowPitchNumOut)
{
  const unsigned char *pInRow = pIn;
  FloatComp_t *pOutRow = pOut;
  int r, c;
  for (r = 0; r < rowNum; ++r) {
    pIn = pInRow;
    pOut = pOutRow;
    for (c = 0; c < colNum; ++c) {
      *pOut++ = (FloatComp_t)(*pIn++);
    }
    pInRow  = (const unsigned char *)( (const char *)pInRow  + rowPitchIn );
    pOutRow += rowPitchNumOut;
  }
}

int YUV420ToYUV444InU8(FloatComp_t *ref_data_cb, FloatComp_t *ref_data_cr, FloatComp_t *temp_data, int src_width, int src_height, int dst_stride, int bdp)
{
  // U
  RowDoublePlnrInU8(ref_data_cb,
    src_height>>1, src_width>>1, src_width>>1, 0, (FloatComp_t)((1 << bdp) - 1),
    Tap_6_FilterUvRowUsHalfSize, (FloatComp_t *)Tap_6_FilterUvRowUs0_m, (FloatComp_t *)Tap_6_FilterUvRowUs1_m,
    temp_data, dst_stride);

  ColDoublePlnr(temp_data,
    src_height, src_width>>1, src_width, 0, (FloatComp_t)((1 << bdp) - 1),
    FilterUvColUsHalfSize, (FloatComp_t *)FilterUvColUs_m,
    ref_data_cb, dst_stride);

  RowDoublePlnrInU8(ref_data_cr,
    src_height>>1, src_width>>1, src_width>>1, 0, (FloatComp_t)((1 << bdp) - 1),
    Tap_6_FilterUvRowUsHalfSize, (FloatComp_t *)Tap_6_FilterUvRowUs0_m, (FloatComp_t *)Tap_6_FilterUvRowUs1_m,
    temp_data, dst_stride);

  ColDoublePlnr(temp_data,
    src_height, src_width>>1, src_width, 0, (FloatComp_t)((1 << bdp) - 1),
    FilterUvColUsHalfSize, (FloatComp_t *)FilterUvColUs_m,
    ref_data_cr, dst_stride);

  return 0;
}

// convert Chroma component to YUV444 from YUV 422p 
// Converts only one component at a time
int YUV422ToYUV444InU8(FloatComp_t *ref_data_cb, FloatComp_t *temp_data, int src_width, int src_height, int dst_stride, int bdp)
{
  const int filterColUsHalfSize = FilterUvColUsHalfSize;
  const FloatComp_t *filterColUs_m = (FloatComp_t *)FilterUvColUs_m;

  const int rowNum = src_height;
  const int colNum = src_width;
  const int rowPitchUyVy = dst_stride;
  const FloatComp_t *pUyVyRow = temp_data, *pUyVy;
  const FloatComp_t *pUyVyRow1st, *pUyVyRowLast;
  const FloatComp_t minVal = 0, maxVal = (FloatComp_t)((1 << bdp) - 1);

  FloatComp_t *pURow = ref_data_cb;
  FloatComp_t *pUv;


  const FloatComp_t *filterColUsQ = filterColUs_m + (filterColUsHalfSize - 1);

  int r, c, k;
  FloatComp_t fSum;

  ////// upsample U
  for (r = 0; r < rowNum; ++r)
  {
    pUyVy = pUyVyRow;
    pUyVyRow1st = pUyVyRow;                   // first U
    pUyVyRowLast = pUyVyRow + (colNum>>1) - 1; // last U

    pUv = pURow;

    //// left boundary condition
    for (c = 0; c < filterColUsHalfSize - 1; c += 1)
    {
      //// U/V
      // copy for the even one at planar
      *pUv++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSize - 1); k < -c; k += 1)
      {
        fSum += filterColUsQ[k];
      }
      fSum *= *pUyVyRow1st; // first U/V

      for (; k <= filterColUsHalfSize; k += 1)
      {
        fSum += filterColUsQ[k] * pUyVy[k];
      }
      *pUv++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

       ++pUyVy;                         // pass U/V to Y
     }

    // sweet spot
     for (; c < (colNum >> 1) - filterColUsHalfSize; c += 1)
     {
       //// U/V
       // copy for the even one at planar
       *pUv++ = (FloatComp_t)(*pUyVy);

       // upsample for the odd one at planar
       fSum = 0;
       for (k = -(filterColUsHalfSize - 1); k <= filterColUsHalfSize; k += 1)
       {
         fSum += filterColUsQ[k] * pUyVy[k];
       }
       *pUv++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

       ++pUyVy;                         // pass U/V to Y

    }

    // right margine
    for (; c < (colNum >> 1); c += 1)
    {
      //// U/V
      // copy for the even one at planar
      *pUv++ = (FloatComp_t)(*pUyVy);

      // upsample for the odd one at planar
      fSum = 0;
      for (k = -(filterColUsHalfSize - 1); k < (colNum >> 1) - c; k += 1)
      {
        fSum += filterColUsQ[k] * pUyVy[k];
      }
      for (; k <= filterColUsHalfSize; k += 1)
      {
        fSum += filterColUsQ[k] * (*pUyVyRowLast); // last U/V
      }
      *pUv++ = (FloatComp_t)CLAMP(fSum, minVal, maxVal);

      // Y: copy
      ++pUyVy;                         // pass U/V to Y
    }

    pUyVyRow = pUyVyRow + (colNum >> 1);

    pURow += dst_stride;
  }

  return 0;
}
