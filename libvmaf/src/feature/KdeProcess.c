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
#include <math.h>  // for pow(,) abs() only
#ifndef WIN32
#include <stdlib.h> // for abs()
#endif
#include <string.h> // for memcpy only

#include "KdeType.h"
#include "KdeUtil.h"
#include "KdeUpSamp.h"
#include "DeApi.h"


#ifndef M_PI
#define M_PI   ((FloatComp_t)3.14159265358979323846)
#endif


//// map xyz to Itp
static void Map2Itp(const HDeKsFlt_t pKs,
                    FloatComp_t *x0, FloatComp_t *x1, FloatComp_t *x2
                   )
{
  FloatComp_t y0, y1, y2;

  M33_T_V3(pKs->m33Xyz2LmsViaR2020, *x0, *x1, *x2, y0, y1, y2);

  L2PQNV3(y0, y1, y2);

  M33_T_V3(pKs->m33Lms2IctcpDm, y0, y1, y2, *x0, *x1, *x2);
}


////// input color space mapping
static void Map2Xyz(const DmKsIMapFlt_t *pKsIMap,
  FloatComp_t *x0, FloatComp_t *x1, FloatComp_t *x2
)
{
  FloatComp_t y0, y1, y2;

  //// YUV =>RGB if necessary
  V3_M_V3(*x0, *x1, *x2,
    pKsIMap->v3Yuv2RgbOff[0], pKsIMap->v3Yuv2RgbOff[1], pKsIMap->v3Yuv2RgbOff[2],
    y0, y1, y2); // x=>y
  M33_T_V3(pKsIMap->m33Yuv2Rgb, y0, y1, y2, *x0, *x1, *x2); // y=>x

  //// normalization
  NV3(pKsIMap->eotfParam.rangeR, pKsIMap->eotfParam.rangeMin, *x0, *x1, *x2);
  CLAMPV3(*x0, *x1, *x2, 0, 1);

  //// processing with normalized input
  if (pKsIMap->eotfParam.eotf == KEotfBt1886) {
    // deGamma
    G2LV3(pKsIMap->eotfParam, *x0, *x1, *x2);
  }
  else if (pKsIMap->eotfParam.eotf == KEotfPower) {
    // deGamma
    P2LV3(pKsIMap->eotfParam, *x0, *x1, *x2);
  }
  // else in PQ
  else if (pKsIMap->eotfParam.eotf == KEotfPq) {
    //// PQ=>L
    PQN2LV3(*x0, *x1, *x2);
  }

  //CLAMPV3(*x0, *x1, *x2, pKsIMap->eotfParam.min, pKsIMap->eotfParam.max); no longer have it

  //// RGB=>white point D65 XYZ: x=>y
  M33_T_V3(pKsIMap->m33Rgb2Xyz, *x0, *x1, *x2, y0, y1, y2);
  V3FIRST_EQ_V3SECOND(*x0, *x1, *x2, y0, y1, y2);
}

#define sind(a)   ((FloatComp_t )sin((a) * M_PI / 180))
#define cosd(a)   ((FloatComp_t )cos((a) * M_PI / 180))

#define SQR(x_)       ((x_)*(x_))
#define POW7(x_, y_)  y_ = (x_)*(x_)*(x_); y_ = y_*y_*(x_)
#define CPOW25_7      (FloatComp_t)(25.0*25.0*25.0*25.0*25.0*25.0*25.0)




int DolbyDe(const HDeKsFlt_t pKs,
            float pMeanDe[2], float maxDe[2], float pSdDe[2],
            unsigned char *outBuf, unsigned char *outBuf1)
{

  // pointer to frame bufs
  FloatComp_t *pIRow[2], *pI[2];
  FloatComp_t *pPRow[2], *pP[2];
  FloatComp_t *pTRow[2], *pT[2];

  int  inIdx, row, col; // loop var
  int byteOff = 0;
  FloatComp_t x0[2], x1[2], x2[2]; // for each of three component

	FloatComp_t de, deTp;
  // stat
  FloatComp_t mnDe = 0, mxDe = 0, sDe = 0;
  FloatComp_t mnDeTp = 0, mxDeTp = 0, sDeTp = 0;
  FloatComp_t *de_buf;
  FloatComp_t *deTp_buf;

  ////// to save some code / cycle, set up row pointer
  for (inIdx = 0; inIdx < 2; ++inIdx) {
    pIRow[inIdx] = pKs->frmBuf0[inIdx] - pKs->rowPitchNum;
    pPRow[inIdx] = pKs->frmBuf1[inIdx] - pKs->rowPitchNum;
    pTRow[inIdx] = pKs->frmBuf2[inIdx] - pKs->rowPitchNum;
  }

  ////// per pixel based process, assuming input are of the same size
  assert(pKs->ksFrmFmt->rowNum == pKs->ksFrmFmt[1].rowNum);
  assert(pKs->ksFrmFmt->colNum == pKs->ksFrmFmt[1].colNum);

  // temp_buff will be used to calculate deitp standard deviation
  memset(outBuf, 0.0, pKs->ksFrmFmt->rowNum * pKs->ksFrmFmt->colNum * sizeof(FloatComp_t));
  memset(outBuf1, 0.0, pKs->ksFrmFmt->rowNum * pKs->ksFrmFmt->colNum * sizeof(FloatComp_t));
  
  de_buf = (FloatComp_t *)outBuf;
  deTp_buf = (FloatComp_t *)outBuf1;

  for (row = 0; row < pKs->ksFrmFmt->rowNum; ++row)
  {
    for (inIdx = 0; inIdx < 2; ++inIdx) {
      pI[inIdx] = pIRow[inIdx] += pKs->rowPitchNum;
      pP[inIdx] = pPRow[inIdx] += pKs->rowPitchNum;
      pT[inIdx] = pTRow[inIdx] += pKs->rowPitchNum;
    }

    for (col = 0; col < pKs->ksFrmFmt->colNum; ++col) {

      //// map to xyz
      for (inIdx = 0; inIdx < 2; ++inIdx) {
        //// get input pxl value
        if (pKs->ksUds[inIdx].chrm != KChrm444) {
          // from up sampled one
          x0[inIdx] = *pI[inIdx]++;
          x1[inIdx] = *pP[inIdx]++;
          x2[inIdx] = *pT[inIdx]++;
        }
        else {
          byteOff = GET_OFFSET(pKs->ksFrmFmt[inIdx], row, col);
         // GET_444_VALUE(pKs->ksFrmFmt[inIdx].dtp, byteOff, inBuf0[inIdx], inBuf1[inIdx], inBuf2[inIdx], x0[inIdx], x1[inIdx], x2[inIdx]);
          GET_444_VALUE(KDtpF32, byteOff, pKs->frmBuf0[inIdx], pKs->frmBuf1[inIdx], pKs->frmBuf2[inIdx], x0[inIdx], x1[inIdx], x2[inIdx]);
        }

        Map2Xyz(pKs->ksIMap + inIdx, x0 + inIdx, x1 + inIdx, x2 + inIdx);
      }


      for (inIdx = 0; inIdx < 2; ++inIdx) {
        //// Xyz=>Itp
        Map2Itp(pKs, x0 + inIdx, x1 + inIdx, x2 + inIdx);
      }

      //// the deItp
        de   = 720 * sqrt(SQR(x0[0] - x0[1]) + SQR(x1[0] - x1[1]) + SQR(x2[0] - x2[1]));
        //deTP considers only color diffrences ignoring intensity variations
        deTp = 720 * sqrt(SQR(x1[0] - x1[1]) + SQR(x2[0] - x2[1]));

        // stats
        mnDe += de; mnDeTp += deTp;
        *de_buf = de; *deTp_buf = deTp;
        de_buf++; deTp_buf++;
        mxDe = MAX(mxDe, de); mxDeTp = MAX(mxDeTp, deTp);

    }
  }
  pMeanDe[0] = (float)(mnDe   * pKs->pxlNumR);
  pMeanDe[1] = (float)(mnDeTp * pKs->pxlNumR);
  maxDe[0] = (float)(mxDe);
  maxDe[1] = (float)(mxDeTp);

  //measure Standard deviation
  de_buf    = (FloatComp_t *)outBuf;
  deTp_buf  = (FloatComp_t *)outBuf1;
  for (row = 0; row < pKs->ksFrmFmt->rowNum; ++row)
  {
    for(col = 0; col < pKs->ksFrmFmt->colNum; ++col)
    {
      sDe += (*de_buf - pMeanDe[0]) * (*de_buf - pMeanDe[0]);
      sDeTp += (*deTp_buf - pMeanDe[1]) * (*deTp_buf - pMeanDe[1]);
      de_buf++; deTp_buf++;
    }
  }
  sDe   = (sDe  * pKs->pxlNumR);
  sDeTp = (sDeTp * pKs->pxlNumR);
  sDe = sqrt(sDe);sDeTp = sqrt(sDeTp);
  pSdDe[0] = sDe;pSdDe[1] = sDeTp;

  return 0;
}

