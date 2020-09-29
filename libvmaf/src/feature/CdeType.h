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
#ifndef C_DE_TYPE_FLT_H
#define C_DE_TYPE_FLT_H

#include <stdint.h>

typedef float FloatComp_t;

/*! @brief configuration parameter structure

    configuration/Ctrol layer parameters

*/
typedef enum {
  CPlatformCpu = 0,
  CPlatformCuda,
  CPlatformOpenCl,
  CPlatformNum
} CPlatform_t;

typedef enum {
  CMsMethodDbEdge = 0,
  CMsMethodOff,
  CMsMethodNum
} CMsMethod_t;

////// pixel definition, force to unsigned char
typedef enum {
  CClrYuv = 0,
  CClrRgb,
  CClrRgba,
  CClrIpt,
  CClrICtCp,
  CClrGrey,
  CClrNum
} CClr_t;

typedef enum {
  CChrm420 = 0,
  CChrm422,
  CChrm444,
  CChrmNum
} CChrm_t;

typedef enum {
  CDtpU16 = 0,
  CDtpU8,
  CDtpF32,
//PxlDtpF16,
  CDtpNum
} CDtp_t;

typedef enum {
  CLocHost = 0,
  CLocDev,
  CLocNum
} CLoc_t;

typedef enum {
  CWeavPlnr = 0,
  CWeavIntl, // per component interleave, YUVYUV..., RGBRGB...
  CWeavUyVy, // UYVY and UYVY only
  CWeavNum
} CWeav_t;

typedef enum {
  CEotfBt1886 = 0,
  CEotfPq,
  CEotfPower,
  CEotfHlg,
  CEotfNum = 4+1, // to make CEotfNum larger by 1 to leave room for EotfGamma in CLI design
} CEotf_t;

typedef enum {
  CRngNarrow = 0, // head
  CRngFull = 1,   // will be the in data type(bits) range
  CRngSdi = 2,    // pq
  CRngNum
} CRng_t;

typedef enum {
  CWpDefD65 = 0,
  CWpDefD60,
  CWpDefNum
} CWpDef_t;

typedef enum {
  CRgbDefP3d65 = 0,
  CRgbDefDci,
  CRgbDefR709,
  CRgbDefR2020,
  CRgbDefAces,
  CRgbDefAlexa,

  CRgbDefNum
} CRgbDef_t;

typedef enum {
  CYuvXferSpecP3d65  = CRgbDefP3d65,
  CYuvXferSpecDci    = CRgbDefDci,
  CYuvXferSpecR709   = CRgbDefR709,
  CYuvXferSpecR2020  = CRgbDefR2020,
  CYuvXferSpecAces   = CRgbDefAces,
  CYuvXferSpecAlexa  = CRgbDefAlexa,
  CYuvXferSpecR601,

  CYuvXferSpecNum
} CYuvXferSpec_t;

#define GET_PXL_COMP_BYTES(cDtp) (  ( (cDtp) == CDtpU16 ) ? 2 :  \
                                    ( (cDtp) == CDtpU8 )  ? 1 :  \
                                    ( (cDtp) == CDtpF32 ) ? 4 :  \
                                    2 )

#define GET_PXL_ECMP_NUM(cClr, cChrm)    ( ( (cClr)  == CClrRgba ) ? 4 : \
                                           ( (cClr)  == CClrGrey ) ? 1 : \
                                           ( (cChrm) == CChrm422 ) ? 2 : \
                                           3  )
// UYVY is interleaved as UY at even, VY at odd, treat as 2 effective component

#define GET_PXL_COL_PITCH(cClr, cChrm, cWeav, cDtp)  (  (cWeav == CWeavPlnr) ?            \
    GET_PXL_COMP_BYTES(cDtp) : GET_PXL_ECMP_NUM(cClr, cChrm) * GET_PXL_COMP_BYTES(cDtp)  )

#endif // C_DE_TYPE_FLT_H
