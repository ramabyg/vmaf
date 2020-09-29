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
#ifndef K_DE_TYPE_FLT_H
#define K_DE_TYPE_FLT_H

#include "CdeType.h"

// up sampling need to know the size and type of data for filtering,
// without template, we have to work on a specific type.
// to avoid duplicated code, we assume the following type for input 420 and 422 data type

//// for now, these table size/array/matrix(in fact the pre-allocated max size) are fixed
// for ms
#define MS_GUSSIAN_FLTR_RADIUS_MAX    5

typedef unsigned char KDtp_t;
typedef unsigned char KClr_t;
typedef unsigned char KChrm_t;
typedef unsigned char KWeav_t;
typedef unsigned char KEotf_t;
typedef unsigned char KLoc_t;
typedef unsigned char KBdp_t;


////// pixel definition, force to unsigned char
// bit depth is of type unsigned char, numerical value.
// MSB aliagned will take the whole bits of an interger data type
#define KDtpU16    ((KDtp_t)0)
#define KDtpU8    ((KDtp_t)1)
#define KDtpF32   ((KDtp_t)2)
#define KDtpNum   ((KDtp_t)3)

#define KClrYuv   ((KClr_t)0)
#define KClrRgb   ((KClr_t)1)
#define KClrRgba  ((KClr_t)2)

#define KClrIpt   ((KClr_t)3)
#define KClrNum   ((KClr_t)4)

#define KChrm420  ((KChrm_t)0)
#define KChrm422  ((KChrm_t)1)
#define KChrm444  ((KChrm_t)2)
#define KChrmNum  ((KChrm_t)3)

#define KWeavPlnr ((KWeav_t)0)
#define KWeavIntl ((KWeav_t)1) // per component interleave, YUVYUV..., RGBRGB...
#define KWeavUyVy ((KWeav_t)2) // UYVY and UYVY only
#define KWeavNum  ((KWeav_t)3)

#define KEotfBt1886   ((KEotf_t)0)
#define KEotfPq       ((KEotf_t)1)
#define KEotfPower    ((KEotf_t)2)
#define KEotfHlg      ((KEotf_t)3)
#define KEotfNum      ((KEotf_t)4)

#define KLocHost    ((KLoc_t)0)
#define KLocDev     ((KLoc_t)1)
#define KLocNum     ((KLoc_t)2)

//// the io frame format
typedef struct DmKsFrmFmt_t_
{
  int rowNum, colNum;  // size
  KDtp_t dtp;
  KWeav_t weav;
  KLoc_t loc;

  //// that for interleaved one or the first component, and alpha of planar one
  int rowPitch;       // row pitch in byte
  int colPitch;       // col pitch in byte: a componet size for planar, whole pixel size for interleaved
  //// that for the second,third(likely Chroma) in case of planar layout
  int rowPitchC;      // row pitch in byte
} DmKsFrmFmt_t;

////// data path/kernel structure: KS
// for pre-processing: up sampling
typedef struct DmKsUdsFlt_t_
{
  //// for UV
  KChrm_t chrm;

  FloatComp_t minUs, maxUs;

  int filterUvRowUsHalfSize;
  FloatComp_t *filterUvRowUs0;
  FloatComp_t *filterUvRowUs1;

  int filterUvColUsHalfSize;
  FloatComp_t *filterUvColUs;
} DmKsUdsFlt_t;

//// E=>O
typedef struct EotfParamFlt_t_
{
  // range
  FloatComp_t rangeMin, rangeR;
  // in linear orgb space
  FloatComp_t min, max;
  // gamma stuff
  KEotf_t eotf;  // KEotfBt1886, KEotfPq, KEotfPower coded
  FloatComp_t a, b, gamma;
} EotfParamFlt_t;

//// input mapping
typedef struct DmKsIMapFlt_t_ {
  // YUV(rgb, ipt, itp)=>RGB(lms, xyz)
  FloatComp_t m33Yuv2Rgb[3][3];
  FloatComp_t v3Yuv2RgbOff[3];
  // RGB(lms, xyz)=>XYZ
  EotfParamFlt_t eotfParam;
  FloatComp_t m33Rgb2Xyz[3][3];
  // reciprocal of white point: de2000 only
  FloatComp_t   ksWpR[3];
} DmKsIMapFlt_t;

//// the ks collection
typedef struct DeKsFlt_t_
{
  //// tmp buf to store (I, P, T), post up sample or pre down sample (Y, U, V)
  int rowPitchNum;   // pitch in index
  FloatComp_t *frmBuf0[2], *frmBuf1[2], *frmBuf2[2]; // [0]/[1] ref/test frame

  //// that for xfer to Lab
  // [0]/[1] ref/test frame
  DmKsUdsFlt_t  ksUds[2];
  DmKsFrmFmt_t  ksFrmFmt[2];
  DmKsIMapFlt_t ksIMap[2];

  FloatComp_t pxlNumR; /* total pxl number for calculate mean */

  // deItp only
  FloatComp_t m33Xyz2LmsViaR2020[3][3];
  FloatComp_t m33Lms2IctcpDm[3][3];

  CEotf_t eotf;
  FloatComp_t gamma;
  CClr_t clr_fmt;
  CRng_t Rng;
  int bpc;
  int Yuv2RgbExt;
  CYuvXferSpec_t YuvXferSpec;
  CRgbDef_t RgbDef; // color space def for RGB<=>LMS 
  int Yuv2RgbOffNormExt;
  int Yuv2RgbOffExt;
  int Rgb2LmsM33Ext;
  int Rgb2LmsRgbwExt;
  int NonStdItp; 

} DmKsFlt_t;

#endif // K_DE_TYPE_FLT_H
