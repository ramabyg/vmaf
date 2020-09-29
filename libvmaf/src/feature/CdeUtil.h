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
#ifndef C_DE_UTIL_FLT_H
#define C_DE_UTIL_FLT_H

//#include <stdio.h>
#include "CdeType.h"
#include "KdeType.h"

#define MAX2S(a_, b_)  (((a_) >= (b_)) ? (a_) : (b_))
#define MIN2S(a_, b_)  (((a_) <= (b_)) ? (a_) : (b_))
#define CLAMPS(a_, mn_, mx_)  ( ((a_) >= (mn_)) ? ( ((a_) <= (mx_)) ? (a_) : (mx_) ) : (mn_) )

#define IS_CCLR_IPT_ICTCP(cclr_)  ( (cclr_) == CClrIpt || (cclr_) == CClrICtCp )
#define IS_CCLR_RGB_RGBA(cclr_)   ( (cclr_) == CClrRgb || (cclr_) == CClrRgba )
#define IS_CCLR_COMPONENT(cclr_)  (!IS_CCLR_RGB_RGBA(cclr_))

#define IS_CT4(m33InvCt_) ( 0 > 3*(int)(m33InvCt_)[0][0] + 97*(int)(m33InvCt_)[0][1] )

#define ZERO_V3(v3_)        { (v3_)[0] = 0; (v3_)[1] = 0; (v3_)[2] = 0; }
#define ONE_V3(v3_)         { (v3_)[0] = 1; (v3_)[1] = 1; (v3_)[2] = 1; }

#define S_TIMES_V3(s_, v3_) { (v3_)[0] *= (s_); (v3_)[1] *= (s_); (v3_)[2] *= (s_); }
#define S_TIMES_M33(s_, m33_) { S_TIMES_V3(s_, (m33_)[0]); S_TIMES_V3(s_, (m33_)[1]); S_TIMES_V3(s_, (m33_)[2]); }

#define DIAG_OTHER_M33(diag_, othr_, m33_) \
{ (m33_)[0][0] = (diag_); (m33_)[0][1] = (othr_); (m33_)[0][2] = (othr_); \
  (m33_)[1][0] = (othr_); (m33_)[1][1] = (diag_); (m33_)[1][2] = (othr_); \
  (m33_)[2][0] = (othr_); (m33_)[2][1] = (othr_); (m33_)[2][2] = (diag_); }

#define IDENTITY_M33(m33_) DIAG_OTHER_M33(1, 0, m33_)

#define DE_CROSSTALK_M33(ct_, m33_) DIAG_OTHER_M33((1-(ct_))/(1-3*(ct_)), -(ct_)/(1- 3*(ct_)), m33_)

/* de-crosstalk to the left matrix */
#define DE_CT2LEFT(m33In_, ct_, m33Ot_, ty_)  \
{ double rowFct[3], d1m3ctR = 1/(1 - 3*(ct_));\
  int i_, j_;                                   \
  for (i_ = 0; i_ < 3; ++i_) rowFct[i_] = (ct_)*((m33In_)[i_][0] + (m33In_)[i_][1] + (m33In_)[i_][2]); \
  for (i_ = 0; i_ < 3; ++i_) for (j_ = 0; j_ < 3; ++j_) (m33Ot_)[i_][j_] = (ty_)(((m33In_)[i_][j_] - rowFct[i_])*d1m3ctR); }


#if defined(c_plusplus) || defined(__cplusplus)
extern "C"
{
#endif

// to avoid c lib memcpy
void MemCpyByte(void *pDst, const void *pSrc, int len);
// to avoid c lib memcmp
int MemEqByte(const void *pD1, const void *pD2, int len);

// assign matrix and vector(to kernel: data type is FloatComp_t)
void AssignM33D2(const double m33In[3][3], FloatComp_t m33Out[3][3]);
void AssignV3D2(const double v3In[3], FloatComp_t v3Out[3]);

// PQ normalized, assuming input is within well defined range
double L2PqNorm (double L);
double PqNorm2L (double PQ);

void GetEotfParams(double lMin, double lMax, double gamma, double *pA, double *pB);

void GetYuvRgbOff(CRng_t rng, int bits, double v3YuvRgbOff[3]);
void GetYuv2RgbM33(CYuvXferSpec_t yuvXferSpec, double m33Yuv2Rgb[3][3]);
void GetYuv2RgbM33Narrow(CYuvXferSpec_t yuvXferSpec, int bits, double m33Yuv2Rgb[3][3]);

void GetRgb2LmsByDefM33(CRgbDef_t rgbDef, double m33Rgb2Lms[3][3]);
void GetRgb2LmsByPrimsM33(double rx, double ry, double gx, double gy, double bx, double by,
                          double wx, double wy, double m33Rgb2Lms[3][3]);

void GetRgb2LmsCtWpM33(const double m33[3][3], double crossTalk, const double v3Wp[3], double m33Rgb2Lms[3][3]);

void GetRgb2XyzFrom2LmsM33(const double m33Rgb2Lms[3][3], double m33Rgb2Xyz[3][3]);

void GetIpt2LmsM33(double m33Ipt2Lms[3][3]);
void GetIpt2LmsM33Narrow(int bits, double m33Ipt2Lms[3][3]);
void GetIctcpIo2LmsM33(double m33Ipt2Lms[3][3]);
void GetIctcpIo2LmsM33Narrow(int bits, double m33IctcpIo2Lms[3][3]);

void GetLms2IctcpDmM33(double m33Ipt2Lms[3][3]);

void GetLms2XyzM33(double m33Lms2Xyz[3][3]);

void GetWpV3(CWpDef_t wpDef, double v3Wp[3]);

// PQ is denormalized
void GaussianFilter(int MSRadius, double MsSigma, double *gaussFltr);

#define EOTF_C2K(c) (((c) == CEotfBt1886) ? KEotfBt1886 : ((c) == CEotfPq) ? KEotfPq \
                                                                           : KEotfPower)

int InitDEITPCore(DmKsFlt_t *hKs);

//// C cant override function so

void GetEotfParamsFc(double lMin, double lMax, double gamma, float *pA, float *pB);

void GetYuvRgbOffFc(CRng_t rng, int bits, float v3YuvRgbOff[3]);
void GetYuv2RgbM33Fc(CYuvXferSpec_t yuvXferSpec, float m33Yuv2Rgb[3][3]);
void GetYuv2RgbM33NarrowFc(CYuvXferSpec_t yuvXferSpec, int bits, float m33Yuv2Rgb[3][3]);

void GetRgb2XyzFrom2LmsM33Fc(const double m33Rgb2Lms[3][3], float m33Rgb2Xyz[3][3]);

void GetIpt2LmsM33Fc(float m33Ipt2Lms[3][3]);
void GetIpt2LmsM33NarrowFc(int bits, float m33Ipt2Lms[3][3]);
void GetIctcpIo2LmsM33Fc(float m33Ipt2Lms[3][3]);
void GetIctcpIo2LmsM33NarrowFc(int bits, float m33IctcpIo2Lms[3][3]);

void GetLms2IctcpDmM33Fc(float m33Lms2IctcpDm[3][3]);

void GetLms2XyzM33Fc(float m33Lms2Xyz[3][3]);

void GetWpV3Fc(CWpDef_t wpDef, float v3Wp[3]);

void GaussianFilterFc(int MSRadius, double MsSigma, float *gaussFltr);

#if defined(c_plusplus) || defined(__cplusplus)
}
#endif

#endif // C_DE_UTIL_FLT_H
