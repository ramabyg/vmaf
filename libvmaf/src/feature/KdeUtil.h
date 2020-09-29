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


#ifndef K_DE_UTIL_FLT_H
#define K_DE_UTIL_FLT_H

#include "KdeType.h"
#include <math.h>  // for pow(,) fabs() only

# define SAVE_2_DBG_BUF444(r_, c_, u_, v_, w_)

# define INIT_SAVE_2_DBG_BUF422()
# define SAVE_2_DBG_BUF422(u_, mn, mx)


//////////////////////// math ////////////////////////
#define POW_FUNC powf
#define ABS_FUNC fabsf

//////////////////////// utilities ////////////////////////
#define MIN(a, b)  (((a) <= (b)) ? (a) : (b))

#define MINV3(x0, x1, x2, mx)\
  x0 = MIN(x0, mx);          \
  x1 = MIN(x1, mx);          \
  x2 = MIN(x2, mx)

#define MAX(a, b)  (((a) >= (b)) ? (a) : (b))

#define MAXV3(x0, x1, x2, mn)\
  x0 = MAX(x0, mn);          \
  x1 = MAX(x1, mn);          \
  x2 = MAX(x2, mn)

#define CLAMP(a, mn, mx)  ( ((a) >= (mn)) ? ( ((a) <= (mx)) ? (a) : (mx) ) : (mn) )

#define CLAMPV3(x0, x1, x2, mn, mx)\
  x0 = CLAMP(x0, mn, mx);          \
  x1 = CLAMP(x1, mn, mx);          \
  x2 = CLAMP(x2, mn, mx)

#define NV3(slopeR, offset, x0, x1, x2) \
  x0 = (x0 - (offset))*(slopeR);        \
  x1 = (x1 - (offset))*(slopeR);        \
  x2 = (x2 - (offset))*(slopeR)

#define DNV3(slope, offset, x0, x1, x2)  \
  x0 = (slope)*x0 + (offset);            \
  x1 = (slope)*x1 + (offset);            \
  x2 = (slope)*x2 + (offset)

#define V3FIRST_EQ_V3SECOND(x0, x1, x2, y0,  y1, y2)  \
  x0 = y0;                                   \
  x1 = y1;                                   \
  x2 = y2


#define V3_ME_V3(x0, x1, x2, y0,  y1, y2)  \
  x0 -= y0;                                \
  x1 -= y1;                                \
  x2 -= y2

#define V3_M_V3(x0, x1, x2, y0,  y1, y2, z0, z1, z2) \
  z0 = (x0) - (y0);                                  \
  z1 = (x1) - (y1);                                  \
  z2 = (x2) - (y2)

#define V3_PE_V3(x0, x1, x2, y0,  y1, y2)  \
  x0 += y0;                                \
  x1 += y1;                                \
  x2 += y2

#define V3_P_V3(x0, x1, x2, y0, y1, y2, z0, z1, z2)  \
  z0 = (x0) + (y0);                                  \
  z1 = (x1) + (y1);                                  \
  z2 = (x2) + (y2)

#define V3_TE_V3(x0, x1, x2, y0,  y1, y2)  \
  x0 *= y0;                                \
  x1 *= y1;                                \
  x2 *= y2

#define V3_T_V3(x0, x1, x2, y0, y1, y2, z0, z1, z2)  \
  z0 = (x0) * (y0);                                  \
  z1 = (x1) * (y1);                                  \
  z2 = (x2) * (y2)

#define S_T_V3(s, x0, x1, x2) \
  x0 *= s;                    \
  x1 *= s;                    \
  x2 *= s

#define M33_T_V3(m33, x0, x1, x2, y0, y1, y2)                       \
  y0 = (m33)[0][0] * (x0) + (m33)[0][1] * (x1) + (m33)[0][2] * (x2);\
  y1 = (m33)[1][0] * (x0) + (m33)[1][1] * (x1) + (m33)[1][2] * (x2);\
  y2 = (m33)[2][0] * (x0) + (m33)[2][1] * (x1) + (m33)[2][2] * (x2)


////// L <=> gamma
// to linear ecoded value
#define G2L(eotfParam, x)  \
  x += (eotfParam).b;      \
  x = MAX(x, 0);           \
  x = (eotfParam).a * POW_FUNC(x, (eotfParam).gamma)

#define G2LV3(eotfParam, x0, x1, x2)  \
  G2L(eotfParam, x0);                 \
  G2L(eotfParam, x1);                 \
  G2L(eotfParam, x2)

// to gamma ecoded value
#define L2G(oetfParam, x) \
  x = MAX(x, 0);          \
  x = POW_FUNC(x*(oetfParam).aR, (oetfParam).gammaR) - (oetfParam).b

#define L2GV3(oetfParam, x0, x1, x2)  \
  L2G(oetfParam, x0);                 \
  L2G(oetfParam, x1);                 \
  L2G(oetfParam, x2)

////// L <=> power
// to linear encoded value
#define P2L(eotfParam, x)  \
  x = MAX(x, 0);           \
  x = (eotfParam).a * POW_FUNC(x, (eotfParam).gamma) + (eotfParam).b

#define P2LV3(eotfParam, x0, x1, x2)  \
  P2L(eotfParam, x0);                 \
  P2L(eotfParam, x1);                 \
  P2L(eotfParam, x2)

// to power encoded value
#define L2P(oetfParam, x) \
  x -= (oetfParam).b;     \
  x = MAX(x, 0);          \
  x = POW_FUNC(x*(oetfParam).aR, (oetfParam).gammaR);

#define L2PV3(oetfParam, x0, x1, x2)  \
  L2P(oetfParam, x0);                 \
  L2P(oetfParam, x1);                 \
  L2P(oetfParam, x2)

////// L <=> PQ
static const FloatComp_t k_l2pq_n  = (FloatComp_t)(2610/16384.0); // 16384 =4096*4
static const FloatComp_t k_l2pq_m  = (FloatComp_t)(2523/32.0); // 32 = 4096/128

static const FloatComp_t k_pq2l_n  = (FloatComp_t)(16384/2610.0); // 16384 =4096*4
static const FloatComp_t k_pq2l_m  = (FloatComp_t)(32/2523.0); // 32 = 4096/128

static const FloatComp_t k_lpq_c1  = (FloatComp_t)(3424/4096.0);
static const FloatComp_t k_lpq_c2  = (FloatComp_t)(2413/128.0); // 128 = 4096/32
static const FloatComp_t k_lpq_c3  = (FloatComp_t)(2392/128.0);


#define PQN2L(x)                              \
  x = MAX(x, 0);                              \
  x = POW_FUNC(x, k_pq2l_m);                  \
  x = (x - k_lpq_c1)/(k_lpq_c2 - k_lpq_c3*x); \
  x = MAX(x, 0);                              \
  x = POW_FUNC(x, k_pq2l_n)*10000

#define PQN2LV3(x0, x1, x2) \
  PQN2L(x0);                \
  PQN2L(x1);                \
  PQN2L(x2)


#define L2PQN(x)                                \
  x /= 10000;                                   \
  x = MAX(x, 0);                                \
  x = POW_FUNC(x, k_l2pq_n);                    \
  x = (k_lpq_c1 + k_lpq_c2*x)/(1 + k_lpq_c3*x); \
  x = POW_FUNC(x, k_l2pq_m)


#define L2PQNV3(x0, x1, x2) \
  L2PQN(x0);                \
  L2PQN(x1);                \
  L2PQN(x2)

////// 1 D Lut look up
#define LUT1D(lut, lutNodesM1, x, y)                    \
  if ((x) <= 0) y = (lut)[0];                           \
  else if ((x) < 1) {                                   \
    int iX;                                             \
    y = (x)*(lutNodesM1);                               \
    iX = (int)y;                                        \
    y = (lut)[iX] + ((lut)[iX+1] - (lut)[iX])*(y - iX); \
  }                                                     \
  else y = (lut)[(lutNodesM1)]

#define LUT1DV3(lut, lutNodesM1, x0, x1, x2)  \
  LUT1D(lut, lutNodesM1, x0, x0);             \
  LUT1D(lut, lutNodesM1, x1, x1);             \
  LUT1D(lut, lutNodesM1, x1, x2)

#define ALPHA_BLENDINGV3(x0, x1, x2, z0, z1, z2, alpha) \
  (x0) = (x0) + (alpha)*((z0) - (x0));                  \
  (x1) = (x1) + (alpha)*((z1) - (x1));                  \
  (x2) = (x2) + (alpha)*((z2) - (x2))


////// Get/Set pixel value for data type F32, U16, U8 and frame arragement of planar or interleaved
#define GET_OFFSET(frmFmt_, row_, col_)   ((row_)*(frmFmt_).rowPitch + (col_)*(frmFmt_).colPitch)

#define GET_444_TYPE_VALUE(dtp_t, byteOff, frmBuf0, frmBuf1, frmBuf2, x0, x1, x2) \
  x0 =  *(dtp_t *)((char *)frmBuf0 + byteOff);                                    \
  x1 =  *(dtp_t *)((char *)frmBuf1 + byteOff);                                    \
  x2 =  *(dtp_t *)((char *)frmBuf2 + byteOff)


#define GET_444_VALUE(dtp, byteOff, frmBuf0, frmBuf1, frmBuf2, x0, x1, x2)            \
  if (dtp == KDtpF32) {                                                               \
    GET_444_TYPE_VALUE(float,         byteOff, frmBuf0, frmBuf1, frmBuf2, x0, x1, x2);\
  }                                                                                   \
  else if (dtp == KDtpU8) {                                                           \
    GET_444_TYPE_VALUE(unsigned char, byteOff, frmBuf0, frmBuf1, frmBuf2, x0, x1, x2);\
  }                                                                                   \
  else {                                                                              \
    GET_444_TYPE_VALUE(unsigned short,byteOff, frmBuf0, frmBuf1, frmBuf2, x0, x1, x2);\
  }

#define SET_444_TYPE_VALUE(dtp_t,  byteOff, x0, x1, x2, frmBuf0, frmBuf1, frmBuf2)  \
  *(dtp_t *)((char *)frmBuf0 + byteOff) = (dtp_t)(x0);                              \
  *(dtp_t *)((char *)frmBuf1 + byteOff) = (dtp_t)(x1);                              \
  *(dtp_t *)((char *)frmBuf2 + byteOff) = (dtp_t)(x2)

#define DT5    ((FloatComp_t)0.5)

#define SET_444_VALUE(dtp, byteOff, x0, x1, x2, frmBuf0, frmBuf1, frmBuf2)                              \
  if (dtp == KDtpF32) {                                                                                 \
    SET_444_TYPE_VALUE(float,         byteOff, x0,       x1,        x2,      frmBuf0, frmBuf1, frmBuf2);\
  }                                                                                                     \
  else if (dtp == KDtpU8) {                                                                             \
    SET_444_TYPE_VALUE(unsigned char, byteOff, (x0)+DT5, (x1)+DT5, (x2)+DT5, frmBuf0, frmBuf1, frmBuf2);\
  }                                                                                                     \
  else {                                                                                                \
    SET_444_TYPE_VALUE(unsigned short,byteOff, (x0)+DT5, (x1)+DT5, (x2)+DT5, frmBuf0, frmBuf1, frmBuf2);\
  }

#define SET_422_TYPE_FLT_VALUE_INC_ADDR(dtp_t, x, mn, mx, pUyVy)                              \
    *(dtp_t *)(pUyVy) = ((x) < (mn)) ? (dtp_t)(mn) : ((x) < (mx)) ? (dtp_t)(x) : (dtp_t)(mx); \
    pUyVy = (unsigned char *)pUyVy + sizeof(dtp_t)

#define SET_422_TYPE_FXP_VALUE_INC_ADDR(dtp_t, x, mn, mx, pUyVy)                                    \
    *(dtp_t *)(pUyVy) = ((x) < (mn)) ? (dtp_t)(mn) : ((x) < (mx)) ? (dtp_t)((x)+DT5) : (dtp_t)(mx); \
    pUyVy = (unsigned char *)pUyVy + sizeof(dtp_t)


#define SET_422_VALUE_INC_ADDR(dtp, x, mn, mx, pUyVy)                  \
  SAVE_2_DBG_BUF422(x, mn, mx);                                        \
  if (dtp == KDtpF32) {                                                \
    SET_422_TYPE_FLT_VALUE_INC_ADDR(float, x, mn, mx, pUyVy);          \
  }                                                                    \
  else if (dtp == KDtpU8) {                                            \
    SET_422_TYPE_FXP_VALUE_INC_ADDR(unsigned char, x, mn, mx, pUyVy);  \
  }                                                                    \
  else {                                                               \
    SET_422_TYPE_FXP_VALUE_INC_ADDR(unsigned short, x, mn, mx, pUyVy); \
  }

#define SET_VALUE(dtp, x, mn, mx, pV)                                  \
  if (dtp == KDtpF32) {                                                \
    *(float *)(pV)          = ((x) < (mn)) ? (float)(mn) : ((x) < (mx)) ? (float)(x) : (float)(mx); \
  }                                                                    \
  else if (dtp == KDtpU8) {                                            \
    *(unsigned char *)(pV)  = ((x) < (mn)) ? (unsigned char)(mn)  : ((x) < (mx)) ? (unsigned char)((x)+DT5)  : (unsigned char)(mx); \
  }                                                                    \
  else {                                                               \
    *(unsigned short *)(pV) = ((x) < (mn)) ? (unsigned short)(mn) : ((x) < (mx)) ? (unsigned short)((x)+DT5) : (unsigned short)(mx);\
  }


#endif // K_DE_UTIL_FLT_H
