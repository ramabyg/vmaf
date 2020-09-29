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
#include <math.h> // for pow, exp, only
#include <assert.h>
#include "KdeType.h"
#include "CdeType.h"
#include "CdeUtil.h"

#ifdef WIN32
#pragma warning( push )
#pragma warning( disable : 4995 4996 )
#endif

                          static void
                          AssignM33D2D(const double m33In[3][3], double m33Out[3][3])
{
  int i, j;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      m33Out[i][j] = m33In[i][j];
}

static void AssignV3D2D(const double v3In[3], double v3Out[3])
{
  v3Out[0] = v3In[0];
  v3Out[1] = v3In[1];
  v3Out[2] = v3In[2];
}

// to avoid to many functions, all op in double
// use the assignment to convert double to single
static void AssignM33D2S(const double m33In[3][3], float m33Out[3][3])
{
  int i, j;
  for (i = 0; i < 3; ++i)
    for (j = 0; j < 3; ++j)
      m33Out[i][j] = (float)m33In[i][j];
}

void AssignV3D2S(const double v3In[3], float v3Out[3])
{
  v3Out[0] = (float)v3In[0];
  v3Out[1] = (float)v3In[1];
  v3Out[2] = (float)v3In[2];
}

void AssignM33D2(const double m33In[3][3], FloatComp_t m33Out[3][3])
{
  AssignM33D2S(m33In, m33Out);
}

void AssignV3D2(const double v3In[3], FloatComp_t v3Out[3])
{
  AssignV3D2S(v3In, v3Out);
}

static void M33TimesM33(const double m33a[3][3], const double m33b[3][3], double m33[3][3])
{
  int i, j, k;
  double acc;
  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      acc = 0;
      for (k = 0; k < 3; ++k) {
        acc += m33a[i][k]*m33b[k][j];
      }
      m33[i][j] = acc;
    }
  }
}

// inv_mat_3x3(): inverts 3x3 matrix
static void InvM33(const double m33In[3][3], double m33Out[3][3])
{
  int i,j;
  double det=0;

  for(i=0; i<3; i++)
    det += m33In[0][i]*(m33In[1][(i+1)%3]*m33In[2][(i+2)%3] -
                        m33In[1][(i+2)%3]*m33In[2][(i+1)%3]);
  for(i=0; i<3; i++)
    for(j=0; j<3; j++)
      m33Out[j][i] =
      (m33In[(i+1)%3][(j+1)%3] * m33In[(i+2)%3][(j+2)%3] -
       m33In[(i+1)%3][(j+2)%3] * m33In[(i+2)%3][(j+1)%3]) / det;
}

/* to get a single precision m33 output from double m33 output function verion */
#define GETM33FC_FROM_DB0(funcDb, out)          \
{                                               \
  double m33Db[3][3];                           \
  funcDb(m33Db);                                \
  AssignM33D2S((const double (*)[3])m33Db, out);\
}

#define GETM33FC_FROM_DB1(funcDb, in1, out)     \
{                                               \
  double m33Db[3][3];                           \
  funcDb(in1, m33Db);                           \
  AssignM33D2S((const double (*)[3])m33Db, out);\
}

#define GETM33FC_FROM_DB2(funcDb, in1, in2, out)  \
{                                                 \
  double m33Db[3][3];                             \
  funcDb(in1, in2, m33Db);                        \
  AssignM33D2S((const double (*)[3])m33Db, out);  \
}

#define GETM33FC_FROM_DB3(funcDb, in1, in2, in3, out)\
{                                                    \
  double m33Db[3][3];                                \
  funcDb(in1, in2, in3, m33Db);                      \
  AssignM33D2S((const double (*)[3])m33Db, out);     \
}

/* to get a reverse m33 matrix from get forward m33 function verion */
#define GETM33_REV1(funcFw, in1, out)     \
{                                         \
  double m33Fw[3][3];                     \
  funcFw(in1, m33Fw);                     \
  InvM33((const double (*)[3])m33Fw, out);\
}

#define GETM33_REV2(funcFw, in1, in2, out)  \
{                                           \
  double m33Fw[3][3];                       \
  funcFw(in1, in2, m33Fw);                  \
  InvM33((const double (*)[3])m33Fw, out);  \
}

#define GETM33_REV3(funcFw, in1, in2, in3, out)  \
{                                                \
  double m33Fw[3][3];                            \
  funcFw(in1, in2, in3, m33Fw);                  \
  InvM33((const double (*)[3])m33Fw, out);       \
}

#define GETM33FC_FROM_DB8(funcDb, in1, in2, in3, in4, in5, in6, in7, in8, out)\
{                                                                             \
  double m33Db[3][3];                                                         \
  funcDb(in1, in2, in3, in4, in5, in6, in7, in8,m33Db);                       \
  AssignM33D2S((const double (*)[3])m33Db, out);                              \
}


void MemCpyByte(void *pDst, const void *pSrc, int len)
{
  char *pD = (char *)pDst;
  const char *pS = (char *)pSrc;
  const char *pSE = pS + len;


  while (pS != pSE) *pD++ = *pS++;
}

int MemEqByte(const void *pD1, const void *pD2, int len)
{
  const char *pB1 = (char *)pD1;
  const char *pB2 = (char *)pD2;
  const char *pB1E = pB1 + len;


  while (pB1 != pB1E && *pB1 == *pB2) {
    ++pB1;
    ++pB2;
  }

  return pB1 == pB1E;
}

static double GetNarrow2FullYRatio(int bits)
{
  return (double)((1 << bits) - 1) / ((235 - 16)*(1 << (bits - 8)));
}

static double GetNarrow2FullUvRatio(int bits)
{
  return (double)((1 << bits) - 1) / ((240 - 16)*(1 << (bits - 8)));
}


////// L <=> PQ
static const double l2pq_c1  = 3424/4096.0;
static const double l2pq_c2  = 2413/128.0;
static const double l2pq_c3  = 2392/128.0;

// L2PQ(): converts luminance to perceptual quantized
// - luminance (0:10000 double precision nits),
// perceptual quantized (0:4095 double precision "code words")
// - L is in nits, PQ is in 12bit code values
// - Conformed parameters to ITU submission, July 2012
// here the PQ is normalized and assuming input is in the well defined range
double L2PqNorm(double L)
{
  static const double l2pq_n  = 2610/16384.0;
  static const double l2pq_m  = 2523/32.0;

  // Calculate PQ code values (no rounding)
  double Yn = pow(L/10000.0, l2pq_n);
  double V = pow((l2pq_c1 + l2pq_c2 * Yn)/(1.0 + l2pq_c3 * Yn), l2pq_m);

  return CLAMPS(V, 0.0, 1.0);
}

// PQ2L(): converts perceptual quantized to luminance
// - perceptual quantized (0:4095 double precision "code words"),
// luminance (0:10000 double precision nits)
// - PQ is in 12bit code values, L is in nits
// - Conformed parameters to ITU submission, July 2012
// here the PQ is normalized and assuming input is in the well defined range
double PqNorm2L(double PQ)
{
  static const double l2pq_n_r  = 16384/2610.0; // 1/l2pq_n;
  static const double l2pq_m_r  = 32/2523.0; // 1/l2pq_m;

  double Vm = pow(PQ, l2pq_m_r);
  double L  = pow( MAX2S(Vm - l2pq_c1, 0.0) / (l2pq_c2 - l2pq_c3 * Vm), l2pq_n_r) * 10000.0;

  return CLAMPS(L, 0.0, 10000.0);;
}


// E O tf
void GetEotfParams(double lMin, double lMax, double gamma,  double *pA, double *pB)
{
  double maxG, minG, rangeG;

  minG = pow(lMin, 1.0/gamma);
  maxG = pow(lMax, 1.0/gamma);
  rangeG = maxG - minG;
  *pA = pow(rangeG, gamma);
  *pB = minG / rangeG;
}
void GetEotfParamsFc(double lMin, double lMax, double gamma,  float *pA, float *pB)
{
  double a, b;

  GetEotfParams(lMin, lMax, gamma, &a, &b);
  *pA = (float)a;
  *pB = (float)b;
}

// Xyy2Xyz(): calculates XYZ tristimulus from xy chromaticity and Y luminance
//
// Xyy2Xyz() can be called in three different ways:
// 1. To calculate XYZ from xyY
// 2. To calculate XYZ from a media white where Y=100
// 3. To calculate XYZ from xy chromaticity values and the corresponding media white xy or XYZ
//
// In case 1, xyy is an nx3 array where the third column is Y luminance
// Example: Xyy2Xyz (xyY, XYZ)
//
// In case 2, xyy is an nx2 array where Y is assumed to be 100
// Example: Xyy2Xyz (xy, XYZ)
//
// In case 3, xyy is a 3x2 array where the rows correspond to the xy values
// of the R, G, and B primaries; and xyzw is the media white. The media white
// can be given as xy or XYZ values
// Example: Xyy2Xyz (xy, xyzw, XYZ)
//
// For our application, only Case 1 is implemented (Y is assumed to be 1.0)
static void Xyy2Xyz(double x, double y, double v3Xyz[3])
{
  v3Xyz[1] = 1;            // Y is 1.0
  v3Xyz[0] = (x/y);        // (x/y)*Y
  v3Xyz[2] = ((1.0-x-y)/y);  // ((1-x-y)/y)*Y
}

static void GetRgb2XyzM33FromPrims(double rx, double ry, double gx, double gy, double bx, double by,
              double wx, double wy, double m33Rgb2Xyz[3][3])
{
  //// caulculate u, v, w
  double D =  (rx - bx)*(gy - by) - (ry - by)*(gx - bx);
  double u = ((wx - bx)*(gy - by) - (wy - by)*(gx - bx))/D;
  double v = ((rx - bx)*(wy - by) - (ry - by)*(wx - bx))/D;
  double w = 1.0 - u - v;

  //// scale u, v, w by wy
  u /= wy;
  v /= wy;
  w /= wy;

  //// get m33Rgb2Xyz
  m33Rgb2Xyz[0][0] = (rx*u);
  m33Rgb2Xyz[0][1] = (gx*v);
  m33Rgb2Xyz[0][2] = (bx*w);

  m33Rgb2Xyz[1][0] = (ry*u);
  m33Rgb2Xyz[1][1] = (gy*v);
  m33Rgb2Xyz[1][2] = (by*w);

  m33Rgb2Xyz[2][0] = ((1.0 - rx - ry)*u); // rz*u
  m33Rgb2Xyz[2][1] = ((1.0 - gx - gy)*v); // gz*u
  m33Rgb2Xyz[2][2] = ((1.0 - bx - by)*w); // bz*u
}


// getRgb2XyzM33(): returns the matrix required to go from the specified RGB space to XYZ
//  - The target space can be one of,
//    'p3d65', 'dci', 'r709', 'r2020',  'aces'
static void GetRgb2XyzM33(CRgbDef_t rgbDef, double m33Rgb2Xyz[3][3])
{
  // temp var use double for accuracy
  double rx,ry,gx,gy,bx,by,wx,wy;

  if (rgbDef == CRgbDefP3d65) {
    rx = 0.6800; ry = 0.3200;
    gx = 0.2650; gy = 0.6900;
    bx = 0.1500; by = 0.0600;
    wx = 0.3127; wy = 0.3290;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else if (rgbDef == CRgbDefDci) {
    rx = 0.6800; ry = 0.3200;
    gx = 0.2650; gy = 0.6900;
    bx = 0.1500; by = 0.0600;
    wx = 0.3140; wy = 0.3510;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else if (rgbDef == CRgbDefR709) {
    rx = 0.6400; ry = 0.3300;
    gx = 0.3000; gy = 0.6000;
    bx = 0.1500; by = 0.0600;
    wx = 0.3127; wy = 0.3290;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else if (rgbDef == CRgbDefR2020) {
    rx = 0.7080; ry = 0.2920;
    gx = 0.1700; gy = 0.7970;
    bx = 0.1310; by = 0.0460;
    wx = 0.3127; wy = 0.3290;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else if (rgbDef == CRgbDefAces) {
    rx = 0.73470; ry = 0.26530;
    gx = 0.00000; gy = 1.00000;
    bx = 0.00010; by = -0.07700;
    wx = 0.32168; wy = 0.33767;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else if (rgbDef == CRgbDefAlexa) {
    rx = 0.6840; ry = 0.3130;
    gx = 0.2210; gy = 0.8480;
    bx = 0.0861; by = -0.1020;
    wx = 0.3127; wy = 0.3290;
    GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  }
  else {
    assert(rgbDef < CRgbDefNum);
  }
}

static void GetRgb2YuvFromRgb2XyzM33(double m33Rgb2Xyz[3][3], double m33Rgb2Yuv[3][3])
{
  m33Rgb2Yuv[0][0] = m33Rgb2Xyz[1][0];
  m33Rgb2Yuv[0][1] = m33Rgb2Xyz[1][1];
  m33Rgb2Yuv[0][2] = m33Rgb2Xyz[1][2];

  m33Rgb2Yuv[1][0] = -(0.5*m33Rgb2Xyz[1][0])/(1 - m33Rgb2Xyz[1][2]);
  m33Rgb2Yuv[1][1] = -(0.5*m33Rgb2Xyz[1][1])/(1 - m33Rgb2Xyz[1][2]);
  m33Rgb2Yuv[1][2] = 0.5;

  m33Rgb2Yuv[2][0] = 0.5;
  m33Rgb2Yuv[2][1] = -(0.5*m33Rgb2Xyz[1][1])/(1 - m33Rgb2Xyz[1][0]);
  m33Rgb2Yuv[2][2] = -(0.5*m33Rgb2Xyz[1][2])/(1 - m33Rgb2Xyz[1][0]);
}

void GetRgb2YuvM33(CYuvXferSpec_t yuvXferSpec, double m33Rgb2Yuv[3][3])
{
  if (yuvXferSpec != CYuvXferSpecR601) {
    double m33Rgb2Xyz[3][3] ={{1, 0, 0}, {0, 1, 0}, {0, 0, 1}}; // to make compiler happy

    GetRgb2XyzM33((CRgbDef_t)yuvXferSpec, m33Rgb2Xyz);
    GetRgb2YuvFromRgb2XyzM33(m33Rgb2Xyz, m33Rgb2Yuv);
  }
  else {
    m33Rgb2Yuv[0][0] =  0.299000000000000;
    m33Rgb2Yuv[0][1] =  0.587000000000000;
    m33Rgb2Yuv[0][2] =  0.114000000000000;

    m33Rgb2Yuv[1][0] = -0.168735891647856;
    m33Rgb2Yuv[1][1] = -0.331264108352145;
    m33Rgb2Yuv[1][2] =  0.500000000000000;

    m33Rgb2Yuv[2][0] =  0.500000000000000;
    m33Rgb2Yuv[2][1] = -0.418687589158345;
    m33Rgb2Yuv[2][2] = -0.0813124108416548;
  }
}
void GetRgb2YuvM33Fc(CYuvXferSpec_t yuvXferSpec, float m33Rgb2Yuv[3][3])
GETM33FC_FROM_DB1(GetRgb2YuvM33, yuvXferSpec, m33Rgb2Yuv)

/* narrow version */
static void Rgb2YuvM33Full2Narrow(double m33Full[3][3], int bits, double m33Narrow[3][3])
{
  // from full to narrow
  double yRation =  ( (235-16)*(1<<(bits - 8)) )/(double)((1<<bits) - 1);
  double uvRation = ( (240-16)*(1<<(bits - 8)) )/(double)((1<<bits) - 1);
  int c;

  for (c = 0; c < 3; ++c) {
    m33Narrow[0][c] = m33Full[0][c]*yRation;
    m33Narrow[1][c] = m33Full[1][c]*uvRation;
    m33Narrow[2][c] = m33Full[2][c]*uvRation;
  }
}

void GetRgb2YuvM33Narrow(CYuvXferSpec_t yuvXferSpec, int bits, double m33Rgb2Yuv[3][3])
{
  double m33Rgb2YuvFull[3][3];

  GetRgb2YuvM33(yuvXferSpec, m33Rgb2YuvFull);
  Rgb2YuvM33Full2Narrow(m33Rgb2YuvFull, bits, m33Rgb2Yuv);
}
void GetRgb2YuvM33NarrowFc(CYuvXferSpec_t yuvXferSpec, int bits, float m33Rgb2Yuv[3][3])
GETM33FC_FROM_DB2(GetRgb2YuvM33Narrow, yuvXferSpec, bits, m33Rgb2Yuv)


void GetYuvRgbOff(CRng_t rng, int bits, double v3YuvRgbOff[3])
{
  v3YuvRgbOff[0] = ((rng == CRngNarrow) ? (1<<(bits - 4)) : 0);
  v3YuvRgbOff[1] = v3YuvRgbOff[2] = (1<<(bits - 1));
}
void GetYuvRgbOffFc(CRng_t rng, int bits, float v3YuvRgbOff[3])
{
  double v3D[3];

  GetYuvRgbOff(rng, bits, v3D);
  AssignV3D2S(v3D, v3YuvRgbOff);
}

static void Yuv2RgbM33Narrow2Full(const double m33Full[3][3], int bits, double m33Narrow[3][3])
{
  // from narrow to full
  double yRation = GetNarrow2FullYRatio(bits);
  double uvRation = GetNarrow2FullUvRatio(bits);
  int r;

  for (r = 0; r < 3; ++r) {
    m33Narrow[r][0] = m33Full[r][0] * yRation;
    m33Narrow[r][1] = m33Full[r][1] * uvRation;
    m33Narrow[r][2] = m33Full[r][2] * uvRation;
  }
}

static const double m33Xyz2LmsFlt[3][3] =
  {  { 0.400238220153002,   0.707593156246992, -0.080805581487871},
    {-0.226298103801373,   1.165315591308231,  0.045700774541672},
    { 0.0,                 0.0,                0.918224951158247} };

static const double m33Lms2XyzFlt[3][3] =
{    {1.859949326962933, -1.129382825168997,  0.219889425257736},
    {0.361192289022606,  0.638816483492423, -0.000008772515029},
    {0.0,                0.0,                1.089057750759878}  };

/* These values must be passed with external inputs, 
   for now set to Zero as external inputs are not supported
*/
static const double V3Yuv2RgbOffNorm[3] =
    {0.0, 0.0, 0.0};
static const double V3Yuv2RgbOff[3] =
    {0.0, 0.0, 0.0};
static const double M33Rgb2Lms[3][3] = {{0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0}};
static const double V8Rgbw[8] = {
    0.0,0.0,0.0,0.0,
    0.0,0.0,0.0,0.0,
};

void
GetRgb2LmsByDefM33(CRgbDef_t rgbDef, double m33Rgb2Lms[3][3])
{
  double m33Rgb2Xyz[3][3];

  GetRgb2XyzM33(rgbDef, m33Rgb2Xyz);
  M33TimesM33((const double (*)[3])m33Xyz2LmsFlt, (const double (*)[3])m33Rgb2Xyz, m33Rgb2Lms);
}
void GetRgb2LmsByDefM33Fc(CRgbDef_t rgbDef, float m33Rgb2Lms[3][3])
GETM33FC_FROM_DB1(GetRgb2LmsByDefM33, rgbDef, m33Rgb2Lms)

void GetRgb2LmsByPrimsM33(double rx, double ry, double gx, double gy, double bx, double by,
                          double wx, double wy, double m33Rgb2Lms[3][3])
{
  double m33Rgb2Xyz[3][3];

  GetRgb2XyzM33FromPrims(rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Xyz);
  M33TimesM33((const double (*)[3])m33Xyz2LmsFlt, (const double (*)[3])m33Rgb2Xyz, m33Rgb2Lms);
}
// void GetRgb2LmsByPrimsM33Fc(CRgbDef_t rgbDef, float m33Rgb2Lms[3][3])
// GETM33FC_FROM_DB8(GetRgb2LmsByPrimsM33, rx,ry,gx,gy,bx,by,wx,wy, m33Rgb2Lms)

void GetRgb2LmsCtWpM33(const double m33[3][3], double crossTalk, const double v3Wp[3], double m33Rgb2Lms[3][3])
{
  double colFct[3], d1m3ct = 1 - 3*crossTalk;
  int i, j;

  for (j = 0; j < 3; ++j) {
    colFct[j] = crossTalk*(m33[0][j] + m33[1][j] + m33[2][j]);
  }

  for (i = 0; i < 3; ++i) {
    for (j = 0; j < 3; ++j) {
      m33Rgb2Lms[i][j] = (d1m3ct*m33[i][j] + colFct[j])/v3Wp[i];
    }
  }
}
void GetRgb2LmsCtWpM33Fc(const double m33[3][3], double crossTalk, const double v3Wp[3], float m33Rgb2Lms[3][3])
GETM33FC_FROM_DB3(GetRgb2LmsCtWpM33, m33, crossTalk, v3Wp, m33Rgb2Lms)

/* RGB<=>LMS with white balance and cross talk*/
void GetRgb2XyzFrom2LmsM33(const double m33Rgb2Lms[3][3], double m33Rgb2Xyz[3][3])
{
  M33TimesM33(m33Lms2XyzFlt, m33Rgb2Lms, m33Rgb2Xyz);
}
void GetRgb2XyzFrom2LmsM33Fc(const double m33Rgb2Lms[3][3], float m33Rgb2Xyz[3][3])
GETM33FC_FROM_DB1(GetRgb2XyzFrom2LmsM33, m33Rgb2Lms, m33Rgb2Xyz)

void GetYuv2RgbM33(CYuvXferSpec_t yuvXferSpec, double m33Yuv2Rgb[3][3])
GETM33_REV1(GetRgb2YuvM33, yuvXferSpec, m33Yuv2Rgb)
void GetYuv2RgbM33Fc(CYuvXferSpec_t yuvXferSpec, float m33Yuv2Rgb[3][3])
GETM33FC_FROM_DB1(GetYuv2RgbM33, yuvXferSpec, m33Yuv2Rgb)

static void Yuv2RgbM33Full2Narrow(double m33Full[3][3], int bits, double m33Narrow[3][3])
{
  // from full to narrow
  double yRation =  (double)((1<<bits) - 1)/( (235-16)*(1<<(bits - 8)) );
  double uvRation = (double)((1<<bits) - 1)/( (240-16)*(1<<(bits - 8)) );
  int r;

  for (r = 0; r < 3; ++r) {
    m33Narrow[r][0] = m33Full[r][0]*yRation;
    m33Narrow[r][1] = m33Full[r][1]*uvRation;
    m33Narrow[r][2] = m33Full[r][2]*uvRation;
  }
}

void GetYuv2RgbM33Narrow(CYuvXferSpec_t yuvXferSpec, int bits, double m33Yuv2Rgb[3][3])
{
  double m33Yuv2RgbFull[3][3];

  GetYuv2RgbM33(yuvXferSpec, m33Yuv2RgbFull);
  Yuv2RgbM33Full2Narrow(m33Yuv2RgbFull, bits, m33Yuv2Rgb);
}
void GetYuv2RgbM33NarrowFc(CYuvXferSpec_t yuvXferSpec, int bits, float m33Yuv2Rgb[3][3])
GETM33FC_FROM_DB2(GetYuv2RgbM33Narrow, yuvXferSpec, bits, m33Yuv2Rgb)


static const double m33Lms2IptFlt[3][3] =
{    {0.4000,  0.4000,  0.2000},
    {4.4550, -4.8510,  0.3960},
    {0.8056,  0.3572, -1.1628}  };

void GetLms2IptM33(double m33Lms2Ipt[3][3])
{
  AssignM33D2D(m33Lms2IptFlt, m33Lms2Ipt);
}
void GetLms2IptM33Fc(float m33Lms2Ipt[3][3])
{
  AssignM33D2S(m33Lms2IptFlt, m33Lms2Ipt);
}

static const double m33Ipt2LmsFlt[3][3] =
{   {1.0,  0.097568930514614,  0.205226433164592},
    {1.0, -0.113876485473147,  0.133217158369998},
    {1.0,  0.032615109917066, -0.676887183069179}  };

void GetIpt2LmsM33(double m33Ipt2Lms[3][3])
{
  AssignM33D2D(m33Ipt2LmsFlt, m33Ipt2Lms);
}
void GetIpt2LmsM33Narrow(int bits, double m33Ipt2Lms[3][3])
{
  double yRatio = GetNarrow2FullYRatio(bits);

  GetIpt2LmsM33(m33Ipt2Lms);
  S_TIMES_M33(yRatio, m33Ipt2Lms);
}
void GetIpt2LmsM33Fc(float m33Ipt2Lms[3][3])
{
  AssignM33D2S(m33Ipt2LmsFlt, m33Ipt2Lms);
}
void GetIpt2LmsM33NarrowFc(int bits, float m33Ipt2Lms[3][3])
GETM33FC_FROM_DB1(GetIpt2LmsM33Narrow, bits, m33Ipt2Lms)

static const double m33IctcpIo2LmsFlt[3][3] =
{
  { 1,  0.008609037037933,  0.111029625003026 },
  { 1, -0.008609037037933, -0.111029625003026 },
  { 1,  0.560031335710679, -0.320627174987319 }
};
void GetIctcpIo2LmsM33(double m33Ipt2Lms[3][3])
{
  AssignM33D2D(m33IctcpIo2LmsFlt, m33Ipt2Lms);
}
void GetIctcp2LmsM33Fc(float m33Ipt2Lms[3][3])
{
  AssignM33D2S(m33IctcpIo2LmsFlt, m33Ipt2Lms);
}
void GetIctcpIo2LmsM33Narrow(int bits, double m33IctcpIo2Lms[3][3])
{
  /* same rule */
  Yuv2RgbM33Narrow2Full(m33IctcpIo2LmsFlt, bits, m33IctcpIo2Lms);
}
void GetIctcpIo2LmsM33NarrowFc(int bits, float m33IctcpIo2Lms[3][3])
{
  GETM33FC_FROM_DB1(GetIctcpIo2LmsM33Narrow, bits, m33IctcpIo2Lms);
}

void GetIctcpIo2LmsM33Fc(float m33Ipt2Lms[3][3])
{
  GETM33FC_FROM_DB0(GetIctcpIo2LmsM33, m33Ipt2Lms);
}

/*** lms <=> ictcpDm ***/
static const double m33Lms2IctcpDmFlt[3][3] =
{
  { 2048.0 / 4096,     2048.0 / 4096,               0 },
  { 6610.0 / 4096 / 2, -13613.0 / 4096 / 2, 7003.0 / 4096 / 2 }, /* m33Lms2IctcpFlt[1][]/2 */
  { 17933.0 / 4096,   -17390.0 / 4096,   -543.0 / 4096 }
};
void GetLms2IctcpDmM33(double m33Lms2IctcpDm[3][3])
{
  AssignM33D2D(m33Lms2IctcpDmFlt, m33Lms2IctcpDm);
}
void GetLms2IctcpDmM33Fc(float m33Lms2IctcpDm[3][3])
{
  AssignM33D2S(m33Lms2IctcpDmFlt, m33Lms2IctcpDm);
}

void GetLms2XyzM33(double m33Lms2Xyz[3][3])
{
  AssignM33D2D(m33Lms2XyzFlt, m33Lms2Xyz);
}
void GetLms2XyzM33Fc(float m33Lms2Xyz[3][3]) {
  AssignM33D2S(m33Lms2XyzFlt, m33Lms2Xyz);
}

void GetWpV3(CWpDef_t wpDef, double v3Wp[3])
{
  if (wpDef == CWpDefD60) {
    Xyy2Xyz(0.32168, 0.33767, v3Wp);
  }
  else if (wpDef == CWpDefD65) {
    Xyy2Xyz(0.3127, 0.3290, v3Wp);
  }
  else {
    assert(wpDef < CWpDefNum);
  }
}
void GetWpV3Fc(CWpDef_t wpDef, float v3Wp[3])
{
  double v3D[3];

  GetWpV3(wpDef, v3D);
  AssignV3D2S(v3D, v3Wp);
}

void GaussianFilter(int MSRadius, double MsSigma, double *gaussFltr)
{
  const double ss2 = 2*MsSigma*MsSigma;
  double sum = 0.0;
  int i;

  for (i = 1; i <= MSRadius; ++i) {
    gaussFltr[i] = exp( -(i*i)/ss2 );
    sum += gaussFltr[i];
  }
  gaussFltr[0] = 1;
  sum = 2*sum + 1;

  for (i = 0; i <= MSRadius; ++i) {
    gaussFltr[i] /=  sum;
  }
}

void GaussianFilterFc(int MSRadius, double MsSigma, float *gaussFltr)
{
# define dbBufLen  16

  double dbFltr[dbBufLen];
  int i;

  assert(dbBufLen > MSRadius);

  GaussianFilter(MSRadius, MsSigma, dbFltr);
  for (i = 0; i <= MSRadius; ++i) {
    gaussFltr[i] = (float)dbFltr[i];
  }
}

static void GetRange(int inBits, CRng_t rng, CEotf_t eotf, FloatComp_t *rangeMin, FloatComp_t *range)
{
  if (rng == CRngNarrow)
  {
    *rangeMin = (unsigned short)(16 * (1 << (inBits - 8)));
    *range = (unsigned short)(235 * (1 << (inBits - 8)) - (*rangeMin));
  }
  else if (rng == CRngSdi)
  {
    *rangeMin = (unsigned short)(1 << (inBits - 8));
    if (eotf == CEotfPq)
    {
      *range = (unsigned short)(1019 * (1 << (inBits - 10)) - *rangeMin); // 1019*(1<< inBits)/1024
    }
    else
    {
      *range = (unsigned short)((1 << inBits) - 1 - 2 * (*rangeMin));
    }
  }
  else
  {
    *rangeMin = 0;
    *range = (unsigned short)((1 << inBits) - 1);
  }
}

static double DeriveHlgGamma(double maxLd)
{
  double hlgGamma;

  if (maxLd >= 400 && maxLd <= 2000)
    hlgGamma = 1.2 + 0.42 * log10(maxLd / 1000);
  else
    hlgGamma = 1.2 * pow(1.111, log2(maxLd / 1000));

  /*hlgGamma = MAX2S(hlgGamma, 1);*/

  return hlgGamma;
}

static double DeriveHlgAlpha(double minLd, double maxLd, double hlgGamma)
{
  double alpha;

  alpha = (FloatComp_t)maxLd;
  (void)minLd;
  (void)hlgGamma;

  return alpha;
}

static double DeriveHlgBeta(double minLd, double maxLd, double hlgGamma)
{
  double beta;

  if (minLd)
    beta = 1.732050807568877 * pow(minLd / maxLd, 1 / (2 * hlgGamma)); /*sqrt(3*pow(minLd/maxLd, 1/hlgGamma));*/
  else
    beta = 0;

  return beta;
}

static void DeriveHlgEotfParam(double minLd, double maxLd, EotfParamFlt_t *pEotfParam)
{
  double hlgGamma;

  hlgGamma = DeriveHlgGamma(maxLd);
  pEotfParam->gamma = (FloatComp_t)hlgGamma;

  pEotfParam->a = (FloatComp_t)DeriveHlgAlpha(minLd, maxLd, hlgGamma);
  pEotfParam->b = (FloatComp_t)DeriveHlgBeta(minLd, maxLd, hlgGamma);

#if EN_HARDWARE_HLG
  pEotfParam->gamma -= 1; /* avoid kernel -1 */
#endif
}

/*
  Update EOTF related parameters           
  Here it is always considered full range as libvmaf currently supports yuv only
  */

static void DeitpUpdateEotfParam(DmKsFlt_t *hKs, CEotf_t eotf, FloatComp_t gamma, int bpc)
{

  CRng_t RngLocal = (IS_CCLR_RGB_RGBA(hKs->clr_fmt)) ? hKs->Rng : CRngFull;

  GetRange(bpc, RngLocal, eotf,
           &hKs->ksIMap[0].eotfParam.rangeMin, &hKs->ksIMap[0].eotfParam.rangeR);
  GetRange(bpc, RngLocal, eotf,
           &hKs->ksIMap[1].eotfParam.rangeMin, &hKs->ksIMap[1].eotfParam.rangeR);

  hKs->ksIMap[0].eotfParam.rangeR = 1 / hKs->ksIMap[0].eotfParam.rangeR;
  hKs->ksIMap[1].eotfParam.rangeR = 1 / hKs->ksIMap[1].eotfParam.rangeR;

  hKs->ksIMap[0].eotfParam.eotf = EOTF_C2K(eotf);
  hKs->ksIMap[1].eotfParam.eotf = EOTF_C2K(eotf);

  if (eotf == CEotfBt1886)
  {
    GetEotfParamsFc(hKs->ksIMap[0].eotfParam.min, hKs->ksIMap[0].eotfParam.max,
                    hKs->ksIMap[0].eotfParam.gamma,
                    &(hKs->ksIMap[0].eotfParam.a), &(hKs->ksIMap[0].eotfParam.b));

    GetEotfParamsFc(hKs->ksIMap[1].eotfParam.min, hKs->ksIMap[1].eotfParam.max,
                    hKs->ksIMap[1].eotfParam.gamma,
                    &(hKs->ksIMap[1].eotfParam.a), &(hKs->ksIMap[1].eotfParam.b));
  }
  else if (eotf == CEotfPower)
  {
    hKs->ksIMap[0].eotfParam.a = (FloatComp_t)(hKs->ksIMap[0].eotfParam.max - hKs->ksIMap[0].eotfParam.min);
    hKs->ksIMap[0].eotfParam.b = (FloatComp_t)hKs->ksIMap[0].eotfParam.min;
    hKs->ksIMap[1].eotfParam.a = (FloatComp_t)(hKs->ksIMap[1].eotfParam.max - hKs->ksIMap[1].eotfParam.min);
    hKs->ksIMap[1].eotfParam.b = (FloatComp_t)hKs->ksIMap[1].eotfParam.min;
  }
  else if (eotf == CEotfHlg)
  {
    DeriveHlgEotfParam(hKs->ksIMap[0].eotfParam.min, hKs->ksIMap[0].eotfParam.max, &hKs->ksIMap[0].eotfParam);
    DeriveHlgEotfParam(hKs->ksIMap[1].eotfParam.min, hKs->ksIMap[1].eotfParam.max, &hKs->ksIMap[1].eotfParam);
  }
  else
  {
    //PQ=>L params are const
    //pKsIMap->eotfParam.gamma = 0;
  }

}

/* Init all DEITP params based on inputs & default params */
int InitDEITPCore(DmKsFlt_t *hKs)
{
  double v3Wp[3] = {1, 1, 1};
  double m33Rgb2Lms[3][3];
  int i;
  int inBits = hKs->bpc;

  // range in rgb/lms linear space
//  pKsIMap->eotfParam.min = pSrcSigEnv->Min;
  //pKsIMap->eotfParam.max = pSrcSigEnv->Max;

  if (!IS_CCLR_IPT_ICTCP(hKs->clr_fmt))
  {
    if (!IS_CCLR_RGB_RGBA(hKs->clr_fmt))
    {
      //// Yuv2Rgb m33 xfer
      if (hKs->Yuv2RgbExt)
      {
        /*For now, let us use CYuvXferSpecR601, full range later can be changed as needed */
        //AssignM33D2((const double(*)[3])pSrcSigEnv->M33Yuv2Rgb, pKsIMap->m33Yuv2Rgb);
        GetYuv2RgbM33Fc(CYuvXferSpecR601, &hKs->ksIMap[0].m33Yuv2Rgb);
        GetYuv2RgbM33Fc(CYuvXferSpecR601, &hKs->ksIMap[1].m33Yuv2Rgb);

      }
      else if (hKs->Rng != CRngNarrow)
      {
        GetYuv2RgbM33Fc(hKs->YuvXferSpec, &hKs->ksIMap[0].m33Yuv2Rgb);
        GetYuv2RgbM33Fc(hKs->YuvXferSpec, &hKs->ksIMap[1].m33Yuv2Rgb);
      }
      else
      {
        GetYuv2RgbM33NarrowFc(hKs->YuvXferSpec, hKs->bpc, &hKs->ksIMap[0].m33Yuv2Rgb);
        GetYuv2RgbM33NarrowFc(hKs->YuvXferSpec, hKs->bpc, &hKs->ksIMap[1].m33Yuv2Rgb);
      }

      if (hKs->Yuv2RgbOffNormExt)
      {
        double v3[3];
        v3[0] = V3Yuv2RgbOffNorm[0] * (1 << hKs->bpc);
        v3[1] = V3Yuv2RgbOffNorm[1] * (1 << hKs->bpc);
        v3[2] = V3Yuv2RgbOffNorm[2] * (1 << hKs->bpc);

        AssignV3D2(v3, &hKs->ksIMap[0].v3Yuv2RgbOff);
        AssignV3D2(v3, &hKs->ksIMap[1].v3Yuv2RgbOff);
      }
      else if (hKs->Yuv2RgbOffExt)
      {
        AssignV3D2(V3Yuv2RgbOff, &hKs->ksIMap[0].v3Yuv2RgbOff);
        AssignV3D2(V3Yuv2RgbOff, &hKs->ksIMap[1].v3Yuv2RgbOff);
      }
      else
      {
        GetYuvRgbOffFc(hKs->Rng, hKs->bpc, &hKs->ksIMap[0].v3Yuv2RgbOff);
        GetYuvRgbOffFc(hKs->Rng, hKs->bpc, &hKs->ksIMap[1].v3Yuv2RgbOff);
      }
        
    }
    else
    {
      // just to be sure
      IDENTITY_M33(hKs->ksIMap[0].m33Yuv2Rgb);
      IDENTITY_M33(hKs->ksIMap[1].m33Yuv2Rgb);
      ZERO_V3(hKs->ksIMap[0].v3Yuv2RgbOff);
      ZERO_V3(hKs->ksIMap[1].v3Yuv2RgbOff);
    }

    DeitpUpdateEotfParam(hKs, hKs->eotf, (FloatComp_t)hKs->gamma, hKs->bpc);

    if(hKs->Rgb2LmsM33Ext)
    {
      MemCpyByte(&m33Rgb2Lms[0][0], &M33Rgb2Lms[0][0], 9 * (int)sizeof(m33Rgb2Lms[0][0]));
    }
    else if (hKs->Rgb2LmsRgbwExt)
    {
      GetRgb2LmsByPrimsM33(V8Rgbw[0], V8Rgbw[1], V8Rgbw[2], V8Rgbw[3],
                           V8Rgbw[4], V8Rgbw[5], V8Rgbw[6], V8Rgbw[7], m33Rgb2Lms);
    }
    else
      GetRgb2LmsByDefM33(hKs->RgbDef, m33Rgb2Lms);
  }
  else
  {
    //// for ipt pq input
    double ctSig;

    /*** m33Yuv2Rgb: Ipt2Lms xfer, no range conversion ***/
    if (hKs->clr_fmt == CClrICtCp)
    {
      ctSig = 0.04;
      if (hKs->Rng != CRngNarrow)
      {
        GetIctcpIo2LmsM33Fc(&hKs->ksIMap[0].m33Yuv2Rgb);
        GetIctcpIo2LmsM33Fc(&hKs->ksIMap[1].m33Yuv2Rgb);
      }
      else
      {
        GetIctcpIo2LmsM33NarrowFc(hKs->bpc, &hKs->ksIMap[0].m33Yuv2Rgb);
        GetIctcpIo2LmsM33NarrowFc(hKs->bpc, &hKs->ksIMap[1].m33Yuv2Rgb);
      }
      if(hKs->NonStdItp)
      {
        S_TIMES_V3(2, &hKs->ksIMap[0].m33Yuv2Rgb[1][0]);
        S_TIMES_V3(2, &hKs->ksIMap[1].m33Yuv2Rgb[1][0]);
      }
    }
    else
    {
      ctSig = 0.02;
      if (hKs->Rng != CRngNarrow)
      {
        GetIpt2LmsM33Fc(&hKs->ksIMap[0].m33Yuv2Rgb);
        GetIpt2LmsM33Fc(&hKs->ksIMap[1].m33Yuv2Rgb);
      }
      else
      {
        GetIpt2LmsM33NarrowFc(hKs->bpc, &hKs->ksIMap[0].m33Yuv2Rgb);
        GetIpt2LmsM33NarrowFc(hKs->bpc, &hKs->ksIMap[1].m33Yuv2Rgb);
      }
    }
    /* same as yuv=>rgb */
    GetYuvRgbOffFc(hKs->Rng, hKs->bpc, &hKs->ksIMap[0].v3Yuv2RgbOff);
    GetYuvRgbOffFc(hKs->Rng, hKs->bpc, &hKs->ksIMap[1].v3Yuv2RgbOff);

    /**** eotf ***/
    /* already in full range */
    hKs->ksIMap[0].eotfParam.rangeR = (FloatComp_t)(1.0 / ((1 << hKs->bpc) - 1));
    hKs->ksIMap[1].eotfParam.rangeR = (FloatComp_t)(1.0 / ((1 << hKs->bpc) - 1));
    hKs->ksIMap[0].eotfParam.rangeMin = 0;
    hKs->ksIMap[1].eotfParam.rangeMin = 0;
    /* ipt must be in pq */
    hKs->ksIMap[0].eotfParam.eotf = KEotfPq;
    hKs->ksIMap[1].eotfParam.eotf = KEotfPq;

    /*** PQ=>L: const ***/
    /*** m33Rgb2Lms: crosstalk handling ***/
    DE_CROSSTALK_M33(ctSig, m33Rgb2Lms);

    /*** L=>PQ const ***/
  }

  GetRgb2XyzFrom2LmsM33Fc((const double(*)[3])m33Rgb2Lms, &hKs->ksIMap[0].m33Rgb2Xyz);
  GetRgb2XyzFrom2LmsM33Fc((const double(*)[3])m33Rgb2Lms, hKs->ksIMap[1].m33Rgb2Xyz);

  /* Q&D: since these two are constant, just assign it
* XYZ2RGB = Dolby_getmatrix('xyz2r2020');
* RGB2LMS = [1688,2146,262; 683,2951,462; 99,309,3688]'./4096;
* xyz2lms = XYZ2RGB*RGB2LMS;
* xyz2lms' =
*/
  hKs->m33Xyz2LmsViaR2020[0][0] = 0.359283259012122;
  hKs->m33Xyz2LmsViaR2020[0][1] = 0.697605114777950;
  hKs->m33Xyz2LmsViaR2020[0][2] = -0.035891593232029;

  hKs->m33Xyz2LmsViaR2020[1][0] = -0.192080846370499;
  hKs->m33Xyz2LmsViaR2020[1][1] = 1.100476797037432;
  hKs->m33Xyz2LmsViaR2020[1][2] = 0.075374865851912;

  hKs->m33Xyz2LmsViaR2020[2][0] = 0.007079784460748;
  hKs->m33Xyz2LmsViaR2020[2][1] = 0.074839666218637;
  hKs->m33Xyz2LmsViaR2020[2][2] = 0.843326545389877;

  /*
* lms2IctcpDm' = [2048,2048,0; 6610,-13613,7003; 17933,-17390,-543]'./4096*diag([1,0.5,1]);
*/
  hKs->m33Lms2IctcpDm[0][0] = 0.500000000000000;
  hKs->m33Lms2IctcpDm[0][1] = 0.500000000000000;
  hKs->m33Lms2IctcpDm[0][2] = 0;

  hKs->m33Lms2IctcpDm[1][0] = 0.806884765625000;
  hKs->m33Lms2IctcpDm[1][1] = -1.661743164062500;
  hKs->m33Lms2IctcpDm[1][2] = 0.854858398437500;

  hKs->m33Lms2IctcpDm[2][0] = 4.378173828125000;
  hKs->m33Lms2IctcpDm[2][1] = -4.245605468750000;
  hKs->m33Lms2IctcpDm[2][2] = -0.132568359375000;


  return 0;
}

#ifdef WIN32
#pragma warning ( pop )
#endif
