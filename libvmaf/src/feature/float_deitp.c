/* Copyright Dolby Laboratories Inc >*/


#include <errno.h>
#include <math.h>
#include <string.h>

#include "feature_collector.h"
#include "feature_extractor.h"
#include "picture_copy.h"

#include "mem.h"

#include "KdeType.h"
#include "CdeType.h"
#include "DeApi.h"
#include "CdeUtil.h"

#define Tap_6_FilterUvRowUsHalfSize    3
static const FloatComp_t Tap_6_FilterUvRowUs0_m[Tap_6_FilterUvRowUsHalfSize<<1] ={
    (FloatComp_t)(2/256.0), -(FloatComp_t)(12/256.0), (FloatComp_t)(65/256.0),
    (FloatComp_t)(222/256.0), -(FloatComp_t)(25/256.0), (FloatComp_t)(4/256.0)
};

static const FloatComp_t Tap_6_FilterUvRowUs1_m[Tap_6_FilterUvRowUsHalfSize<<1] ={
    (FloatComp_t)(4/256.0), -(FloatComp_t)(25/256.0), (FloatComp_t)(222/256.0),
    (FloatComp_t)(65/256.0), -(FloatComp_t)(12/256.0), (FloatComp_t)(2/256.0)
};

#define FilterUvColUsHalfSize    4
static const FloatComp_t FilterUvColUs_m[FilterUvColUsHalfSize<<1] ={
    (FloatComp_t)(22/4096.0), (FloatComp_t)(94/4096.0), -(FloatComp_t)(524/4096.0), (FloatComp_t)(2456/4096.0),
    (FloatComp_t)(2456/4096.0), -(FloatComp_t)(524/4096.0), (FloatComp_t)(94/4096.0), (FloatComp_t)(22/4096.0)
};

#define CHRM_C2K(c) (  ((c) == CChrm420) ? KChrm420 :          \
                       ((c) == CChrm422) ? KChrm422 :  KChrm444 )

#define DTP_C2K(c) (  ((c) == CDtpU16) ? KDtpU16 :          \
                     ((c) == CDtpU8)  ? KDtpU8  :  KDtpF32 )

#define WEAV_C2K(c) (  ((c) == CWeavPlnr) ? KWeavPlnr :            \
                      ((c) == CWeavIntl) ? KWeavIntl :  KWeavUyVy )

#define LOC_C2K(c)  (  ((c) == CLocHost) ?  KLocHost : KLocDev  )


typedef struct DeitpState {

    DmKsFlt_t dEITPParams;

    /* Optional input parameters can be set from configuration */
    int eotf;
    int clrFmt; //YUV, RGB, IPT, ICtCp
    int RgbDef;
    double gamma;
    double minLum; 
    double maxLum;
    bool Yuv2RgbExt;
    int YuvXferSpec;
    int Rng;
    

 /* Create reference picture & distorted picture buffers
    each picture can accomidate YUV 444 16 bit pixels
    */
    VmafPicture ref;
    VmafPicture dis;
    /* Temporary buffers used to hold intermediate DE 
       scores for each pixel location
       temp_buf also used for upsampling */
    FloatComp_t *temp_buf; 
    FloatComp_t *temp_buf1;

} DeitpState;

static const VmafOption options[] = {
    {
        .name = "eotf",
        .help = "eoft mode: set pq for ipt, default is bt1886",
        .offset = offsetof(DeitpState, eotf),
        .type = VMAF_OPT_TYPE_INT,
        .default_val.i = CEotfBt1886,
        .min = 0,
        .max = CEotfNum,
    },
    {
        .name = "clrFmt",
        .help = "input color format: DEITP can accept YUV, RGB, Ipt, ICtCp",
        .offset = offsetof(DeitpState, clrFmt),
        .type = VMAF_OPT_TYPE_INT,
        .default_val.i = CClrYuv,
        .min = CClrYuv,
        .max = CClrNum,
    },
    {
        .name = "gamma",
        .help = "gamma value used to conver YUV to XYZ format, default is set to 2.4 (bt1886)",
        .offset = offsetof(DeitpState, gamma),
        .type = VMAF_OPT_TYPE_DOUBLE,
        .default_val.d = 2.4,
        .min = 0,
        .max = 5.0, /*Set some peak value*/
    },
    {
        .name = "minLum",
        .help = "Minimum Luminance Value ",
        .offset = offsetof(DeitpState, minLum),
        .type = VMAF_OPT_TYPE_DOUBLE,
        .default_val.d = 0.005,
        .min = 0.005,
        .max = 100.0, /*Set some peak value*/
    },
    {
        .name = "maxLum",
        .help = "Maximum Luminance Value ",
        .offset = offsetof(DeitpState, maxLum),
        .type = VMAF_OPT_TYPE_DOUBLE,
        .default_val.d = 100.0,
        .min = 0.005,
        .max = 100.0, /*Set some peak value*/
    },
    {
        .name = "Yuv2RgbExt",
        .help = "Enable User defined YUV to RGB Matrices. At present there is no method to give as input, \
                  but need to hard coded in the code.",
        .offset = offsetof(DeitpState, Yuv2RgbExt),
        .type = VMAF_OPT_TYPE_BOOL,
        .default_val.b = false,
    },
    {
        .name = "YuvXferSpec",
        .help = "YUV to RGB specification, default is CYuvXferSpecR709.",
        .offset = offsetof(DeitpState, YuvXferSpec),
        .type = VMAF_OPT_TYPE_INT,
        .default_val.i = CYuvXferSpecR709,
        .min = CYuvXferSpecP3d65,
        .max = CYuvXferSpecNum,
    },
    {
        .name = "Rng",
        .help = "Rgb to YUV range CRngFull/Narrow/Sdi, default is Narrow.",
        .offset = offsetof(DeitpState, Rng),
        .type = VMAF_OPT_TYPE_INT,
        .default_val.i = CRngNarrow,
        .min = CRngNarrow,
        .max = CRngNum,
    },
    {
        .name = "RgbDef",
        .help = "Color space defenition for RGB->LMS,default is 709.",
        .offset = offsetof(DeitpState, RgbDef),
        .type = VMAF_OPT_TYPE_INT,
        .default_val.i = CRgbDefR709,
        .min = CRgbDefP3d65,
        .max = CRgbDefNum,
    },

    {NULL}};

static int init(VmafFeatureExtractor *fex, enum VmafPixelFormat pix_fmt,
                unsigned bpc, unsigned w, unsigned h)
{
  int err = 0;
  //printf(" DEITP Init  function called \n" );

  DeitpState *s = fex->priv;
  HDeKsFlt_t hKs = &(s->dEITPParams);

  int stride = w * sizeof(FloatComp_t);
  int buf_size = stride * h; // per channel buffer
  
  hKs->rowPitchNum = w;
  hKs->frmBuf0[0] = (FloatComp_t *)aligned_malloc(buf_size, 32); // ref buf 
  hKs->frmBuf1[0] = (FloatComp_t *)aligned_malloc(buf_size, 32); //ref_buf_cb
  hKs->frmBuf2[0] = (FloatComp_t *)aligned_malloc(buf_size, 32); //ref_buf_cr
  hKs->frmBuf0[1] = (FloatComp_t *)aligned_malloc(buf_size, 32); //dis_buf
  hKs->frmBuf1[1] = (FloatComp_t *)aligned_malloc(buf_size, 32); //dis_buf_cb
  hKs->frmBuf2[1] = (FloatComp_t *)aligned_malloc(buf_size, 32); //dis_buf_cr
  s->temp_buf     = (FloatComp_t *)aligned_malloc(buf_size, 32);
  s->temp_buf1    = (FloatComp_t *)aligned_malloc(buf_size, 32);

  if (!hKs->frmBuf0[0] || !hKs->frmBuf1[0] || !hKs->frmBuf2[0] ||
      !hKs->frmBuf0[1] || !hKs->frmBuf1[1] || !hKs->frmBuf2[1] ||
      !s->temp_buf || !s->temp_buf1)
      goto fail;

  if (pix_fmt == VMAF_PIX_FMT_YUV420P)
  {
      hKs->ksUds[0].chrm = CHRM_C2K(CChrm420);
      hKs->ksUds[1].chrm = CHRM_C2K(CChrm420);
  }
  else if (pix_fmt == VMAF_PIX_FMT_YUV422P)
  {
      hKs->ksUds[0].chrm = CHRM_C2K(CChrm422);
      hKs->ksUds[1].chrm = CHRM_C2K(CChrm422);
  }
  else if(pix_fmt == VMAF_PIX_FMT_YUV444P)
  {
    hKs->ksUds[0].chrm = CHRM_C2K(CChrm444);
    hKs->ksUds[1].chrm = CHRM_C2K(CChrm444);
  }
  else
  {
    goto fail; // Other picture formats are not supported
  }
  
  
  hKs->ksUds[0].minUs = 0;
  hKs->ksUds[0].maxUs = (FloatComp_t)((1 << bpc) - 1);
  hKs->ksUds[0].filterUvRowUsHalfSize = Tap_6_FilterUvRowUsHalfSize;
  hKs->ksUds[0].filterUvRowUs0 = (FloatComp_t *)Tap_6_FilterUvRowUs0_m;
  hKs->ksUds[0].filterUvRowUs1 = (FloatComp_t *)Tap_6_FilterUvRowUs1_m;
  hKs->ksUds[0].filterUvColUsHalfSize = FilterUvColUsHalfSize;
  hKs->ksUds[0].filterUvColUs = (FloatComp_t *)FilterUvColUs_m;

  hKs->ksUds[1].minUs = 0;
  hKs->ksUds[1].maxUs = (FloatComp_t)((1 << bpc) - 1);
  hKs->ksUds[1].filterUvRowUsHalfSize = Tap_6_FilterUvRowUsHalfSize;
  hKs->ksUds[1].filterUvRowUs0 = (FloatComp_t *)Tap_6_FilterUvRowUs0_m;
  hKs->ksUds[1].filterUvRowUs1 = (FloatComp_t *)Tap_6_FilterUvRowUs1_m;
  hKs->ksUds[1].filterUvColUsHalfSize = FilterUvColUsHalfSize;
  hKs->ksUds[1].filterUvColUs = (FloatComp_t *)FilterUvColUs_m;

  /*When input is 8 bit, assume data type is always 1 byte or U8 */
  /*VMAF supports only Planar as of now */
  hKs->ksFrmFmt[0].rowPitch = w * sizeof(FloatComp_t);
  hKs->ksFrmFmt[1].rowPitch = w * sizeof(FloatComp_t);
  hKs->ksFrmFmt[0].colPitch = sizeof(FloatComp_t);
  hKs->ksFrmFmt[1].colPitch = sizeof(FloatComp_t);

  if (bpc == 8)
  {
    hKs->ksFrmFmt[0].dtp = DTP_C2K(CDtpU8);
    hKs->ksFrmFmt[1].dtp = DTP_C2K(CDtpU8);
  }
  else
  {
    hKs->ksFrmFmt[0].dtp = DTP_C2K(CDtpU16);
    hKs->ksFrmFmt[1].dtp = DTP_C2K(CDtpU16);
  }

  if (pix_fmt == VMAF_PIX_FMT_YUV420P)
  {
    hKs->ksFrmFmt[0].rowPitchC = hKs->ksFrmFmt[0].rowPitch >> 1;
    hKs->ksFrmFmt[1].rowPitchC = hKs->ksFrmFmt[1].rowPitch >> 1;
  }
  else
  {
    hKs->ksFrmFmt[0].rowPitchC = hKs->ksFrmFmt[0].rowPitch;
    hKs->ksFrmFmt[1].rowPitchC = hKs->ksFrmFmt[1].rowPitch;
  }
  
  hKs->ksFrmFmt[0].rowNum = h;
  hKs->ksFrmFmt[0].colNum = w;
  hKs->ksFrmFmt[0].weav = WEAV_C2K(CWeavPlnr); 
  hKs->ksFrmFmt[0].loc = LOC_C2K(CLocHost);

  hKs->ksFrmFmt[1].rowNum = h;
  hKs->ksFrmFmt[1].colNum = w;
  hKs->ksFrmFmt[1].weav = WEAV_C2K(CWeavPlnr);
  hKs->ksFrmFmt[1].loc = LOC_C2K(CLocHost);

  /* Set default & input params to DEITP context structure */
  hKs->eotf = (CEotf_t)s->eotf;
  hKs->gamma = (FloatComp_t)s->gamma;
  hKs->ksIMap[0].eotfParam.min = (FloatComp_t)s->minLum;
  hKs->ksIMap[0].eotfParam.max = (FloatComp_t)s->maxLum;
  hKs->ksIMap[1].eotfParam.min = (FloatComp_t)s->minLum;
  hKs->ksIMap[1].eotfParam.max = (FloatComp_t)s->maxLum;
  hKs->ksIMap[0].eotfParam.gamma = (FloatComp_t)s->gamma;
  hKs->ksIMap[1].eotfParam.gamma = (FloatComp_t)s->gamma;

  hKs->ksIMap[0].eotfParam.eotf = EOTF_C2K(s->eotf);
  hKs->ksIMap[1].eotfParam.eotf = EOTF_C2K(s->eotf);
  hKs->clr_fmt = (CClr_t)s->clrFmt;
  hKs->bpc = bpc;
  hKs->Yuv2RgbExt = (int)s->Yuv2RgbExt;
  hKs->YuvXferSpec = (CYuvXferSpec_t)s->YuvXferSpec;
  hKs->Rng = (CRng_t)s->Rng;
  hKs->RgbDef = (CRgbDef_t)s->RgbDef;
  hKs->Yuv2RgbOffNormExt = 0; //for now, no external offsets
  hKs->Yuv2RgbOffExt = 0;     //for now, no external offsets
  hKs->Rgb2LmsM33Ext = 0;     //for now, no external offsets
  hKs->Rgb2LmsRgbwExt = 0;    //for now, no external offsets
  hKs->NonStdItp = 0;


  err = InitDEITPCore(hKs);

  hKs->pxlNumR = w * h;
  hKs->pxlNumR = w * h;
  hKs->pxlNumR = 1 / hKs->pxlNumR;


  if(err)
    goto fail;
  else
  {
      return 0;
  }

  fail:
    free(hKs->frmBuf0[0]);
    free(hKs->frmBuf0[1]);
    free(hKs->frmBuf1[0]);
    free(hKs->frmBuf1[1]);
    free(hKs->frmBuf2[0]);
    free(hKs->frmBuf2[1]);
    free(s->temp_buf);
    free(s->temp_buf1);
    return err;

}

static int extract(VmafFeatureExtractor *fex,
                   VmafPicture *ref_pic, VmafPicture *ref_pic_90,
                   VmafPicture *dist_pic, VmafPicture *dist_pic_90,
                   unsigned index, VmafFeatureCollector *feature_collector)
{
  int err = 0;

  //printf(" DEITP Extract function called \n" );

  DeitpState *s = fex->priv;
  HDeKsFlt_t hKs = &(s->dEITPParams);
  
  size_t float_stride = ALIGN_CEIL(hKs->ksFrmFmt[0].colNum * sizeof(float));

  picture_copy_no_scale(hKs->frmBuf0[0], float_stride, ref_pic, 0, ref_pic->bpc);
  picture_copy_no_scale(hKs->frmBuf0[1], float_stride, dist_pic, 0, dist_pic->bpc);

  if (ref_pic->pix_fmt == VMAF_PIX_FMT_YUV420P)
  {
    picture_copy_chroma_no_scale(hKs->frmBuf1[0], ref_pic, 0, ref_pic->bpc, VMAF_Cb);
    picture_copy_chroma_no_scale(hKs->frmBuf1[1], dist_pic, 0, dist_pic->bpc, VMAF_Cb);

    picture_copy_chroma_no_scale(hKs->frmBuf2[0], ref_pic, 0, ref_pic->bpc, VMAF_Cr);
    picture_copy_chroma_no_scale(hKs->frmBuf2[1], dist_pic, 0, dist_pic->bpc, VMAF_Cr);

    YUV420ToYUV444InU8(hKs->frmBuf1[0], hKs->frmBuf2[0], s->temp_buf, ref_pic->w[0], ref_pic->h[0], ref_pic->w[0], ref_pic->bpc);
    YUV420ToYUV444InU8(hKs->frmBuf1[1], hKs->frmBuf2[1], s->temp_buf, dist_pic->w[0], dist_pic->h[0], dist_pic->w[0], dist_pic->bpc);

  }else if (ref_pic->pix_fmt == VMAF_PIX_FMT_YUV422P)
  {
    picture_copy_chroma_no_scale(s->temp_buf, ref_pic, 0, ref_pic->bpc, VMAF_Cb);
    YUV422ToYUV444InU8(hKs->frmBuf1[0], s->temp_buf, ref_pic->w[0], ref_pic->h[0], ref_pic->w[0], ref_pic->bpc);

    picture_copy_chroma_no_scale(s->temp_buf, ref_pic, 0, ref_pic->bpc, VMAF_Cr);
    YUV422ToYUV444InU8(hKs->frmBuf2[0], s->temp_buf, ref_pic->w[0], ref_pic->h[0], ref_pic->w[0], ref_pic->bpc);

    picture_copy_chroma_no_scale(s->temp_buf, dist_pic, 0, dist_pic->bpc, VMAF_Cb);
    YUV422ToYUV444InU8(hKs->frmBuf1[1], s->temp_buf, dist_pic->w[0], dist_pic->h[0], dist_pic->w[0], dist_pic->bpc);

    picture_copy_chroma_no_scale(s->temp_buf, dist_pic, 0, dist_pic->bpc, VMAF_Cr);
    YUV422ToYUV444InU8(hKs->frmBuf2[1], s->temp_buf, dist_pic->w[0], dist_pic->h[0], dist_pic->w[0], dist_pic->bpc);
  }
  else if(ref_pic->pix_fmt == VMAF_PIX_FMT_YUV444P)
  {
    picture_copy_no_scale(hKs->frmBuf0[0], float_stride, ref_pic, 0, ref_pic->bpc);
    picture_copy_no_scale(hKs->frmBuf0[1], float_stride, dist_pic, 0, dist_pic->bpc);

    picture_copy_chroma_no_scale(hKs->frmBuf1[0], ref_pic, 0, ref_pic->bpc, VMAF_Cb);
    picture_copy_chroma_no_scale(hKs->frmBuf1[1], dist_pic, 0, dist_pic->bpc, VMAF_Cb);

    picture_copy_chroma_no_scale(hKs->frmBuf2[0], ref_pic, 0, ref_pic->bpc, VMAF_Cr);
    picture_copy_chroma_no_scale(hKs->frmBuf2[1], dist_pic, 0, dist_pic->bpc, VMAF_Cr);

  } 
  else
  {
    printf("This format is not yet supported \n");
  }

  /*Measure mean, max, standard deviation DE ITP and DE TP scores */
  float meanDe[2];
  float maxDe[2];
  float sdDe[2]; // Standard Deviation

  err = DolbyDe(hKs, meanDe, maxDe, sdDe, s->temp_buf, s->temp_buf1);
  if (err) goto fail;

 // printf("%6d %11s %-7.03f %-7.03f %-7.03f\n", index, "no sm", meanDe, maxDe[0], sdDe);
  printf("%6d %11s %-7.03f %-7.03f\n", index, "no sm", meanDe[0], maxDe[0], sdDe[0]);

  err = vmaf_feature_collector_append(feature_collector, "deitp_max", maxDe[0], index);
  if (err) goto fail;
  err = vmaf_feature_collector_append(feature_collector, "deitp_mean", meanDe[0], index);
  if (err) goto fail;
  err = vmaf_feature_collector_append(feature_collector, "deitp_sd", sdDe[0], index);

  err = vmaf_feature_collector_append(feature_collector, "deTP_max", maxDe[1], index);
  if (err) goto fail;
  err = vmaf_feature_collector_append(feature_collector, "deTP_mean", meanDe[1], index);
  if (err) goto fail;
  err = vmaf_feature_collector_append(feature_collector, "deTP_sd", sdDe[1], index);

  if (err)
    goto fail;

  return 0;

  fail:
    return err;
}


static int close(VmafFeatureExtractor *fex)
{
    DeitpState *s = fex->priv;
    HDeKsFlt_t hKs = &(s->dEITPParams);


    if (hKs->frmBuf0[0]) aligned_free(hKs->frmBuf0[0]);
    if (hKs->frmBuf0[1]) aligned_free(hKs->frmBuf0[1]);
    if (hKs->frmBuf1[0]) aligned_free(hKs->frmBuf1[0]);
    if (hKs->frmBuf1[1]) aligned_free(hKs->frmBuf1[1]);
    if (hKs->frmBuf2[0]) aligned_free(hKs->frmBuf2[0]);
    if (hKs->frmBuf2[1]) aligned_free(hKs->frmBuf2[1]);
    if (s->temp_buf) aligned_free(s->temp_buf);
    if (s->temp_buf1) aligned_free(s->temp_buf1);

}

static const char *provided_features[] = {
    "deitp_mean", "deitp_max", "deitp_sd",
    "deTP_mean", "deTP_max", "deTP_sd",
    NULL};

VmafFeatureExtractor vmaf_fex_float_deitp = {
    .name = "deitp",
    .init = init,
    .extract = extract,
    .close = close,
    .priv_size = sizeof(DeitpState),
    .provided_features = provided_features,
    .options = options,
    //.flags = VMAF_FEATURE_EXTRACTOR_TEMPORAL,
};
