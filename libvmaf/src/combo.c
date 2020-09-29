/**
 *
 *  Copyright 2016-2020 Netflix, Inc.
 *
 *     Licensed under the BSD+Patent License (the "License");
 *     you may not use this file except in compliance with the License.
 *     You may obtain a copy of the License at
 *
 *         https://opensource.org/licenses/BSDplusPatent
 *
 *     Unless required by applicable law or agreed to in writing, software
 *     distributed under the License is distributed on an "AS IS" BASIS,
 *     WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 *     See the License for the specific language governing permissions and
 *     limitations under the License.
 *
 */

#include <limits.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "mem.h"
#include "motion_tools.h"
#include "common/convolution.h"
#include "common/convolution_internal.h"
#include "iqa/ssim_tools.h"
#include "darray.h"
#include "vif_options.h"
#include "adm_options.h"
#include "combo.h"
#include "debug.h"
#include "psnr_tools.h"
#include "offset.h"

#include "KdeType.h"
#include "CdeType.h"
#include "DeApi.h"
#include "CdeUtil.h"

#include "common/blur_array.h"
#include "cpu_info.h"

#define convolution_f32_c  convolution_f32_c_s
#define offset_image       offset_image_s
#define FILTER_5           FILTER_5_s
int compute_adm(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double border_factor, double adm_enhn_gain_limit);
#ifdef COMPUTE_ANSNR
int compute_ansnr(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_psnr, double peak, double psnr_max);
#endif
int compute_vif(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double *score_num, double *score_den, double *scores, double vif_enhn_gain_limit);
int compute_motion(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score);
int compute_psnr(const float *ref, const float *dis, int w, int h, int ref_stride, int dis_stride, double *score, double peak, double psnr_max);
int compute_ssim(const float *ref, const float *cmp, int w, int h, int ref_stride, int cmp_stride, double *score, double *l_score, double *c_score, double *s_score);
int compute_ms_ssim(const float *ref, const float *cmp, int w, int h, int ref_stride, int cmp_stride, double *score, double* l_scores, double* c_scores, double* s_scores);

#define MIN(x, y) (((x) < (y)) ? (x) : (y))


/* Below is updated even in read_frame.c file we need to have common place holder*/
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

#define CHRM_C2K(c) (  ((c) == CChrm420) ? KChrm420 :          \
                       ((c) == CChrm422) ? KChrm422 :  KChrm444 )

#define DTP_C2K(c) (  ((c) == CDtpU16) ? KDtpU16 :          \
                     ((c) == CDtpU8)  ? KDtpU8  :  KDtpF32 )

#define WEAV_C2K(c) (  ((c) == CWeavPlnr) ? KWeavPlnr :            \
                      ((c) == CWeavIntl) ? KWeavIntl :  KWeavUyVy )

#define EOTF_C2K(c) (  ((c) == CEotfBt1886) ?  KEotfBt1886 :              \
                      ((c) == CEotfPq)     ?  KEotfPq     : KEotfPower  )

#define LOC_C2K(c)  (  ((c) == CLocHost) ?  KLocHost : KLocDev  )

void* combo_threadfunc(void* vmaf_thread_data)
{
    // this is our shared thread data
    VMAF_THREAD_STRUCT* thread_data = (VMAF_THREAD_STRUCT*)vmaf_thread_data;

    // set the local variables from the thread shared data
    size_t data_sz = thread_data->data_sz;
    int stride = thread_data->stride;
    double peak = thread_data->peak;
    double psnr_max = thread_data->psnr_max;
    int w = thread_data->w;
    int h = thread_data->h;
    char* errmsg = thread_data->errmsg;
    void* user_data = thread_data->user_data;
    const char* fmt = thread_data->fmt;
    int n_subsample = thread_data->n_subsample;

    double score = 0;
    double score2 = 0;
    double scores[4*2];
    double score_num = 0;
    double score_den = 0;
    double l_score = 0, c_score = 0, s_score = 0;
    double l_scores[SCALES], c_scores[SCALES], s_scores[SCALES];

#ifdef COMPUTE_ANSNR
    double score_psnr = 0;
#endif

    float *ref_buf = 0;
    float *dis_buf = 0;
    float *prev_blur_buf = 0;
    float *blur_buf = 0;
    float *next_ref_buf = 0;
    float *next_dis_buf = 0;
    float *next_blur_buf = 0;
    float *temp_buf = 0;

    /* buffer pointers for color channels*/
    float *ref_buf_no_conv = 0;
    float *next_ref_buf_no_conv = 0;
    float *ref_buf_cb = 0;
    float *ref_buf_cr = 0;
    float *dis_buf_cb = 0;
    float *dis_buf_cr = 0;
    float *next_ref_buf_cb = 0;
    float *next_ref_buf_cr = 0;
    float *next_dis_buf_cb = 0;
    float *next_dis_buf_cr = 0;

    int ret = 0;
    bool next_frame_read;

    bool offset_flag = false;

    // use temp_buf for convolution_f32_c, and fread u and v
    if (!(temp_buf = aligned_malloc(data_sz * 2, MAX_ALIGN)))
    {
        sprintf(errmsg, "aligned_malloc failed for temp_buf.\n");
        goto fail_or_end;
    }

    int frm_idx = -1;

    while (1)
    {

        pthread_mutex_lock(&thread_data->mutex_readframe);

        if (thread_data->stop_threads)
        {
            // this is the signal that another thread has reached the end of the input file, so we all quit
            pthread_mutex_unlock(&thread_data->mutex_readframe);
            goto fail_or_end;
        }
        // the next frame

        frm_idx = thread_data->frm_idx;
        thread_data->frm_idx++;

        if (frm_idx == 0)
        {
            // Allocating the free buffers from buffer array
            blur_buf    = get_free_blur_buf_slot(&thread_data->blur_buf_array, frm_idx);
            ref_buf     = get_free_blur_buf_slot(&thread_data->ref_buf_array, frm_idx);
            dis_buf     = get_free_blur_buf_slot(&thread_data->dis_buf_array, frm_idx);

            if((NULL == blur_buf) || (NULL == ref_buf) || (NULL == dis_buf))
            {
                thread_data->stop_threads = 1;
                sprintf(errmsg, "No free slot found for buffer allocation.\n");
                pthread_mutex_unlock(&thread_data->mutex_readframe);
                goto fail_or_end;
            }

            /* Get extra buffers when DEITP or DESITP is enabled */
            if (thread_data->dlb_mean_deitp_array != NULL)
            {
                ref_buf_no_conv = get_free_blur_buf_slot(&thread_data->ref_buf_no_conv_array, frm_idx);
                ref_buf_cb = get_free_blur_buf_slot(&thread_data->ref_buf_cb_array, frm_idx);
                ref_buf_cr = get_free_blur_buf_slot(&thread_data->ref_buf_cr_array, frm_idx);
                dis_buf_cb = get_free_blur_buf_slot(&thread_data->dis_buf_cb_array, frm_idx);
                dis_buf_cr = get_free_blur_buf_slot(&thread_data->dis_buf_cr_array, frm_idx);

                if ((ref_buf_cb == NULL) || (ref_buf_cr == NULL) || (dis_buf_cb == NULL) || (ref_buf_cb == NULL) || (ref_buf_no_conv == NULL))
                {
                    thread_data->stop_threads = 1;
                    sprintf(errmsg, "No free slot found for buffer allocation.\n");
                    pthread_mutex_unlock(&thread_data->mutex_readframe);
                    goto fail_or_end;
                }
            }

                // read frame from file
            if (thread_data->dlb_mean_deitp_array != NULL)
            {
                // read frame from file and convert to 444
              //  ret = thread_data->read_frame_cb_cr(ref_buf, dis_buf, ref_buf_cb, dis_buf_cb,
                //                                    ref_buf_cr, dis_buf_cr, stride, user_data);
                ret = thread_data->read_frame_cb_cr(ref_buf, dis_buf, ref_buf_cb, dis_buf_cb,
                                                    ref_buf_cr, dis_buf_cr, temp_buf, stride, user_data);
                //copy entire ref buf to extra buffer for deitp
                int size = h * stride;
                memcpy(ref_buf_no_conv,ref_buf,size);
            }
            else
            {
                ret = thread_data->read_frame(ref_buf, dis_buf, temp_buf, stride, user_data);
            }

            if (ret == 1)
            {
                thread_data->stop_threads = 1;
                pthread_mutex_unlock(&thread_data->mutex_readframe);
                goto fail_or_end;
            }
            if (ret == 2)
            {
                thread_data->stop_threads = 1;
                pthread_mutex_unlock(&thread_data->mutex_readframe);
                goto fail_or_end;
            }

            // ===============================================================
            // offset pixel by OPT_RANGE_PIXEL_OFFSET
            // ===============================================================
            offset_image(ref_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_image(dis_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);

            // ===============================================================
            // filter
            // apply filtering (to eliminate effects film grain)
            // stride input to convolution_f32_c is in terms of (sizeof(float) bytes)
            // since stride = ALIGN_CEIL(w * sizeof(float)), stride divides sizeof(float)
            // ===============================================================
            convolution_f32_c(FILTER_5, 5, ref_buf, blur_buf, temp_buf, w, h, stride / sizeof(float), stride / sizeof(float));

        }
        else
        {
            // retrieve from buffer array
            ref_buf     = get_blur_buf(&thread_data->ref_buf_array, frm_idx);
            dis_buf     = get_blur_buf(&thread_data->dis_buf_array, frm_idx);
            blur_buf    = get_blur_buf(&thread_data->blur_buf_array, frm_idx);

            if((NULL == ref_buf) || (NULL == dis_buf) || (NULL == blur_buf))
            {
                thread_data->stop_threads = 1;
                sprintf(errmsg, "Data not available.\n");
                pthread_mutex_unlock(&thread_data->mutex_readframe);
                goto fail_or_end;
            }

            if (thread_data->dlb_mean_deitp_array != NULL)
            {
                ref_buf_no_conv = get_blur_buf(&thread_data->ref_buf_no_conv_array, frm_idx);
                ref_buf_cb = get_blur_buf(&thread_data->ref_buf_cb_array, frm_idx);
                ref_buf_cr = get_blur_buf(&thread_data->ref_buf_cr_array, frm_idx);
                dis_buf_cb = get_blur_buf(&thread_data->dis_buf_cb_array, frm_idx);
                dis_buf_cr = get_blur_buf(&thread_data->dis_buf_cr_array, frm_idx);

                if ((ref_buf_cb == NULL) || (ref_buf_cr == NULL) || (dis_buf_cb == NULL) || (ref_buf_cb == NULL) || (ref_buf_no_conv == NULL)) {
                    thread_data->stop_threads = 1;
                    sprintf(errmsg, "Data not available.\n");
                    pthread_mutex_unlock(&thread_data->mutex_readframe);
                    goto fail_or_end;
                }
            }
        }

        // Allocate free buffer from the buffer array for next frame index
        next_ref_buf 	= get_free_blur_buf_slot(&thread_data->ref_buf_array, frm_idx + 1);
        next_dis_buf 	= get_free_blur_buf_slot(&thread_data->dis_buf_array, frm_idx + 1);
        if((NULL == next_ref_buf) || (NULL == next_dis_buf))
        {
            thread_data->stop_threads = 1;
            sprintf(errmsg, "No free slot found for next buffer.\n");
            pthread_mutex_unlock(&thread_data->mutex_readframe);
            goto fail_or_end;
        }

        if (thread_data->dlb_mean_deitp_array != NULL)
        {
            next_ref_buf_no_conv = get_free_blur_buf_slot(&thread_data->ref_buf_no_conv_array, frm_idx + 1);
            next_ref_buf_cb = get_free_blur_buf_slot(&thread_data->ref_buf_cb_array, frm_idx + 1);
            next_ref_buf_cr = get_free_blur_buf_slot(&thread_data->ref_buf_cr_array, frm_idx + 1);
            next_dis_buf_cb = get_free_blur_buf_slot(&thread_data->dis_buf_cb_array, frm_idx + 1);
            next_dis_buf_cr = get_free_blur_buf_slot(&thread_data->dis_buf_cr_array, frm_idx + 1);

            if ((NULL == next_ref_buf_cb) || (NULL == next_ref_buf_cr) || (NULL == next_dis_buf_cb) || (NULL == next_dis_buf_cr) || (next_ref_buf_no_conv == NULL))
            {
                thread_data->stop_threads = 1;
                sprintf(errmsg, "No free slot found for next buffer.\n");
                pthread_mutex_unlock(&thread_data->mutex_readframe);
                goto fail_or_end;
            }
        }
        if (thread_data->dlb_mean_deitp_array != NULL)
        {
            // read frame from file and convert to 444
            ret = thread_data->read_frame_cb_cr(next_ref_buf, next_dis_buf, next_ref_buf_cb, next_dis_buf_cb,
                                                next_ref_buf_cr, next_dis_buf_cr, temp_buf, stride, user_data);
            //copy entire ref buf to extra buffer for deitp
            int size = h * stride;
            memcpy(next_ref_buf_no_conv, next_ref_buf, size);
        }
        else
        {
            ret = thread_data->read_frame(next_ref_buf, next_dis_buf, temp_buf, stride, user_data);
        }
        if (ret == 1)
        {
            thread_data->stop_threads = 1;
            pthread_mutex_unlock(&thread_data->mutex_readframe);
            goto fail_or_end;
        }
        if (ret == 2)
        {
            thread_data->stop_threads = 1;
            next_frame_read = false;
        }
        else
        {
            next_frame_read = true;
        }

        if (next_frame_read)
        {
            next_blur_buf     = get_free_blur_buf_slot(&thread_data->blur_buf_array, frm_idx + 1);
            if(NULL == next_blur_buf)
            {
                thread_data->stop_threads = 1;
                sprintf(errmsg, "No free slot found for blur buffer.\n");
                goto fail_or_end;
            }
            // ===============================================================
            // offset pixel by OPT_RANGE_PIXEL_OFFSET
            // ===============================================================
            offset_image(next_ref_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_image(next_dis_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);

            // ===============================================================
            // filter
            // apply filtering (to eliminate effects film grain)
            // stride input to convolution_f32_c is in terms of (sizeof(float) bytes)
            // since stride = ALIGN_CEIL(w * sizeof(float)), stride divides sizeof(float)
            // ===============================================================
            convolution_f32_c(FILTER_5, 5, next_ref_buf, next_blur_buf, temp_buf, w, h, stride / sizeof(float), stride / sizeof(float));

        }

        // release ref and dis buffer references after blur buf computation
        release_blur_buf_reference(&thread_data->ref_buf_array, frm_idx + 1);
        release_blur_buf_reference(&thread_data->dis_buf_array, frm_idx + 1);

        if (thread_data->dlb_mean_deitp_array != NULL)
        {
            release_blur_buf_reference(&thread_data->ref_buf_no_conv_array, frm_idx + 1);
            release_blur_buf_reference(&thread_data->ref_buf_cb_array, frm_idx + 1);
            release_blur_buf_reference(&thread_data->dis_buf_cb_array, frm_idx + 1);
            release_blur_buf_reference(&thread_data->ref_buf_cr_array, frm_idx + 1);
            release_blur_buf_reference(&thread_data->dis_buf_cr_array, frm_idx + 1);
        }

        pthread_mutex_unlock(&thread_data->mutex_readframe);

        dbg_printf("frame: %d, \n", frm_idx);

        // ===============================================================
        // for the PSNR, SSIM and MS-SSIM, offset are 0. Since in prev read
        // step they have been offset by OPT_RANGE_PIXEL_OFFSET, now
        // offset them back.
        // ===============================================================

        // offset back the buffers only if required
        if (frm_idx % n_subsample == 0 && ( (thread_data->psnr_array != NULL) || (thread_data->ssim_array != NULL) || (thread_data->ms_ssim_array != NULL) ))
        {
            offset_image(ref_buf, -OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_image(dis_buf, -OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_flag = true;
		}

        if (frm_idx % n_subsample == 0 && thread_data->psnr_array != NULL)
        {
            /* =========== psnr ============== */
            ret = compute_psnr(ref_buf, dis_buf, w, h, stride, stride, &score, peak, psnr_max);

            if (ret)
            {
                sprintf(errmsg, "compute_psnr failed.\n");
                goto fail_or_end;
            }

            dbg_printf("psnr: %.3f, ", score);

            insert_array_at(thread_data->psnr_array, score, frm_idx);
        }

        if (frm_idx % n_subsample == 0 && thread_data->ssim_array != NULL)
        {

            /* =========== ssim ============== */
            if ((ret = compute_ssim(ref_buf, dis_buf, w, h, stride, stride, &score, &l_score, &c_score, &s_score)))
            {
                sprintf(errmsg, "compute_ssim failed.\n");
                goto fail_or_end;
            }

            dbg_printf("ssim: %.3f, ", score);

            insert_array_at(thread_data->ssim_array, score, frm_idx);
        }

        if (frm_idx % n_subsample == 0 && thread_data->ms_ssim_array != NULL)
        {
            /* =========== ms-ssim ============== */
            if ((ret = compute_ms_ssim(ref_buf, dis_buf, w, h, stride, stride, &score, l_scores, c_scores, s_scores)))
            {
                sprintf(errmsg, "compute_ms_ssim failed.\n");
                goto fail_or_end;
            }

            dbg_printf("ms_ssim: %.3f, ", score);

            insert_array_at(thread_data->ms_ssim_array, score, frm_idx);
        }

        // ===============================================================
        // for the rest, offset pixel by OPT_RANGE_PIXEL_OFFSET
        // ===============================================================

        if (offset_flag)
        {
            offset_image(ref_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_image(dis_buf, OPT_RANGE_PIXEL_OFFSET, w, h, stride);
            offset_flag = false;
		}

        /* =========== adm ============== */
        if (frm_idx % n_subsample == 0)
        {
            if ((ret = compute_adm(ref_buf, dis_buf, w, h, stride, stride, &score, &score_num, &score_den, scores,
                    ADM_BORDER_FACTOR, DEFAULT_ADM_ENHN_GAIN_LIMIT)))
            {
                sprintf(errmsg, "compute_adm failed.\n");
                goto fail_or_end;
            }

            dbg_printf("adm: %.3f, ", score);
            dbg_printf("adm_num: %.3f, ", score_num);
            dbg_printf("adm_den: %.3f, ", score_den);
            dbg_printf("adm_num_scale0: %.3f, ", scores[0]);
            dbg_printf("adm_den_scale0: %.3f, ", scores[1]);
            dbg_printf("adm_num_scale1: %.3f, ", scores[2]);
            dbg_printf("adm_den_scale1: %.3f, ", scores[3]);
            dbg_printf("adm_num_scale2: %.3f, ", scores[4]);
            dbg_printf("adm_den_scale2: %.3f, ", scores[5]);
            dbg_printf("adm_num_scale3: %.3f, ", scores[6]);
            dbg_printf("adm_den_scale3: %.3f, ", scores[7]);

            insert_array_at(thread_data->adm_num_array, score_num, frm_idx);
            insert_array_at(thread_data->adm_den_array, score_den, frm_idx);
            insert_array_at(thread_data->adm_num_scale0_array, scores[0], frm_idx);
            insert_array_at(thread_data->adm_den_scale0_array, scores[1], frm_idx);
            insert_array_at(thread_data->adm_num_scale1_array, scores[2], frm_idx);
            insert_array_at(thread_data->adm_den_scale1_array, scores[3], frm_idx);
            insert_array_at(thread_data->adm_num_scale2_array, scores[4], frm_idx);
            insert_array_at(thread_data->adm_den_scale2_array, scores[5], frm_idx);
            insert_array_at(thread_data->adm_num_scale3_array, scores[6], frm_idx);
            insert_array_at(thread_data->adm_den_scale3_array, scores[7], frm_idx);
        }
#ifdef COMPUTE_ANSNR

        if (frm_idx % n_subsample == 0)
        {

            /* =========== ansnr ============== */
            if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
            {
                // max psnr 60.0 for 8-bit per Ioannis
                ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.0, 60.0);
            }
            else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
            {
                // 10 bit gets normalized to 8 bit, peak is 1023 / 4.0 = 255.75
                // max psnr 72.0 for 10-bit per Ioannis
                ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.75, 72.0);
            }
            else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
            {
                // 12 bit gets normalized to 8 bit, peak is (2^12 - 1) / 16.0 = 255.9375
                // max psnr 84.0 for 12-bit
                ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.9375, 84.0);
            }
            else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
            {
                // 16 bit gets normalized to 8 bit, peak is (2^16 - 1)) / 256.0 = 255.99609375
                // max psnr 108.0 for 16-bit
                ret = compute_ansnr(ref_buf, dis_buf, w, h, stride, stride, &score, &score_psnr, 255.99609375, 108.0);
            }
            else
            {
                sprintf(errmsg, "unknown format %s.\n", fmt);
                goto fail_or_end;
            }
            if (ret)
            {
                sprintf(errmsg, "compute_ansnr failed.\n");
                goto fail_or_end;
            }

            dbg_printf("ansnr: %.3f, ", score);
            dbg_printf("anpsnr: %.3f, ", score_psnr);
        }

#endif

        /* =========== motion ============== */

        if (frm_idx % n_subsample == 0)
        {

            // compute
            if (frm_idx == 0)
            {
                score = 0.0;
                score2 = 0.0;
            }
            else
            {
                // avoid multiple memory copies
                prev_blur_buf = get_blur_buf(&thread_data->blur_buf_array, frm_idx - 1);
                if(NULL == prev_blur_buf)
                {
                    thread_data->stop_threads = 1;
                    sprintf(errmsg, "Data not available for prev_blur_buf.\n");
                    goto fail_or_end;
                }
                if ((ret = compute_motion(prev_blur_buf, blur_buf, w, h, stride, stride, &score)))
                {
                    sprintf(errmsg, "compute_motion (prev) failed.\n");
                    goto fail_or_end;
                }
                release_blur_buf_reference(&thread_data->blur_buf_array, frm_idx - 1);

                if (next_frame_read)
                {
                    if ((ret = compute_motion(blur_buf, next_blur_buf, w, h, stride, stride, &score2)))
                    {
                        sprintf(errmsg, "compute_motion (next) failed.\n");
                        goto fail_or_end;
                    }
                    score2 = MIN(score, score2);
                }
                else
                {
                    score2 = score;
                }
            }

            dbg_printf("motion: %.3f, ", score);
            dbg_printf("motion2: %.3f, ", score2);

            insert_array_at(thread_data->motion_array, score, frm_idx);
            insert_array_at(thread_data->motion2_array, score2, frm_idx);

        }
        /* Indicate that motion score computation for this frame is complete */
        insert_array_at(thread_data->motion_score_compute_flag_array, 1.0, frm_idx);
        release_blur_buf_reference(&thread_data->blur_buf_array, frm_idx + 1);

        /* =========== vif ============== */

        if (frm_idx % n_subsample == 0)
        {
            if ((ret = compute_vif(ref_buf, dis_buf, w, h, stride, stride,
                    &score, &score_num, &score_den, scores, DEFAULT_VIF_ENHN_GAIN_LIMIT)))
            {
                sprintf(errmsg, "compute_vif failed.\n");
                goto fail_or_end;
            }

            // dbg_printf("vif_num: %.3f, ", score_num);
            // dbg_printf("vif_den: %.3f, ", score_den);
            dbg_printf("vif_num_scale0: %.3f, ", scores[0]);
            dbg_printf("vif_den_scale0: %.3f, ", scores[1]);
            dbg_printf("vif_num_scale1: %.3f, ", scores[2]);
            dbg_printf("vif_den_scale1: %.3f, ", scores[3]);
            dbg_printf("vif_num_scale2: %.3f, ", scores[4]);
            dbg_printf("vif_den_scale2: %.3f, ", scores[5]);
            dbg_printf("vif_num_scale3: %.3f, ", scores[6]);
            dbg_printf("vif_den_scale3: %.3f, ", scores[7]);
            dbg_printf("vif: %.3f, ", score);

            insert_array_at(thread_data->vif_num_scale0_array, scores[0], frm_idx);
            insert_array_at(thread_data->vif_den_scale0_array, scores[1], frm_idx);
            insert_array_at(thread_data->vif_num_scale1_array, scores[2], frm_idx);
            insert_array_at(thread_data->vif_den_scale1_array, scores[3], frm_idx);
            insert_array_at(thread_data->vif_num_scale2_array, scores[4], frm_idx);
            insert_array_at(thread_data->vif_den_scale2_array, scores[5], frm_idx);
            insert_array_at(thread_data->vif_num_scale3_array, scores[6], frm_idx);
            insert_array_at(thread_data->vif_den_scale3_array, scores[7], frm_idx);
            insert_array_at(thread_data->vif_array, score, frm_idx);
        }

        /* Dolby Feature extractors */
        /*======================DEITF==================*/
        if ((frm_idx % n_subsample == 0) && (thread_data->dlb_mean_deitp_array != NULL))
        {
            DmKsFlt_t dEITPParams;
            HDeKsFlt_t hKs = &dEITPParams;
            unsigned char *frmBuf0s[2];
            unsigned char *frmBuf1s[2];
            unsigned char *frmBuf2s[2];
            /* Glu logic to use DEITF, eventually we need to change below code
           and use vmaf code
           */
            dEITPParams.rowPitchNum = w;
            dEITPParams.frmBuf0[0] = (FloatComp_t *)ref_buf_no_conv;
            dEITPParams.frmBuf1[0] = (FloatComp_t *)ref_buf_cb;
            dEITPParams.frmBuf2[0] = (FloatComp_t *)ref_buf_cr;
            dEITPParams.frmBuf0[1] = (FloatComp_t *)dis_buf;
            dEITPParams.frmBuf1[1] = (FloatComp_t *)dis_buf_cb;
            dEITPParams.frmBuf2[1] = (FloatComp_t *)dis_buf_cr;

           // dEITPParams.frmBufDe[0] =
           // dEITPParams.frmBufDe[1] =

            dEITPParams.ksUds[0].chrm = CHRM_C2K(CChrm420); //at present only YUV 420
            dEITPParams.ksUds[0].minUs = 0;
            dEITPParams.ksUds[0].maxUs = (FloatComp_t)((1 << 8) - 1);
            dEITPParams.ksUds[0].filterUvRowUsHalfSize = Tap_6_FilterUvRowUsHalfSize;
            dEITPParams.ksUds[0].filterUvRowUs0 = (FloatComp_t *)Tap_6_FilterUvRowUs0_m;
            dEITPParams.ksUds[0].filterUvRowUs1 = (FloatComp_t *)Tap_6_FilterUvRowUs1_m;
            dEITPParams.ksUds[0].filterUvColUsHalfSize = FilterUvColUsHalfSize;
            dEITPParams.ksUds[0].filterUvColUs = (FloatComp_t *)FilterUvColUs_m;
            dEITPParams.ksUds[1].chrm = CHRM_C2K(CChrm420); //at present only YUV 420
            dEITPParams.ksUds[1].minUs = 0;
            dEITPParams.ksUds[1].maxUs = (FloatComp_t)((1 << 8) - 1);
            dEITPParams.ksUds[1].filterUvRowUsHalfSize = Tap_6_FilterUvRowUsHalfSize;
            dEITPParams.ksUds[1].filterUvRowUs0 = (FloatComp_t *)Tap_6_FilterUvRowUs0_m;
            dEITPParams.ksUds[1].filterUvRowUs1 = (FloatComp_t *)Tap_6_FilterUvRowUs1_m;
            dEITPParams.ksUds[1].filterUvColUsHalfSize = FilterUvColUsHalfSize;
            dEITPParams.ksUds[1].filterUvColUs = (FloatComp_t *)FilterUvColUs_m;

            dEITPParams.ksFrmFmt[0].rowNum = h;
            dEITPParams.ksFrmFmt[0].colNum = w;
            dEITPParams.ksFrmFmt[0].dtp = DTP_C2K(CDtpU8);
            dEITPParams.ksFrmFmt[0].weav = WEAV_C2K(CWeavPlnr);
            dEITPParams.ksFrmFmt[0].loc = LOC_C2K(CLocHost);
            dEITPParams.ksFrmFmt[0].rowPitch = w;
            dEITPParams.ksFrmFmt[0].colPitch = GET_PXL_COL_PITCH(CClrYuv, CChrm420, CWeavPlnr, CDtpU8);;
            dEITPParams.ksFrmFmt[0].rowPitchC = w>>1;

            dEITPParams.ksFrmFmt[1].rowNum = h;
            dEITPParams.ksFrmFmt[1].colNum = w;
            dEITPParams.ksFrmFmt[1].dtp = DTP_C2K(CDtpU8);
            dEITPParams.ksFrmFmt[1].weav = WEAV_C2K(CWeavPlnr);
            dEITPParams.ksFrmFmt[1].loc = LOC_C2K(CLocHost);
            dEITPParams.ksFrmFmt[1].rowPitch = w;
            dEITPParams.ksFrmFmt[1].colPitch = GET_PXL_COL_PITCH(CClrYuv, CChrm420, CWeavPlnr, CDtpU8);;
            dEITPParams.ksFrmFmt[1].rowPitchC = w>>1;

            dEITPParams.ksIMap[0].eotfParam.min = 0.0050000000000000001;
            dEITPParams.ksIMap[0].eotfParam.max = 100;
            dEITPParams.ksIMap[0].eotfParam.rangeMin = 0;
            FloatComp_t range = ((1<< 8) - 1);
            dEITPParams.ksIMap[0].eotfParam.rangeR = 1/range;
            dEITPParams.ksIMap[0].eotfParam.eotf = EOTF_C2K(CEotfBt1886);
            dEITPParams.ksIMap[0].eotfParam.a = 96.17006396051292;
            dEITPParams.ksIMap[0].eotfParam.b = 0.016404797534949767;
            dEITPParams.ksIMap[0].eotfParam.gamma = 2.3999999999999999;

            dEITPParams.ksIMap[1].eotfParam.min = 0.0050000000000000001;
            dEITPParams.ksIMap[1].eotfParam.max = 100;
            dEITPParams.ksIMap[1].eotfParam.rangeMin = 0;
            dEITPParams.ksIMap[1].eotfParam.rangeR = 1/range;
            dEITPParams.ksIMap[1].eotfParam.eotf = EOTF_C2K(CEotfBt1886);
            dEITPParams.ksIMap[1].eotfParam.a = 96.17006396051292;
            dEITPParams.ksIMap[1].eotfParam.b = 0.016404797534949767;
            dEITPParams.ksIMap[1].eotfParam.gamma = 2.3999999999999999;

            /* Q&D: since these two are constant, just assign it
     * XYZ2RGB = Dolby_getmatrix('xyz2r2020');
     * RGB2LMS = [1688,2146,262; 683,2951,462; 99,309,3688]'./4096;
     * xyz2lms = XYZ2RGB*RGB2LMS;
     * xyz2lms' =
    */
            dEITPParams.m33Xyz2LmsViaR2020[0][0] = 0.359283259012122;
            dEITPParams.m33Xyz2LmsViaR2020[0][1] = 0.697605114777950;
            dEITPParams.m33Xyz2LmsViaR2020[0][2] = -0.035891593232029;

            dEITPParams.m33Xyz2LmsViaR2020[1][0] = -0.192080846370499;
            dEITPParams.m33Xyz2LmsViaR2020[1][1] = 1.100476797037432;
            dEITPParams.m33Xyz2LmsViaR2020[1][2] = 0.075374865851912;

            dEITPParams.m33Xyz2LmsViaR2020[2][0] = 0.007079784460748;
            dEITPParams.m33Xyz2LmsViaR2020[2][1] = 0.074839666218637;
            dEITPParams.m33Xyz2LmsViaR2020[2][2] = 0.843326545389877;

            /*
    * lms2IctcpDm' = [2048,2048,0; 6610,-13613,7003; 17933,-17390,-543]'./4096*diag([1,0.5,1]);
    */
            dEITPParams.m33Lms2IctcpDm[0][0] = 0.500000000000000;
            dEITPParams.m33Lms2IctcpDm[0][1] = 0.500000000000000;
            dEITPParams.m33Lms2IctcpDm[0][2] = 0;

            dEITPParams.m33Lms2IctcpDm[1][0] = 0.806884765625000;
            dEITPParams.m33Lms2IctcpDm[1][1] = -1.661743164062500;
            dEITPParams.m33Lms2IctcpDm[1][2] = 0.854858398437500;

            dEITPParams.m33Lms2IctcpDm[2][0] = 4.378173828125000;
            dEITPParams.m33Lms2IctcpDm[2][1] = -4.245605468750000;
            dEITPParams.m33Lms2IctcpDm[2][2] = -0.132568359375000;


            double m33Rgb2Lms[3][3];
            GetRgb2LmsByDefM33(CRgbDefR709,m33Rgb2Lms);
            GetRgb2XyzFrom2LmsM33Fc((const double(*)[3])m33Rgb2Lms, &dEITPParams.ksIMap[0].m33Rgb2Xyz);
            GetRgb2XyzFrom2LmsM33Fc((const double(*)[3])m33Rgb2Lms, &dEITPParams.ksIMap[1].m33Rgb2Xyz);

            GetYuv2RgbM33Fc(CYuvXferSpecR709, &dEITPParams.ksIMap[0].m33Yuv2Rgb);
            GetYuv2RgbM33Fc(CYuvXferSpecR709, &dEITPParams.ksIMap[1].m33Yuv2Rgb);
            GetYuvRgbOffFc(CRngFull, 8, &dEITPParams.ksIMap[0].v3Yuv2RgbOff);
            GetYuvRgbOffFc(CRngFull, 8, &dEITPParams.ksIMap[1].v3Yuv2RgbOff);



            dEITPParams.pxlNumR = w * h;
            dEITPParams.pxlNumR = w * h;
            dEITPParams.pxlNumR = 1 / dEITPParams.pxlNumR;

            float meanDe[2];
            float maxDe[2];
            float sdDe[2];

                //offset_image(ref_buf, -OPT_RANGE_PIXEL_OFFSET, w, h, stride);
                offset_image(dis_buf, -OPT_RANGE_PIXEL_OFFSET, w, h, stride);

                ret = DolbyDe(hKs, meanDe, maxDe, sdDe, ref_buf, temp_buf);
                insert_array_at(thread_data->dlb_mean_deitp_array, meanDe[0], frm_idx);
                insert_array_at(thread_data->dlb_max_deitp_array, maxDe[0], frm_idx);
                insert_array_at(thread_data->dlb_sd_deitp_array, sdDe[0], frm_idx);

                printf("%11s %5.03f %5.03f\n", "MeanDe and MaxDe are: ", meanDe[0], maxDe[0]);
        }

        dbg_printf("\n");

        //Release references to reference and distorted buffers
        release_blur_buf_reference(&thread_data->ref_buf_array, frm_idx);
        release_blur_buf_reference(&thread_data->dis_buf_array, frm_idx);
        release_blur_buf_reference(&thread_data->blur_buf_array, frm_idx);

        if (thread_data->dlb_mean_deitp_array != NULL)
        {
            release_blur_buf_reference(&thread_data->ref_buf_no_conv_array, frm_idx);
            release_blur_buf_reference(&thread_data->ref_buf_cb_array, frm_idx);
            release_blur_buf_reference(&thread_data->dis_buf_cb_array, frm_idx);
            release_blur_buf_reference(&thread_data->ref_buf_cr_array, frm_idx);
            release_blur_buf_reference(&thread_data->dis_buf_cr_array, frm_idx);
        }
        /*Loop through the slots and release slots if there are no more
          reference till the current index. Not releasing next frame as
          it may be required for the next loop						   */
        for(int i = 0; i <= frm_idx; i++)
        {
            int ref_reference_count = get_blur_buf_reference_count(&thread_data->ref_buf_array, i);
            int dis_reference_count = get_blur_buf_reference_count(&thread_data->dis_buf_array, i);

            if ((ref_reference_count == 0) && (dis_reference_count == 0))
            {
                release_blur_buf_slot(&thread_data->ref_buf_array, i);
                release_blur_buf_slot(&thread_data->dis_buf_array, i);
            }

            if (thread_data->dlb_mean_deitp_array != NULL)
            {

                int ref_buf_no_conv_reference_count =
                    get_blur_buf_reference_count(&thread_data->ref_buf_no_conv_array, i);
                int ref_cb_reference_count = get_blur_buf_reference_count(&thread_data->ref_buf_cb_array, i);
                int dis_cb_reference_count = get_blur_buf_reference_count(&thread_data->dis_buf_cb_array, i);

                int ref_cr_reference_count = get_blur_buf_reference_count(&thread_data->ref_buf_cr_array, i);
                int dis_cr_reference_count = get_blur_buf_reference_count(&thread_data->dis_buf_cr_array, i);

                if ((ref_cb_reference_count == 0) && (dis_cb_reference_count == 0) &&
                    (ref_cr_reference_count == 0) && (dis_cr_reference_count == 0) &&
                    (ref_buf_no_conv_reference_count == 0))
                {
                    release_blur_buf_slot(&thread_data->ref_buf_no_conv_array,i);
                    release_blur_buf_slot(&thread_data->ref_buf_cb_array, i);
                    release_blur_buf_slot(&thread_data->dis_buf_cb_array, i);

                    release_blur_buf_slot(&thread_data->ref_buf_cr_array, i);
                    release_blur_buf_slot(&thread_data->dis_buf_cr_array, i);
                }
            }
        }

        /* Loop through the blur buffer array and release slots only till current index - 1 */
        /* Only for those whose reference counter is zero */
        for(int i = 0; i <= (frm_idx - 1); i++)
        {
            int reference_count = get_blur_buf_reference_count(&thread_data->blur_buf_array, i);
            if(reference_count == 0)
            {
                /* Release buffer only if motion score is computed for current, previous and next frame */
                if(
                    (get_at(thread_data->motion_score_compute_flag_array, i)) &&
                    (get_at(thread_data->motion_score_compute_flag_array, i + 1)) &&
                    ((i == 0) || (get_at(thread_data->motion_score_compute_flag_array, i - 1)))
                    )
                {
                    release_blur_buf_slot(&thread_data->blur_buf_array, i);
                }
            }
        }

        /* If this is the last frame then release any subsequent slots */
        if (!next_frame_read)
        {
            release_blur_buf_slot(&thread_data->ref_buf_array, frm_idx + 1);
            release_blur_buf_slot(&thread_data->dis_buf_array, frm_idx + 1);
            release_blur_buf_slot(&thread_data->blur_buf_array, frm_idx);

            if (thread_data->dlb_mean_deitp_array != NULL)
            {
                release_blur_buf_slot(&thread_data->ref_buf_no_conv_array,frm_idx + 1);
                release_blur_buf_slot(&thread_data->ref_buf_cb_array, frm_idx + 1);
                release_blur_buf_slot(&thread_data->dis_buf_cb_array, frm_idx + 1);
                release_blur_buf_slot(&thread_data->ref_buf_cr_array, frm_idx + 1);
                release_blur_buf_slot(&thread_data->dis_buf_cr_array, frm_idx + 1);
            }
        }

        if (!next_frame_read)
        {
            thread_data->stop_threads = 1;
            goto fail_or_end;
        }

    }

fail_or_end:

    aligned_free(temp_buf);

    // when one thread ends we signal all other threads to also stop
    thread_data->stop_threads = 1;
    thread_data->ret = ret;
    pthread_exit(&ret);

}

int combo(int (*read_frame)(float *ref_data, float *main_data, float *temp_data, int stride, void *user_data),
          int (*read_frame_cb_cr)(float *ref_data, float *main_data, float *ref_data_cb, float *main_data_cb, float *ref_data_cr, float *main_data_cr, float *temp_data, int stride, void *user_data),
          void *user_data, int w, int h, const char *fmt,
          DArray *adm_num_array,
          DArray *adm_den_array,
          DArray *adm_num_scale0_array,
          DArray *adm_den_scale0_array,
          DArray *adm_num_scale1_array,
          DArray *adm_den_scale1_array,
          DArray *adm_num_scale2_array,
          DArray *adm_den_scale2_array,
          DArray *adm_num_scale3_array,
          DArray *adm_den_scale3_array,
          DArray *motion_array,
          DArray *motion2_array,
          DArray *vif_num_scale0_array,
          DArray *vif_den_scale0_array,
          DArray *vif_num_scale1_array,
          DArray *vif_den_scale1_array,
          DArray *vif_num_scale2_array,
          DArray *vif_den_scale2_array,
          DArray *vif_num_scale3_array,
          DArray *vif_den_scale3_array,
          DArray *vif_array,
          DArray *psnr_array,
          DArray *ssim_array,
          DArray *ms_ssim_array,
          DArray *dlb_mean_deitp_array,
          DArray *dlb_max_deitp_array,
          DArray *dlb_sd_deitp_array,
          char *errmsg,
          int n_thread,
          int n_subsample)
{
    // init shared thread data
    VMAF_THREAD_STRUCT combo_thread_data;
    combo_thread_data.read_frame = read_frame;
    combo_thread_data.read_frame_cb_cr = read_frame_cb_cr;
    combo_thread_data.user_data = user_data;
    combo_thread_data.w = w;
    combo_thread_data.h = h;
    combo_thread_data.fmt = fmt;
    combo_thread_data.adm_num_array = adm_num_array;
    combo_thread_data.adm_den_array = adm_den_array;
    combo_thread_data.adm_num_scale0_array = adm_num_scale0_array;
    combo_thread_data.adm_den_scale0_array = adm_den_scale0_array;
    combo_thread_data.adm_num_scale1_array = adm_num_scale1_array;
    combo_thread_data.adm_den_scale1_array = adm_den_scale1_array;
    combo_thread_data.adm_num_scale2_array = adm_num_scale2_array;
    combo_thread_data.adm_den_scale2_array = adm_den_scale2_array;
    combo_thread_data.adm_num_scale3_array = adm_num_scale3_array;
    combo_thread_data.adm_den_scale3_array = adm_den_scale3_array;
    combo_thread_data.motion_array = motion_array;
    combo_thread_data.motion2_array = motion2_array;
    combo_thread_data.vif_num_scale0_array = vif_num_scale0_array;
    combo_thread_data.vif_den_scale0_array = vif_den_scale0_array;
    combo_thread_data.vif_num_scale1_array = vif_num_scale1_array;
    combo_thread_data.vif_den_scale1_array = vif_den_scale1_array;
    combo_thread_data.vif_num_scale2_array = vif_num_scale2_array;
    combo_thread_data.vif_den_scale2_array = vif_den_scale2_array;
    combo_thread_data.vif_num_scale3_array = vif_num_scale3_array;
    combo_thread_data.vif_den_scale3_array = vif_den_scale3_array;
    combo_thread_data.vif_array = vif_array;
    combo_thread_data.psnr_array = psnr_array;
    combo_thread_data.ssim_array = ssim_array;
    combo_thread_data.ms_ssim_array = ms_ssim_array;

    combo_thread_data.dlb_mean_deitp_array = dlb_mean_deitp_array;
    combo_thread_data.dlb_max_deitp_array = dlb_max_deitp_array;
    combo_thread_data.dlb_sd_deitp_array = dlb_sd_deitp_array;

    combo_thread_data.errmsg = errmsg;
    combo_thread_data.frm_idx = 0;
    combo_thread_data.stop_threads = 0;
    combo_thread_data.n_subsample = n_subsample;

    DArray	motion_score_compute_flag_array;
    init_array(&motion_score_compute_flag_array, 1000);
    combo_thread_data.motion_score_compute_flag_array = &motion_score_compute_flag_array;

    // sanity check for width/height
    if (w <= 0 || h <= 0 || (size_t)w > ALIGN_FLOOR(INT_MAX) / sizeof(float))
    {
        sprintf(errmsg, "wrong width %d or height %d.\n", w, h);
        return -1;
    }

    // calculate stride and data size
    combo_thread_data.stride = ALIGN_CEIL(w * sizeof(float));
    if ((size_t)h > SIZE_MAX / combo_thread_data.stride)
    {
        sprintf(errmsg, "height %d too large.\n", h);
        return -1;
    }

    if (psnr_constants(fmt, &combo_thread_data.peak, &combo_thread_data.psnr_max))
    {
        sprintf(errmsg, "unknown format %s.\n", fmt);
        return -1;
    }

    combo_thread_data.data_sz = (size_t)combo_thread_data.stride * h;

    if (n_thread == 0)
    {
        combo_thread_data.thread_count = getNumCores();
    }
    else
    {
        combo_thread_data.thread_count = MIN(getNumCores(), n_thread);
    }

    // for motion analysis we compare to previous buffer and next buffer
    /*
     *	In the multi-thread mode, allocate a fixed size buffer pool for the reference, distorted and blur buffers.
     *	At any point, the no. of required ref and dis buffers is 1 more than the total no. of allotted threads,
        to accomodate reading the next frame index.
     *	At any point, one thread operates on the current, previous and next blur buffers, and hence, the no. of
        required blur buffers will be three times the total no. of allotted threads.
     */
    init_blur_array(&combo_thread_data.ref_buf_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
    init_blur_array(&combo_thread_data.dis_buf_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
    init_blur_array(&combo_thread_data.blur_buf_array, MIN(3 * (combo_thread_data.thread_count), MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);

    /*
     * For DEITP and DEISTP calculations, we need to process color buffers too
     * So create buffers reference and distorted color buffers
     * Size of the buffer is same as Luma, as while reading the frame we also to upsampling/coversion to YUV 444
     * Create one extra reference buffer as DEITP does no need convoluted frame. But VMAF seems to do convolution to it.
     * */
    if (dlb_mean_deitp_array != NULL)
    {
        init_blur_array(&combo_thread_data.ref_buf_no_conv_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
        init_blur_array(&combo_thread_data.ref_buf_cb_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
        init_blur_array(&combo_thread_data.dis_buf_cb_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
        init_blur_array(&combo_thread_data.ref_buf_cr_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
        init_blur_array(&combo_thread_data.dis_buf_cr_array, MIN(combo_thread_data.thread_count + 1, MAX_NUM_THREADS), combo_thread_data.data_sz, MAX_ALIGN);
    }
        // initialize the mutex that protects the read_frame function
        pthread_mutex_init(&combo_thread_data.mutex_readframe, NULL);

    // create a joinable thread
    pthread_attr_t attr;
    pthread_attr_init(&attr);
    pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

    // start threads
    int t;
    int numThread = combo_thread_data.thread_count;
    pthread_t* thread = (pthread_t*)calloc(numThread, sizeof(pthread_t));
    memset(thread, 0, numThread * sizeof(pthread_t));

    for (t=0; t < combo_thread_data.thread_count; t++)
    {
        pthread_create(&thread[t], &attr, combo_threadfunc, &combo_thread_data);
    }

    pthread_attr_destroy(&attr);

    // wait for all threads to finish
    for (t=0; t < combo_thread_data.thread_count; t++)
    {
        void* thread_ret;
        int rc = pthread_join(thread[t], &thread_ret);

        if (rc)
        {
            printf("ERROR; return code from pthread_join() for thread[%d] is %d\n", t, (long)thread_ret);
            return -1;
        }
    }

    free_blur_buf(&combo_thread_data.ref_buf_array);
    free_blur_buf(&combo_thread_data.dis_buf_array);
    free_blur_buf(&combo_thread_data.blur_buf_array);

    /* If the color buffers if we have created them */
    if (dlb_mean_deitp_array != NULL)
    {
        free_blur_buf(&combo_thread_data.ref_buf_no_conv_array);
        free_blur_buf(&combo_thread_data.ref_buf_cb_array);
        free_blur_buf(&combo_thread_data.dis_buf_cb_array);
        free_blur_buf(&combo_thread_data.ref_buf_cr_array);
        free_blur_buf(&combo_thread_data.dis_buf_cr_array);
    }
    free_array(&motion_score_compute_flag_array);

    free(thread);

    return 0;
}
