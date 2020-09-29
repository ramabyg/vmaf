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

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "file_io.h"
#include "read_frame.h"

#include "KdeUpSamp.h"


/**
 * Note: stride is in terms of bytes
 */
static int read_image(FILE *rfile, void *buf, int width, int height, int stride, int elem_size)
{
	char *byte_ptr = buf;
	int i;
	int ret = 1;

	if (width <= 0 || height <= 0 || elem_size <= 0)
	{
		goto fail_or_end;
	}

	for (i = 0; i < height; ++i)
	{
		if (fread(byte_ptr, elem_size, width, rfile) != (size_t)width)
		{
			goto fail_or_end;
		}

		byte_ptr += stride;
	}

	ret = 0;

fail_or_end:
	return ret;
}

/**
 * Note: stride is in terms of bytes
 */
static int read_image_b(FILE * rfile, float *buf, float off, int width, int height, int stride)
{
	char *byte_ptr = (char *)buf;
	unsigned char *tmp_buf = 0;
	int i, j;
	int ret = 1;

	if (width <= 0 || height <= 0)
	{
		goto fail_or_end;
	}

	if (!(tmp_buf = malloc(width)))
	{
		goto fail_or_end;
	}

	for (i = 0; i < height; ++i)
	{
		float *row_ptr = (float *)byte_ptr;

		if (fread(tmp_buf, 1, width, rfile) != (size_t)width)
		{
			goto fail_or_end;
		}

		for (j = 0; j < width; ++j)
		{
			row_ptr[j] = tmp_buf[j] + off;
		}

		byte_ptr += stride;
	}

	ret = 0;

fail_or_end:
	free(tmp_buf);
	return ret;
}

/**
 * Note: stride is in terms of bytes; image is 10-bit little-endian
 */
static int read_image_w(FILE * rfile, float *buf, float off, int width, int height, int stride, float scaler)
{
	// make sure unsigned short is 2 bytes
	assert(sizeof(unsigned short) == 2);

	char *byte_ptr = (char *)buf;
	unsigned short *tmp_buf = 0;
	int i, j;
	int ret = 1;

	if (width <= 0 || height <= 0)
	{
		goto fail_or_end;
	}

	if (!(tmp_buf = malloc(width * 2))) // '*2' to accommodate words
	{
		goto fail_or_end;
	}

	for (i = 0; i < height; ++i)
	{
		float *row_ptr = (float *)byte_ptr;

		if (fread(tmp_buf, 2, width, rfile) != (size_t)width) // '2' for word
		{
			goto fail_or_end;
		}

		for (j = 0; j < width; ++j)
		{
			row_ptr[j] = tmp_buf[j] / scaler + off; // '/4' to convert from x-bit to 8-bit
        }

		byte_ptr += stride;
	}

	ret = 0;

fail_or_end:
	free(tmp_buf);
	return ret;
}

static int completed_frames = 0;

int read_frame(float *ref_data, float *dis_data, float *temp_data, int stride_byte, void *s)
{
    struct data *user_data = (struct data *)s;
    char *fmt = user_data->format;
    int w = user_data->width;
    int h = user_data->height;
    int ret;

    // read ref y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->ref_rfile, ref_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 4.0f);
    }
    else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 16.0f);
    }
    else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 256.0f);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->ref_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // read dis y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->dis_rfile, dis_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 4.0f);
    }
    else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 16.0f);
    }
    else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 256.0f);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->dis_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // ref skip u and v
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        if (fread(temp_data, 1, user_data->offset, user_data->ref_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "ref fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le") ||
             !strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le") ||
             !strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le")
            )
    {
        if (fread(temp_data, 2, user_data->offset, user_data->ref_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "ref fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        goto fail_or_end;
    }

    // dis skip u and v
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        if (fread(temp_data, 1, user_data->offset, user_data->dis_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "dis fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le") ||
             !strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le") ||
             !strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le")
            )
    {
        if (fread(temp_data, 2, user_data->offset, user_data->dis_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "dis fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        goto fail_or_end;
    }

    fprintf(stderr, "Frame: %d/%d\r", completed_frames++, user_data->num_frames);


fail_or_end:
    return ret;
}
// Read frame with color components and if need, do the upsampling to 444
int read_frame_cb_cr(float *ref_data, float *dis_data,
                    float *ref_data_cb, float *dis_data_cb,
                    float *ref_data_cr, float *dis_data_cr,
                    float *temp_data, int stride_byte, void *s)
{
    struct data *user_data = (struct data *)s;
    char *fmt = user_data->format;
    int w = user_data->width;
    int h = user_data->height;
    int ret;

    // read ref y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->ref_rfile, ref_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 4.0f);
    }
    else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 16.0f);
    }
    else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte, 256.0f);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->ref_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // read dis y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->dis_rfile, dis_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 4.0f);
    }
    else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 16.0f);
    }
    else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 256.0f);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->dis_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // ref read cb buffer
    if(!strcmp(fmt, "yuv420p")) {
        /* First read into original buffer, use tmp_buf as intermediate buffer during upsampling*/
        int size = (w >> 1) * (h >> 1);
        //read entire cb plane for the frame
        if (fread(ref_data_cb, 1, size, user_data->ref_rfile) != (size_t)size)
        {
                fprintf(stderr, "ref fread u and v failed.\n");
                goto fail_or_end;
        }
        if (fread(ref_data_cr, 1, size, user_data->ref_rfile) != (size_t)size)
        {
                fprintf(stderr, "ref fread u and v failed.\n");
                goto fail_or_end;
        }

        YUV420ToYUV444InU8(ref_data_cb, ref_data_cr, temp_data, w, h, w, 8);


        //read entire cb plane for the frame
        if (fread(dis_data_cb, 1, size, user_data->dis_rfile) != (size_t)size)
        {
            fprintf(stderr, "ref fread u and v failed.\n");
            goto fail_or_end;
        }
        if (fread(dis_data_cr, 1, size, user_data->dis_rfile) != (size_t)size)
        {
            fprintf(stderr, "ref fread u and v failed.\n");
            goto fail_or_end;
        }

        YUV420ToYUV444InU8(dis_data_cb, dis_data_cr, temp_data, w, h, w, 8);

			// //420 p is assumed to be 8 bit per pixel
			// RowDoublePlnrInU8(ref_data_cb,
	// 	h>>1, w>>1, w>>1, 0, (FloatComp_t)((1 << 8) - 1),
	// 	Tap_6_FilterUvRowUsHalfSize, (FloatComp_t *)Tap_6_FilterUvRowUs0_m, (FloatComp_t *)Tap_6_FilterUvRowUs1_m,
	// 	temp_data, stride_byte);
      //
		// ColDoublePlnr(temp_data,
	// 	h, w>>1, stride_byte, 0, (FloatComp_t)((8) - 1),
	// 	FilterUvColUsHalfSize, (FloatComp_t *)FilterUvColUs_m,
	// 	ref_data_cb, stride_byte);

    } else {
        //at this moment other formats are not supported
        return -1;
    }

    fprintf(stderr, "Frame: %d/%d\r", completed_frames++, user_data->num_frames);

fail_or_end:
    return ret;
}

#if 0
// Read frame with color components and if need, do the upsampling to 444
int read_frame_cb_cr(float *ref_data, float *dis_data,
                    float *ref_data_cb, float *dis_data_cb,
                    float *ref_data_cr, float *dis_data_cr,
                    float *temp_data, int stride_byte, void *s)
{
    struct data *user_data = (struct data *)s;
    char *fmt = user_data->format;
    int w = user_data->width;
    int h = user_data->height;
    int ret;

    // read ref y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->ref_rfile, ref_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->ref_rfile, ref_data, 0, w, h, stride_byte);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->ref_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // read dis y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->dis_rfile, dis_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->dis_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

		// ref read cb buffer
		if(!strcmp(fmt, "yuv420p")) {
			/* First read into original buffer, use tmp_buf as intermediate buffer during upsampling*/
			int size = (w >> 1) * (h >> 1);
			//read entire cb plane for the frame
			if (fread(ref_data_cb, 1, size, user_data->ref_rfile) != (size_t)size)
			{
					fprintf(stderr, "ref fread u and v failed.\n");
					goto fail_or_end;
			}
      if (fread(ref_data_cr, 1, size, user_data->ref_rfile) != (size_t)size)
			{
					fprintf(stderr, "ref fread u and v failed.\n");
					goto fail_or_end;
			}

      YUV420ToYUV444InU8(ref_data_cb, ref_data_cr, temp_data, w, h, w);


      //read entire cb plane for the frame
      if (fread(dis_data_cb, 1, size, user_data->dis_rfile) != (size_t)size)
      {
          fprintf(stderr, "ref fread u and v failed.\n");
          goto fail_or_end;
      }
      if (fread(dis_data_cr, 1, size, user_data->dis_rfile) != (size_t)size)
      {
          fprintf(stderr, "ref fread u and v failed.\n");
          goto fail_or_end;
      }

      YUV420ToYUV444InU8(dis_data_cb, dis_data_cr, temp_data, w, h, w);

			// //420 p is assumed to be 8 bit per pixel
			// RowDoublePlnrInU8(ref_data_cb,
	// 	h>>1, w>>1, w>>1, 0, (FloatComp_t)((1 << 8) - 1),
	// 	Tap_6_FilterUvRowUsHalfSize, (FloatComp_t *)Tap_6_FilterUvRowUs0_m, (FloatComp_t *)Tap_6_FilterUvRowUs1_m,
	// 	temp_data, stride_byte);
      //
		// ColDoublePlnr(temp_data,
	// 	h, w>>1, stride_byte, 0, (FloatComp_t)((8) - 1),
	// 	FilterUvColUsHalfSize, (FloatComp_t *)FilterUvColUs_m,
	// 	ref_data_cb, stride_byte);

		} else {
			//at this moment other formats are not supported
			return -1;
		}

    fprintf(stderr, "Frame: %d/%d\r", completed_frames++, user_data->num_frames);

fail_or_end:
    return ret;
}

#endif

int read_noref_frame(float *dis_data, float *temp_data, int stride_byte, void *s)
{
    struct noref_data *user_data = (struct noref_data *)s;
    char *fmt = user_data->format;
    int w = user_data->width;
    int h = user_data->height;
    int ret;

    // read dis y
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        ret = read_image_b(user_data->dis_rfile, dis_data, 0, w, h, stride_byte);
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 4.0f);
    }
    else if (!strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 16.0f);
    }
    else if (!strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le"))
    {
        ret = read_image_w(user_data->dis_rfile, dis_data, 0, w, h, stride_byte, 256.0f);
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    if (ret)
    {
        if (feof(user_data->dis_rfile))
        {
            ret = 2; // OK if end of file
        }
        return ret;
    }

    // dis skip u and v
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv444p"))
    {
        if (fread(temp_data, 1, user_data->offset, user_data->dis_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "dis fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else if (!strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv444p10le") ||
            !strcmp(fmt, "yuv420p12le") || !strcmp(fmt, "yuv422p12le") || !strcmp(fmt, "yuv444p12le") ||
            !strcmp(fmt, "yuv420p16le") || !strcmp(fmt, "yuv422p16le") || !strcmp(fmt, "yuv444p16le")
             )
    {
        if (fread(temp_data, 2, user_data->offset, user_data->dis_rfile) != (size_t)user_data->offset)
        {
            fprintf(stderr, "dis fread u and v failed.\n");
            goto fail_or_end;
        }
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        goto fail_or_end;
    }


fail_or_end:
    return ret;
}

int get_frame_offset(const char *fmt, int w, int h, size_t *offset)
{
    if (!strcmp(fmt, "yuv420p") || !strcmp(fmt, "yuv420p10le") || !strcmp(fmt, "yuv420p12le"))
    {
        if ((w * h) % 2 != 0)
        {
            fprintf(stderr, "(width * height) %% 2 != 0, width = %d, height = %d.\n", w, h);
            return 1;
        }
        *offset = w * h / 2;
    }
    else if (!strcmp(fmt, "yuv422p") || !strcmp(fmt, "yuv422p10le") || !strcmp(fmt, "yuv422p12le"))
    {
        *offset = w * h;
    }
    else if (!strcmp(fmt, "yuv444p") || !strcmp(fmt, "yuv444p10le") || !strcmp(fmt, "yuv444p12le"))
    {
        *offset = w * h * 2;
    }
    else
    {
        fprintf(stderr, "unknown format %s.\n", fmt);
        return 1;
    }
    return 0;
}
