/* GStCepstrum
 * Copyright (C) <2020> Deji Aribuki <daribuki@ketulabs.ch>
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Library General Public
 * License as published by the Free Software Foundation; either
 * version 2 of the License, or (at your option) any later version.
 *
 * This library is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Library General Public License for more details.
 *
 * You should have received a copy of the GNU Library General Public
 * License along with this library; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA 02110-1301, USA.
 */


#ifndef __GST_CEPSTRUM_H__
#define __GST_CEPSTRUM_H__

#include <gst/gst.h>
#include <gst/audio/gstaudiofilter.h>

#ifdef HAVE_LIBFFTW
#include <fftw3.h>
#else
#include <gst/fft/gstfftf32.h>
#endif


G_BEGIN_DECLS

#define GST_TYPE_CEPSTRUM            (gst_cepstrum_get_type())
#define GST_CEPSTRUM(obj)            (G_TYPE_CHECK_INSTANCE_CAST((obj),GST_TYPE_CEPSTRUM,GstCepstrum))
#define GST_IS_CEPSTRUM(obj)         (G_TYPE_CHECK_INSTANCE_TYPE((obj),GST_TYPE_CEPSTRUM))
#define GST_CEPSTRUM_CLASS(klass)    (G_TYPE_CHECK_CLASS_CAST((klass), GST_TYPE_CEPSTRUM,GstCepstrumClass))
#define GST_IS_CEPSTRUM_CLASS(klass) (G_TYPE_CHECK_CLASS_TYPE((klass), GST_TYPE_CEPSTRUM))
typedef struct _GstCepstrum GstCepstrum;
typedef struct _GstCepstrumClass GstCepstrumClass;
typedef struct _GstCepstrumChannel GstCepstrumChannel;

typedef void (*GstCepstrumInputData)(const guint8 * in, gfloat * out,
    guint len, guint channels, gfloat max_value, guint op, guint nfft);

struct _GstCepstrumChannel
{
  gfloat *input;
  gfloat *input_tmp;
#ifdef HAVE_LIBFFTW
  fftw_complex *fftdata;
  fftw_plan fftplan;
#else
  GstFFTF32Complex *fftdata;
  GstFFTF32 *fft_ctx;
#endif
  gfloat *spect_magnitude;
  gfloat *mfcc;
};

struct _GstCepstrum
{
  GstAudioFilter parent;

  /* properties */
  gint num_coeffs;              /* number of mfcc coefficients */
  gint num_filters;             /* number of Mel filter banks */

  gint sample_rate;             /* sampling rate of the audio signal */
  gint fft_size;
  gint win_size;                /* hamming filter window size */
  gint hop_size;                /* hop size */
  gboolean use_preemphasis;     /* whether or not to use preemphasis filter */
  float preemphasis_coeff;      /* filter coefficient */

  gboolean post_messages;       /* whether or not to post messages */
  guint64 interval;             /* how many nanoseconds between emits */
  guint64 frames_per_interval;  /* how many frames per interval */
  guint64 frames_todo;
  gint threshold;               /* energy level threshold */
  gboolean multi_channel;       /* send separate channel results */

  guint64 num_frames;           /* frame count (1 sample per channel)
                                 * since last emit */
  guint64 num_fft;              /* number of FFTs since last emit */
  GstClockTime message_ts;      /* starttime for next message */

  /* <private> */
  GstCepstrumChannel *channel_data;
  guint num_channels;

  guint input_pos;
  guint64 error_per_interval;
  guint64 accumulated_error;

  gfloat **filter_bank;

  GMutex lock;

  GstCepstrumInputData input_data;
};

struct _GstCepstrumClass
{
  GstAudioFilterClass parent_class;
};

GType gst_cepstrum_get_type (void);

GST_ELEMENT_REGISTER_DECLARE (cepstrum);

G_END_DECLS

#endif /* __GST_CEPSTRUM_H__ */
