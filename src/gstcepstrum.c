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

/**
 * SECTION:element-cepstrum
 * @title: cepstrum
 *
 * The Cepstrum element computes the Mel-Frequency Cepstrum Coeff. of an audio signal.
 * If the #GstCepstrum:post-messages property is %TRUE, it sends analysis results
 * as element messages named
 * `cepstrum` after each interval of time given
 * by the #GstCepstrum:interval property.
 *
 * The message's structure contains some combination of these fields:
 *
 * * #GstClockTime `timestamp`: the timestamp of the buffer that triggered the message.
 * * #GstClockTime `stream-time`: the stream time of the buffer.
 * * #GstClockTime `running-time`: the running_time of the buffer.
 * * #GstClockTime `duration`: the duration of the buffer.
 * * #GstClockTime `endtime`: the end time of the buffer that triggered the message as stream time (this
 *   is deprecated, as it can be calculated from stream-time + duration)
 * * A #GST_TYPE_LIST value of #gfloat `magnitude`: the level for each frequency band in dB.
 *   All values below the value of the
 *   #GstSpectrum:threshold property will be set to the threshold. Only present
 *   if the #GstSpectrum:message-magnitude property is %TRUE.
 * * A #GST_TYPE_LIST of #gfloat `phase`: The phase for each frequency band. The value is between -pi and pi. Only
 *   present if the #GstSpectrum:message-phase property is %TRUE.
 *
 * If #GstCepstrum:multi-channel property is set to true. magnitude and phase
 * fields will be each a nested #GST_TYPE_ARRAY value. The first dimension are the
 * channels and the second dimension are the values.
 *
 * ## Example application
 *
 * {{ tests/examples/cepstrum/cepstrum-example.c }}
 *
 */

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <fftw3.h>
#include "gstcepstrum.h"

GST_DEBUG_CATEGORY_STATIC (gst_cepstrum_debug);
#define GST_CAT_DEFAULT gst_cepstrum_debug

/* elementfactory information */
#if G_BYTE_ORDER == G_LITTLE_ENDIAN
# define FORMATS "{ S16LE, S24LE, S32LE, F32LE, F64LE }"
#else
# define FORMATS "{ S16BE, S24BE, S32BE, F32BE, F64BE }"
#endif

#define ALLOWED_CAPS \
  GST_AUDIO_CAPS_MAKE (FORMATS) ", " \
  "layout = (string) interleaved"

/* properties */
#define DEFAULT_POST_MESSAGES	    TRUE
#define DEFAULT_MULTI_CHANNEL     FALSE
#define DEFAULT_INTERVAL		      (GST_SECOND / 10)
#define DEFAULT_NUM_COEFFS        13
#define DEFAULT_SAMPLE_RATE       16000
#define DEFAULT_FFT_SIZE          512
#define DEFAULT_WINDOW_SIZE       512
#define DEFAULT_HOP_SIZE          256
#define DEFAULT_USE_PREEMPHASIS   TRUE
#define DEFAULT_PREEMPHASIS_COEFF 0.97 


enum
{
  PROP_0,
  PROP_POST_MESSAGES,
  PROP_INTERVAL,
  PROP_NUM_CEPSTRAL_COEFFS,
  PROP_SAMPLE_RATE,
  PROP_FFT_SIZE,
  PROP_WINDOW_SIZE,
  PROP_HOP_SIZE,
  PROP_USE_PREEMPHASIS,
  PROP_PREEMPHASIS_COEFF,
  PROP_MULTI_CHANNEL
};

#define gst_cepstrum_parent_class parent_class
G_DEFINE_TYPE (GstCepstrum, gst_cepstrum, GST_TYPE_AUDIO_FILTER);
GST_ELEMENT_REGISTER_DEFINE (cepstrum, "cepstrum", GST_RANK_NONE,
    GST_TYPE_CEPSTRUM);

static void gst_cepstrum_finalize (GObject * object);
static void gst_cepstrum_set_property (GObject * object, guint prop_id,
    const GValue * value, GParamSpec * pspec);
static void gst_cepstrum_get_property (GObject * object, guint prop_id,
    GValue * value, GParamSpec * pspec);
static gboolean gst_cepstrum_start (GstBaseTransform * trans);
static gboolean gst_cepstrum_stop (GstBaseTransform * trans);
static GstFlowReturn gst_cepstrum_transform_ip (GstBaseTransform * trans,
    GstBuffer * in);
static gboolean gst_cepstrum_setup (GstAudioFilter * base,
    const GstAudioInfo * info);
static void alloc_mel_filterbank (gfloat **fbank, gint nfilts,
          gint sample_rate, gint nfft);
static void free_mel_filterbank (gfloat **fbank, gint nfilts);

#define PI 3.14159265359
#define HZ_TO_MEL (hz) (2595 * log10 (1 + hz / 700))
#define MEL_TO_HZ (mel) (pow (10, mel / 2595) - 1)

static void
gst_cepstrum_class_init (GstCepstrumClass * klass)
{
  GObjectClass *gobject_class = G_OBJECT_CLASS (klass);
  GstElementClass *element_class = GST_ELEMENT_CLASS (klass);
  GstBaseTransformClass *trans_class = GST_BASE_TRANSFORM_CLASS (klass);
  GstAudioFilterClass *filter_class = GST_AUDIO_FILTER_CLASS (klass);
  GstCaps *caps;

  gobject_class->set_property = gst_cepstrum_set_property;
  gobject_class->get_property = gst_cepstrum_get_property;
  gobject_class->finalize = gst_cepstrum_finalize;

  trans_class->start = GST_DEBUG_FUNCPTR (gst_cepstrum_start);
  trans_class->stop = GST_DEBUG_FUNCPTR (gst_cepstrum_stop);
  trans_class->transform_ip = GST_DEBUG_FUNCPTR (gst_cepstrum_transform_ip);
  trans_class->passthrough_on_same_caps = TRUE;

  filter_class->setup = GST_DEBUG_FUNCPTR (gst_cepstrum_setup);

  g_object_class_install_property (gobject_class, PROP_POST_MESSAGES,
      g_param_spec_boolean ("post-messages", "Post Messages",
          "Whether to post a 'cepstrum' element message on the bus for each "
          "passed interval", DEFAULT_POST_MESSAGES,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_MULTI_CHANNEL,
      g_param_spec_boolean ("multi-channel", "Multichannel results",
          "Send separate results for each channel",
          DEFAULT_MULTI_CHANNEL, G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_INTERVAL,
      g_param_spec_uint64 ("interval", "Interval",
          "Interval of time between message posts (in nanoseconds)",
          1, G_MAXUINT64, DEFAULT_INTERVAL,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_NUM_CEPSTRAL_COEFFS,
      g_param_spec_uint ("num-coeffs", "Number of MFCC coefficients",
          "Number of MFCC coefficients to compute",
          1, 512, DEFAULT_NUM_COEFFS,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_SAMPLE_RATE,
      g_param_spec_int ("sample-rate", "Sample rate",
          "Audio sample rate",
          0, 92000, DEFAULT_SAMPLE_RATE,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_FFT_SIZE,
      g_param_spec_int ("fft-size", "FFT size",
          "FFT size for MFCC computation",
          0, 4096, DEFAULT_FFT_SIZE,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_WINDOW_SIZE,
      g_param_spec_int ("window-size", "Window size",
          "Window size for MFCC computation",
          0, 4096, DEFAULT_WINDOW_SIZE,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_HOP_SIZE,
      g_param_spec_int ("hop-size", "Hop size",
          "Hop size for MFCC computation",
          0, 4096, DEFAULT_HOP_SIZE,
          G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property (gobject_class, PROP_USE_PREEMPHASIS,
      g_param_spec_boolean ("use-preemphasis", "Use Pre-emphasis",
          "Whether to apply pre-emphasis filter for MFCC computation",
          DEFAULT_USE_PREEMPHASIS, G_PARAM_READWRITE | G_PARAM_STATIC_STRINGS));

  g_object_class_install_property(gobject_class, PROP_PREEMPHASIS_COEFF,
      g_param_spec_float("preemphasis-coeff", "Pre-emphasis Coefficient",
      "Coefficient for the pre-emphasis filter",
      0.0, 1.0, DEFAULT_PREEMPHASIS_COEFF, G_PARAM_READWRITE));

  GST_DEBUG_CATEGORY_INIT (gst_cepstrum_debug, "cepstrum", 0,
      "audio cepstrum analyser element");

  gst_element_class_set_static_metadata (element_class, "Cepstrum analyzer",
      "Filter/Analyzer/Audio",
      "Run MFCC on the audio signal, output cepstrum data",
      "Deji Aribuki <deji.aribuki@ketulabs.ch>, <deji.aribuki@gmail.com>");

  caps = gst_caps_from_string (ALLOWED_CAPS);
  gst_audio_filter_class_add_pad_templates (filter_class, caps);
  gst_caps_unref (caps);
}

static void
gst_cepstrum_init (GstCepstrum * cepstrum)
{
  cepstrum->post_messages = DEFAULT_POST_MESSAGES;
  cepstrum->multi_channel = DEFAULT_MULTI_CHANNEL;
  cepstrum->interval = DEFAULT_INTERVAL;
  cepstrum->num_coeffs = DEFAULT_NUM_COEFFS;
  cepstrum->sample_rate = DEFAULT_SAMPLE_RATE;
  cepstrum->fft_size = DEFAULT_FFT_SIZE;
  cepstrum->win_size = DEFAULT_WINDOW_SIZE;
  cepstrum->hop_size = DEFAULT_HOP_SIZE;
  cepstrum->use_preemphasis = DEFAULT_USE_PREEMPHASIS;
  cepstrum->preemphasis_coeff = DEFAULT_PREEMPHASIS_COEFF;

  g_mutex_init (&cepstrum->lock);
}

static void
gst_cepstrum_alloc_channel_data (GstCepstrum * cepstrum)
{
  gint i;
  GstCepstrumChannel *cd;
  guint fft_size = cepstrum->fft_size;
  guint num_coeffs = cepstrum->num_coeffs;
  guint nfilts = cepstrum->num_filters;
  guint nfft = 2 * fft_size - 2;
  guint sample_rate = cepstrum->sample_rate;

  g_assert (cepstrum->channel_data == NULL);

  cepstrum->num_channels = (cepstrum->multi_channel) ?
      GST_AUDIO_FILTER_CHANNELS (cepstrum) : 1;

  GST_DEBUG_OBJECT (cepstrum, "allocating data for %d channels",
      cepstrum->num_channels);

  cepstrum->channel_data = g_new (GstCepstrumChannel, cepstrum->num_channels);
  cepstrum->filter_bank = g_malloc0 (sizeof (gfloat*) * nfilts);
  for (i = 0; i < cepstrum->num_channels; i++) {
    cd = &cepstrum->channel_data[i];
    cd->input = g_new0 (gfloat, nfft);
    cd->input_tmp = g_new0 (gfloat, nfft);
#ifdef HAVE_LIBFFTW
    cd->fftdata = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * nfft);
    cd->fftplan = fftw_plan_dft_r2c_1d(nfft, (double *) cd->input_tmp,
                        cd->fftdata, FFTW_ESTIMATE);
#else
    cd->fftdata = (fftw_complex*) fftw_malloc (sizeof(fftw_complex) * fft_size);
#endif
    cd->spect_magnitude = g_new0 (gfloat, fft_size);
    cd->mfcc = g_new0 (gfloat, nfilts);
  }
  alloc_mel_filterbank (cepstrum->filter_bank, nfilts, sample_rate, nfft);
}

static void
gst_cepstrum_free_channel_data (GstCepstrum * cepstrum)
{
  if (cepstrum->channel_data) {
    gint i;
    GstCepstrumChannel *cd;

    GST_DEBUG_OBJECT (cepstrum, "freeing data for %d channels",
        cepstrum->num_channels);

    for (i = 0; i < cepstrum->num_channels; i++) {
      cd = &cepstrum->channel_data[i];
  #ifdef HAVE_LIBFFTW
      if (cd->fftplan)
        fftw_destroy_plan(cd->fftplan);
      if (cd->fftdata)
        fftw_free(cd->fftdata);
  #else
      if (cd->fft_ctx)
        gst_fft_f32_free (cd->fft_ctx);
      g_free (cd->fftdata);
  #endif
      g_free (cd->input);
      g_free (cd->input_tmp);
      g_free (cd->mfcc);
      g_free (cd->spect_magnitude);
    }
    free_mel_filterbank (cepstrum->filter_bank, cepstrum->num_filters);
    g_free (cepstrum->filter_bank);
    g_free (cepstrum->channel_data);
    cepstrum->channel_data = NULL;
  }
}

static void
gst_cepstrum_flush (GstCepstrum * cepstrum)
{
  cepstrum->num_frames = 0;
  cepstrum->num_fft = 0;

  cepstrum->accumulated_error = 0;
}

static void
gst_cepstrum_reset_state (GstCepstrum * cepstrum)
{
  GST_DEBUG_OBJECT (cepstrum, "resetting state");

  gst_cepstrum_free_channel_data (cepstrum);
  gst_cepstrum_flush (cepstrum);
}

static void
gst_cepstrum_finalize (GObject * object)
{
  GstCepstrum *cepstrum = GST_CEPSTRUM (object);

  gst_cepstrum_reset_state (cepstrum);
  g_mutex_clear (&cepstrum->lock);

  G_OBJECT_CLASS (parent_class)->finalize (object);
}

static void
gst_cepstrum_set_property (GObject * object, guint prop_id,
    const GValue * value, GParamSpec * pspec)
{
  GstCepstrum *filter = GST_CEPSTRUM (object);

  switch (prop_id) {
    case PROP_POST_MESSAGES:
      filter->post_messages = g_value_get_boolean (value);
      break;
    case PROP_INTERVAL:{
      guint64 interval = g_value_get_uint64 (value);
      g_mutex_lock (&filter->lock);
      if (filter->interval != interval) {
        filter->interval = interval;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    case PROP_NUM_CEPSTRAL_COEFFS:{
      guint num_coeffs = g_value_get_uint (value);
      g_mutex_lock (&filter->lock);
      if (filter->num_coeffs != num_coeffs) {
        filter->num_coeffs = num_coeffs;
        filter->num_filters = 2 * num_coeffs;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    case PROP_SAMPLE_RATE:
      filter->sample_rate = g_value_get_int (value);
      break;
    case PROP_FFT_SIZE:{
      guint fft_size = g_value_get_uint (value);
      g_mutex_lock (&filter->lock);
      if (filter->fft_size != fft_size) {
        filter->fft_size = fft_size;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    case PROP_WINDOW_SIZE:{
      guint win_size = g_value_get_uint (value);
      g_mutex_lock (&filter->lock);
      if (filter->win_size != win_size) {
        filter->win_size = win_size;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    case PROP_HOP_SIZE:{
      guint hop_size = g_value_get_uint (value);
      g_mutex_lock (&filter->lock);
      if (filter->hop_size != hop_size) {
        filter->hop_size = hop_size;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    case PROP_USE_PREEMPHASIS:
      filter->use_preemphasis = g_value_get_boolean (value);
      break;
    case PROP_PREEMPHASIS_COEFF:
      filter->preemphasis_coeff = g_value_get_float (value);
      break;
    case PROP_MULTI_CHANNEL:{
      gboolean multi_channel = g_value_get_boolean (value);
      g_mutex_lock (&filter->lock);
      if (filter->multi_channel != multi_channel) {
        filter->multi_channel = multi_channel;
        gst_cepstrum_reset_state (filter);
      }
      g_mutex_unlock (&filter->lock);
      break;
    }
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static void
gst_cepstrum_get_property (GObject * object, guint prop_id,
    GValue * value, GParamSpec * pspec)
{
  GstCepstrum *filter = GST_CEPSTRUM (object);

  switch (prop_id) {
    case PROP_POST_MESSAGES:
      g_value_set_boolean (value, filter->post_messages);
      break;
    case PROP_INTERVAL:
      g_value_set_uint64 (value, filter->interval);
      break;
    case PROP_NUM_CEPSTRAL_COEFFS:
      g_value_set_uint (value, filter->num_coeffs);
      break;
    case PROP_SAMPLE_RATE:
      g_value_set_int (value, filter->sample_rate);
      break;
    case PROP_FFT_SIZE:
      g_value_set_int (value, filter->fft_size);
      break;
    case PROP_WINDOW_SIZE:
      g_value_set_int (value, filter->win_size);
      break;
    case PROP_HOP_SIZE:
      g_value_set_int (value, filter->hop_size);
      break;
    case PROP_USE_PREEMPHASIS:
      g_value_set_boolean (value, filter->use_preemphasis);
      break;
    case PROP_PREEMPHASIS_COEFF:
      g_value_set_float (value, filter->preemphasis_coeff);
      break;
    case PROP_MULTI_CHANNEL:
      g_value_set_boolean (value, filter->multi_channel);
      break;
    default:
      G_OBJECT_WARN_INVALID_PROPERTY_ID (object, prop_id, pspec);
      break;
  }
}

static gboolean
gst_cepstrum_start (GstBaseTransform * trans)
{
  GstCepstrum *cepstrum = GST_CEPSTRUM (trans);

  gst_cepstrum_reset_state (cepstrum);

  return TRUE;
}

static gboolean
gst_cepstrum_stop (GstBaseTransform * trans)
{
  GstCepstrum *cepstrum = GST_CEPSTRUM (trans);

  gst_cepstrum_reset_state (cepstrum);

  return TRUE;
}

/* mixing data readers */

static void
input_data_mixed_float (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint i, j, ip = 0;
  gfloat v;
  gfloat *in = (gfloat *) _in;

  for (j = 0; j < len; j++) {
    v = in[ip++];
    for (i = 1; i < channels; i++)
      v += in[ip++];
    out[op] = v / channels;
    op = (op + 1) % nfft;
  }
}

static void
input_data_mixed_double (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint i, j, ip = 0;
  gfloat v;
  gdouble *in = (gdouble *) _in;

  for (j = 0; j < len; j++) {
    v = in[ip++];
    for (i = 1; i < channels; i++)
      v += in[ip++];
    out[op] = v / channels;
    op = (op + 1) % nfft;
  }
}

static void
input_data_mixed_int32_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint i, j, ip = 0;
  gint32 *in = (gint32 *) _in;
  gfloat v;

  for (j = 0; j < len; j++) {
    v = in[ip++] / max_value;
    for (i = 1; i < channels; i++)
      v += in[ip++] / max_value;
    out[op] = v / channels;
    op = (op + 1) % nfft;
  }
}

static void
input_data_mixed_int24_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint i, j;
  gfloat v = 0.0;

  for (j = 0; j < len; j++) {
    for (i = 0; i < channels; i++) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
      gint32 value = GST_READ_UINT24_BE (_in);
#else
      gint32 value = GST_READ_UINT24_LE (_in);
#endif
      if (value & 0x00800000)
        value |= 0xff000000;
      v += value / max_value;
      _in += 3;
    }
    out[op] = v / channels;
    op = (op + 1) % nfft;
  }
}

static void
input_data_mixed_int16_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint i, j, ip = 0;
  gint16 *in = (gint16 *) _in;
  gfloat v;

  for (j = 0; j < len; j++) {
    v = in[ip++] / max_value;
    for (i = 1; i < channels; i++)
      v += in[ip++] / max_value;
    out[op] = v / channels;
    op = (op + 1) % nfft;
  }
}

/* non mixing data readers */

static void
input_data_float (const guint8 * _in, gfloat * out, guint len, guint channels,
    gfloat max_value, guint op, guint nfft)
{
  guint j, ip;
  gfloat *in = (gfloat *) _in;

  for (j = 0, ip = 0; j < len; j++, ip += channels) {
    out[op] = in[ip];
    op = (op + 1) % nfft;
  }
}

static void
input_data_double (const guint8 * _in, gfloat * out, guint len, guint channels,
    gfloat max_value, guint op, guint nfft)
{
  guint j, ip;
  gdouble *in = (gdouble *) _in;

  for (j = 0, ip = 0; j < len; j++, ip += channels) {
    out[op] = in[ip];
    op = (op + 1) % nfft;
  }
}

static void
input_data_int32_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint j, ip;
  gint32 *in = (gint32 *) _in;

  for (j = 0, ip = 0; j < len; j++, ip += channels) {
    out[op] = in[ip] / max_value;
    op = (op + 1) % nfft;
  }
}

static void
input_data_int24_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint j;

  for (j = 0; j < len; j++) {
#if G_BYTE_ORDER == G_BIG_ENDIAN
    gint32 v = GST_READ_UINT24_BE (_in);
#else
    gint32 v = GST_READ_UINT24_LE (_in);
#endif
    if (v & 0x00800000)
      v |= 0xff000000;
    _in += 3 * channels;
    out[op] = v / max_value;
    op = (op + 1) % nfft;
  }
}

static void
input_data_int16_max (const guint8 * _in, gfloat * out, guint len,
    guint channels, gfloat max_value, guint op, guint nfft)
{
  guint j, ip;
  gint16 *in = (gint16 *) _in;

  for (j = 0, ip = 0; j < len; j++, ip += channels) {
    out[op] = in[ip] / max_value;
    op = (op + 1) % nfft;
  }
}

static gboolean
gst_cepstrum_setup (GstAudioFilter * base, const GstAudioInfo * info)
{
  GstCepstrum *cepstrum = GST_CEPSTRUM (base);
  gboolean multi_channel = cepstrum->multi_channel;
  GstCepstrumInputData input_data = NULL;

  g_mutex_lock (&cepstrum->lock);
  switch (GST_AUDIO_INFO_FORMAT (info)) {
    case GST_AUDIO_FORMAT_S16:
      input_data =
          multi_channel ? input_data_int16_max : input_data_mixed_int16_max;
      break;
    case GST_AUDIO_FORMAT_S24:
      input_data =
          multi_channel ? input_data_int24_max : input_data_mixed_int24_max;
      break;
    case GST_AUDIO_FORMAT_S32:
      input_data =
          multi_channel ? input_data_int32_max : input_data_mixed_int32_max;
      break;
    case GST_AUDIO_FORMAT_F32:
      input_data = multi_channel ? input_data_float : input_data_mixed_float;
      break;
    case GST_AUDIO_FORMAT_F64:
      input_data = multi_channel ? input_data_double : input_data_mixed_double;
      break;
    default:
      g_assert_not_reached ();
      break;
  }
  cepstrum->input_data = input_data;

  gst_cepstrum_reset_state (cepstrum);
  g_mutex_unlock (&cepstrum->lock);

  return TRUE;
}

static GValue *
gst_cepstrum_message_add_container (GstStructure * s, GType type,
    const gchar * name)
{
  GValue v = { 0, };

  g_value_init (&v, type);
  /* will copy-by-value */
  gst_structure_set_value (s, name, &v);
  g_value_unset (&v);
  return (GValue *) gst_structure_get_value (s, name);
}

static void
gst_cepstrum_message_add_list (GValue * cv, gfloat * data, guint num_values)
{
  GValue v = { 0, };
  guint i;

  g_value_init (&v, G_TYPE_FLOAT);
  for (i = 0; i < num_values; i++) {
    g_value_set_float (&v, data[i]);
    gst_value_list_append_value (cv, &v);       /* copies by value */
  }
  g_value_unset (&v);
}

static void
gst_cepstrum_message_add_array (GValue * cv, gfloat * data, guint num_values)
{
  GValue v = { 0, };
  GValue a = { 0, };
  guint i;

  g_value_init (&a, GST_TYPE_ARRAY);

  g_value_init (&v, G_TYPE_FLOAT);
  for (i = 0; i < num_values; i++) {
    g_value_set_float (&v, data[i]);
    gst_value_array_append_value (&a, &v);      /* copies by value */
  }
  g_value_unset (&v);

  gst_value_array_append_value (cv, &a);        /* copies by value */
  g_value_unset (&a);
}

static GstMessage *
gst_cepstrum_message_new (GstCepstrum * cepstrum, GstClockTime timestamp,
    GstClockTime duration)
{
  GstBaseTransform *trans = GST_BASE_TRANSFORM_CAST (cepstrum);
  GstCepstrumChannel *cd;
  GstStructure *s;
  GValue *mcv = NULL, *pcv = NULL;
  GstClockTime endtime, running_time, stream_time;

  GST_DEBUG_OBJECT (cepstrum,
      "preparing message, coeffs =%d bands =%d",
      cepstrum->num_coeffs, cepstrum->fft_size);

  running_time = gst_segment_to_running_time (&trans->segment, GST_FORMAT_TIME,
      timestamp);
  stream_time = gst_segment_to_stream_time (&trans->segment, GST_FORMAT_TIME,
      timestamp);
  /* endtime is for backwards compatibility */
  endtime = stream_time + duration;

  s = gst_structure_new ("cepstrum",
      "endtime", GST_TYPE_CLOCK_TIME, endtime,
      "timestamp", G_TYPE_UINT64, timestamp,
      "stream-time", G_TYPE_UINT64, stream_time,
      "running-time", G_TYPE_UINT64, running_time,
      "duration", G_TYPE_UINT64, duration, NULL);

  if (!cepstrum->multi_channel) {
    cd = &cepstrum->channel_data[0];

      /* FIXME 0.11: this should be an array, not a list */
      mcv = gst_cepstrum_message_add_container (s, GST_TYPE_LIST, "magnitude");
      gst_cepstrum_message_add_list (mcv, cd->mfcc, cepstrum->num_coeffs);
  } else {
    guint c;
    guint channels = GST_AUDIO_FILTER_CHANNELS (cepstrum);

    mcv = gst_cepstrum_message_add_container (s, GST_TYPE_ARRAY, "magnitude");
    for (c = 0; c < channels; c++) {
      cd = &cepstrum->channel_data[c];
      gst_cepstrum_message_add_array (mcv, cd->mfcc,
          cepstrum->num_coeffs);
    }
  }
  return gst_message_new_element (GST_OBJECT (cepstrum), s);
}

static void
pre_emphasis (gfloat *data, gint size, float alpha)
{
  int i;
  for (i = size - 1; i > 0; i--) {
      data[i] = data[i] - alpha * data[i - 1];
  } 
}

static void
hamming_window (float *data, guint size)
{
  for (guint i = 0; i < size; ++i) {
      data[i] *= (0.54 - 0.46 * cos((2 * PI * i) / (size - 1)));
  }
}

static void
compute_dct (float *in, float *out, guint size)
{
  for (guint k = 0; k < size; k++) {
      out[k] = 0.0;
      for (guint n = 0; n < size; n++) {
          out[k] += in[n] * cos(PI * k * (n + 0.5) / size);
      }
  }
}

static inline gfloat hz_to_mel(gfloat hz)
{
  return 2595.0 * log10f(1.0 + hz / 700.0);
}

static inline gfloat mel_to_hz(gfloat mel)
{
  return 700.0 * (powf(10.0, mel / 2595.0) - 1.0);
}

static void
alloc_mel_filterbank (gfloat **fbank, gint nfilts,
                gint sample_rate, gint nfft)
{
  gfloat *bin = g_malloc0 ((nfilts + 2) * sizeof(gfloat));

  gfloat lowmel = hz_to_mel (0.0);
  gfloat highmel = hz_to_mel (sample_rate / 2.0);
  gfloat mel_step = (highmel - lowmel) / (nfilts + 1);

  /* calculate Mel center frequencies and convert to FFT bin numbers */
  for (guint i = 0; i <= nfilts + 1; i++) {
    bin[i] = mel_to_hz (lowmel + i * mel_step);
    bin[i] = floor ((nfft + 1) * bin[i] / sample_rate);
  }

  /* create triangular filters */
  for (gint i = 1; i <= nfilts; i++) {
    fbank[i] = g_malloc0 (nfft * sizeof (gfloat));
    for (gint k = 0; k < nfft; k++)
      fbank[i][k] = 0.0f;
    for (guint k = bin[i-1]; k < bin[i]; k++)
        fbank[i][k] = (k - bin[i-1]) / (bin[i] - bin[i-1]);
    for (guint k = bin[i]; k < bin[i+1]; k++)
        fbank[i][k] = (bin[i+1] - k) / (bin[i+1] - bin[i]);
  }

  g_free(bin);
}

static void
free_mel_filterbank (gfloat **fbank, gint nfilts)
{
    for (guint i = 0; i < nfilts; i++) {
        g_free(fbank[i]);
    }
}

static void
compute_mel_filterbank (gfloat *in, gfloat *out, gfloat **fbank, guint nfilts,
          guint nfft)
{
    for (guint i = 0; i < nfilts; i++) {
        out[i] = 0.0;
        for (guint j = 0; j < nfft / 2; j++) {
            out[i] += in[j] * fbank[i][j];
        }
        out[i] = log(in[i] + 1e-10);  /* take log for stability */
    }
}

#ifdef HAVE_LIBFFTW
static void
gst_cepstrum_fft (GstCepstrum * cepstrum, GstCepstrumChannel * cd)
{
  guint fft_size = cepstrum->fft_size;
  guint nfft = 2 * fft_size - 2;
  gfloat *spect_magnitude = cd->spect_magnitude;
  gdouble val;
  fftw_complex *fftdata = cd->fftdata;
  fftw_plan fftplan = cd->fftplan;

  fftw_execute (fftplan);

  /* compute power spectrum */
  for (guint i = 0; i < fft_size; i++) {
    val = fftdata[i][0] * fftdata[i][0];
    val += fftdata[i][1] * fftdata[i][1];
    val /= nfft * nfft;
    spect_magnitude[i] += val;
  }
}
#else
static void
gst_cepstrum_fft (GstCepstrum * cepstrum, GstCepstrumChannel * cd)
{
  guint fft_size = cepstrum->fft_size;
  guint nfft = 2 * fft_size - 2;
  gfloat *spect_magnitude = cd->spect_magnitude;
  gdouble val;
  gfloat *input_tmp = cd->input_tmp;
  GstFFTF32Complex *fftdata = cd->fftdata;
  GstFFTF32 *fft_ctx = cd->fft_ctx;

  gst_fft_f32_window (fft_ctx, input_tmp, GST_FFT_WINDOW_HAMMING);
  gst_fft_f32_fft (fft_ctx, input_tmp, freqdata);

  /* compute power spectrum */
  for (guint i = 0; i < fft_size; i++) {
    val = fftdata[i].r * fftdata[i].r;
    val += fftdata[i].i * fftdata[i].i;
    val /= nfft * nfft;
    spect_magnitude[i] += val;
  }
}
#endif

static void
gst_cepstrum_run_mfcc (GstCepstrum *cepstrum, GstCepstrumChannel *cd,
    guint input_pos)
{
  guint i;
  guint frame_size = cepstrum->win_size;
  guint fft_size = cepstrum->fft_size;
  guint nfft = 2 * fft_size - 2;
  guint nfilts = cepstrum->num_filters;
  gfloat *input = cd->input;
  gfloat *input_tmp = cd->input_tmp;
  gfloat *spect_magnitude = cd->spect_magnitude;
  gfloat *mfcc = cd->mfcc;
  gfloat **fbank = cepstrum->filter_bank;
  guint numcoeffs = cepstrum->num_coeffs;
  gfloat alpha = cepstrum->preemphasis_coeff;
  gboolean use_preemphasis = cepstrum->use_preemphasis;

  for (i = 0; i < frame_size; i++)
    input_tmp[i] = input[(input_pos + i) % frame_size];

  if (use_preemphasis)
    pre_emphasis (input_tmp, frame_size, alpha);

  /* apply hamming window to input data */
  hamming_window(input_tmp, frame_size);

  /* run FFT */
  gst_cepstrum_fft (cepstrum, cd);

  /* apply Mel filterbank */
  compute_mel_filterbank (spect_magnitude, mfcc, fbank, nfilts, nfft);

  /* apply DCT to Mel coefficients to get MFCCs */
  compute_dct(mfcc, mfcc, numcoeffs);
}

static void
gst_cepstrum_prepare_message_data (GstCepstrum * cepstrum,
    GstCepstrumChannel * cd)
{
  guint i;
  guint fft_size = cepstrum->fft_size;
  guint num_fft = cepstrum->num_fft;
  gfloat *spect_magnitude = cd->spect_magnitude;

  /* Calculate average */
  for (i = 0; i < fft_size; i++)
    spect_magnitude[i] /= num_fft;
}

static void
gst_cepstrum_reset_message_data (GstCepstrum * cepstrum,
    GstCepstrumChannel * cd)
{
  guint mfcc_size = cepstrum->num_coeffs;
  gfloat *mfcc = cd->mfcc;

  /* reset accumulators */
  memset (mfcc, 0, mfcc_size * sizeof (gfloat));
}

static GstFlowReturn
gst_cepstrum_transform_ip (GstBaseTransform * trans, GstBuffer * buffer)
{
  GstCepstrum *cepstrum = GST_CEPSTRUM (trans);
  guint rate = GST_AUDIO_FILTER_RATE (cepstrum);
  guint channels = GST_AUDIO_FILTER_CHANNELS (cepstrum);
  guint bps = GST_AUDIO_FILTER_BPS (cepstrum);
  guint bpf = GST_AUDIO_FILTER_BPF (cepstrum);
  guint output_channels = cepstrum->multi_channel ? channels : 1;
  guint c;
  gfloat max_value = (1UL << ((bps << 3) - 1)) - 1;
  guint fft_size = cepstrum->fft_size;
  guint nfft = 2 * fft_size - 2;
  guint input_pos;
  gfloat *input;
  GstMapInfo map;
  const guint8 *data;
  gsize size;
  guint fft_todo, msg_todo, block_size;
  gboolean have_full_interval;
  GstCepstrumChannel *cd;
  GstCepstrumInputData input_data;

  g_mutex_lock (&cepstrum->lock);
  gst_buffer_map (buffer, &map, GST_MAP_READ);
  data = map.data;
  size = map.size;

  GST_LOG_OBJECT (cepstrum, "input size: %" G_GSIZE_FORMAT " bytes", size);

  if (GST_BUFFER_IS_DISCONT (buffer)) {
    GST_DEBUG_OBJECT (cepstrum, "Discontinuity detected -- flushing");
    gst_cepstrum_flush (cepstrum);
  }

  /* If we don't have a FFT context yet (or it was reset due to parameter
   * changes) get one and allocate memory for everything
   */
  if (cepstrum->channel_data == NULL) {
    GST_DEBUG_OBJECT (cepstrum, "allocating for bands %u", fft_size);

    gst_cepstrum_alloc_channel_data (cepstrum);

    /* number of sample frames we process before posting a message
     * interval is in ns */
    cepstrum->frames_per_interval =
        gst_util_uint64_scale (cepstrum->interval, rate, GST_SECOND);
    cepstrum->frames_todo = cepstrum->frames_per_interval;
    /* rounding error for frames_per_interval in ns,
     * aggregated it in accumulated_error */
    cepstrum->error_per_interval = (cepstrum->interval * rate) % GST_SECOND;
    if (cepstrum->frames_per_interval == 0)
      cepstrum->frames_per_interval = 1;

    GST_INFO_OBJECT (cepstrum, "interval %" GST_TIME_FORMAT ", fpi %"
        G_GUINT64_FORMAT ", error %" GST_TIME_FORMAT,
        GST_TIME_ARGS (cepstrum->interval), cepstrum->frames_per_interval,
        GST_TIME_ARGS (cepstrum->error_per_interval));

    cepstrum->input_pos = 0;

    gst_cepstrum_flush (cepstrum);
  }

  if (cepstrum->num_frames == 0)
    cepstrum->message_ts = GST_BUFFER_TIMESTAMP (buffer);

  input_pos = cepstrum->input_pos;
  input_data = cepstrum->input_data;

  while (size >= bpf) {
    /* run input_data for a chunk of data */
    fft_todo = nfft - (cepstrum->num_frames % nfft);
    msg_todo = cepstrum->frames_todo - cepstrum->num_frames;
    GST_LOG_OBJECT (cepstrum,
        "message frames todo: %u, fft frames todo: %u, input frames %"
        G_GSIZE_FORMAT, msg_todo, fft_todo, (size / bpf));
    block_size = msg_todo;
    if (block_size > (size / bpf))
      block_size = (size / bpf);
    if (block_size > fft_todo)
      block_size = fft_todo;

    for (c = 0; c < output_channels; c++) {
      cd = &cepstrum->channel_data[c];
      input = cd->input;
      /* Move the current frames into our ringbuffers */
      input_data (data + c * bps, input, block_size, channels, max_value,
          input_pos, nfft);
    }
    data += block_size * bpf;
    size -= block_size * bpf;
    input_pos = (input_pos + block_size) % nfft;
    cepstrum->num_frames += block_size;

    have_full_interval = (cepstrum->num_frames == cepstrum->frames_todo);

    GST_LOG_OBJECT (cepstrum,
        "size: %" G_GSIZE_FORMAT ", do-fft = %d, do-message = %d", size,
        (cepstrum->num_frames % nfft == 0), have_full_interval);

    /* If we have enough frames for an FFT or we have all frames required for
     * the interval and we haven't run a FFT, then run an FFT */
    if ((cepstrum->num_frames % nfft == 0) ||
        (have_full_interval && !cepstrum->num_fft)) {
      for (c = 0; c < output_channels; c++) {
        cd = &cepstrum->channel_data[c];
        gst_cepstrum_run_mfcc (cepstrum, cd, input_pos);
      }
      cepstrum->num_fft++;
    }

    /* Do we have the FFTs for one interval? */
    if (have_full_interval) {
      GST_DEBUG_OBJECT (cepstrum, "nfft: %u frames: %" G_GUINT64_FORMAT
          " fpi: %" G_GUINT64_FORMAT " error: %" GST_TIME_FORMAT, nfft,
          cepstrum->num_frames, cepstrum->frames_per_interval,
          GST_TIME_ARGS (cepstrum->accumulated_error));

      cepstrum->frames_todo = cepstrum->frames_per_interval;
      if (cepstrum->accumulated_error >= GST_SECOND) {
        cepstrum->accumulated_error -= GST_SECOND;
        cepstrum->frames_todo++;
      }
      cepstrum->accumulated_error += cepstrum->error_per_interval;

      if (cepstrum->post_messages) {
        GstMessage *m;

        for (c = 0; c < output_channels; c++) {
          cd = &cepstrum->channel_data[c];
          gst_cepstrum_prepare_message_data (cepstrum, cd);
        }

        m = gst_cepstrum_message_new (cepstrum, cepstrum->message_ts,
            cepstrum->interval);

        gst_element_post_message (GST_ELEMENT (cepstrum), m);
      }

      if (GST_CLOCK_TIME_IS_VALID (cepstrum->message_ts))
        cepstrum->message_ts +=
            gst_util_uint64_scale (cepstrum->num_frames, GST_SECOND, rate);

      for (c = 0; c < output_channels; c++) {
        cd = &cepstrum->channel_data[c];
        gst_cepstrum_reset_message_data (cepstrum, cd);
      }
      cepstrum->num_frames = 0;
      cepstrum->num_fft = 0;
    }
  }

  cepstrum->input_pos = input_pos;

  gst_buffer_unmap (buffer, &map);
  g_mutex_unlock (&cepstrum->lock);

  g_assert (size == 0);

  return GST_FLOW_OK;
}

static gboolean
plugin_init (GstPlugin * plugin)
{

  return GST_ELEMENT_REGISTER (cepstrum, plugin);
}

GST_PLUGIN_DEFINE (GST_VERSION_MAJOR,
    GST_VERSION_MINOR,
    cepstrum,
    "Run MFCC on the audio signal, output cepstrum coefficients",
    plugin_init,
    VERSION, "LGPL", "GStreamer",
    "https://github.com/deji-aribuki/gst-cepstrum");
