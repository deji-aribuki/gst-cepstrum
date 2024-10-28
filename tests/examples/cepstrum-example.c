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


#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gst/gst.h>

static guint spect_bands = 20;

#define AUDIOFREQ 16000

/* receive spectral data from element message */
static gboolean
message_handler (GstBus * bus, GstMessage * message, gpointer data)
{
  if (message->type == GST_MESSAGE_ELEMENT) {
    const GstStructure *s = gst_message_get_structure (message);
    const gchar *name = gst_structure_get_name (s);
    GstClockTime endtime;

    if (strcmp (name, "cepstrum") == 0) {
      const GValue *coeffs;
      const GValue *coeff;
      gdouble freq;
      guint i;

      if (!gst_structure_get_clock_time (s, "endtime", &endtime))
        endtime = GST_CLOCK_TIME_NONE;

      g_print ("New cepstrum message, endtime %" GST_TIME_FORMAT "\n",
          GST_TIME_ARGS (endtime));

      coeffs = gst_structure_get_value (s, "coeffs");

      for (i = 0; i < gst_value_list_get_size (coeffs); ++i) {
        coeff = gst_value_list_get_value (coeffs, i);

        if (coeff != NULL) {
          g_print ("band %d mfcc %f\n", i, g_value_get_float (coeff));
        }
      }
      g_print ("\n");
    }
  }
  return TRUE;
}

int
main (int argc, char *argv[])
{
  GstElement *bin;
  GstElement *src, *audioconvert, *cepstrum, *sink;
  GstBus *bus;
  GstCaps *caps;
  GMainLoop *loop;

  gst_init (&argc, &argv);

  bin = gst_pipeline_new ("bin");

  src = gst_element_factory_make ("audiotestsrc", "src");
  g_object_set (G_OBJECT (src), "wave", 0, "freq", 300.0, NULL);
  audioconvert = gst_element_factory_make ("audioconvert", NULL);
  g_assert (audioconvert);

  cepstrum = gst_element_factory_make ("cepstrum", "cepstrum");
  g_object_set (G_OBJECT (cepstrum), "fft-size", 512, "window-size", 512,
      "hop-size", 256, "use-preemphasis", TRUE, "preemphasis-coeff", 0.97,
      "sample-rate", AUDIOFREQ, "num-coeffs", 13, "post-messages", TRUE, NULL);

  sink = gst_element_factory_make ("fakesink", "sink");
  g_object_set (G_OBJECT (sink), "sync", TRUE, NULL);

  gst_bin_add_many (GST_BIN (bin), src, audioconvert, cepstrum, sink, NULL);

  caps = gst_caps_new_simple ("audio/x-raw",
      "rate", G_TYPE_INT, AUDIOFREQ, NULL);

  if (!gst_element_link (src, audioconvert) ||
      !gst_element_link_filtered (audioconvert, cepstrum, caps) ||
      !gst_element_link (cepstrum, sink)) {
    fprintf (stderr, "can't link elements\n");
    exit (1);
  }
  gst_caps_unref (caps);

  bus = gst_element_get_bus (bin);
  gst_bus_add_watch (bus, message_handler, NULL);
  gst_object_unref (bus);

  gst_element_set_state (bin, GST_STATE_PLAYING);

  /* we need to run a GLib main loop to get the messages */
  loop = g_main_loop_new (NULL, FALSE);
  g_main_loop_run (loop);

  gst_element_set_state (bin, GST_STATE_NULL);

  gst_object_unref (bin);

  return 0;
}
