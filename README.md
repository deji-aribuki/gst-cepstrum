
# GStreamer Cepstrum Plugin

This project implements a GStreamer plugin for calculating the **cepstrum** and **MFCC** (Mel Frequency Cepstral Coefficients) from audio input. The plugin is designed to integrate seamlessly with GStreamer pipelines and allows for real-time audio feature extraction.

## Features
- **MFCC computation**: Computes MFCCs using **libFFTW** for FFT calculations and Mel filterbank operations.
- **Real-time processing**: The plugin processes audio in real-time and outputs the calculated coefficients.
- **Flexible configuration**: Parameters such as FFT size, window size, hop size, and sample rate can be configured.

## Requirements
- **GStreamer**: Version 1.0 or later.
- **libFFTW**: The Fastest Fourier Transform in the West (FFTW) library for efficient FFT computations.
- **GLib**: Required by GStreamer.

Install the dependencies on a system like Ubuntu:
```bash
sudo apt-get install gstreamer1.0 gstreamer1.0-plugins-base libfftw3-dev libglib2.0-dev meson
```

## Building the Plugin

The project uses the **Meson** build system. To build the plugin, follow these steps:

1. Clone the repository:
   ```bash
   git clone https://github.com/deji-aribuki/gst-cepstrum.git
   cd repo
   ```

2. Configure the build environment and compile the plugin:
   ```bash
   meson builddir
   ninja -C builddir
   ```

3. Install the plugin:
   ```bash
   sudo ninja -C builddir install
   ```

## Usage

The plugin provides an element named `cepstrum` which can be used within GStreamer pipelines.

### Example Pipeline

The following example computes MFCCs from a WAV file and outputs the result:

```bash
gst-launch-1.0 filesrc location=audio.wav ! decodebin ! audioconvert ! cepstrum ! fakesink
```

### Configurable Parameters

- **FFT size**: The size of the FFT window (default: 512).
- **Sample rate**: Sample rate of the audio input (default: 16000 Hz).
- **Hop size**: The number of samples between successive windows (default: 160 samples).
- **Number of filters**: The number of Mel filters in the filterbank (default: 26).
- **Number of MFCCs**: The number of MFCC coefficients to compute (default: 13).

## License

This project uses several open-source components, including GStreamer, libFFTW, and Meson. For more details on the licensing of these components, please refer to the `NOTICE` file.
