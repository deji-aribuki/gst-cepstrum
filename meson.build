project('gst-cepstrum', 'c')

cc = meson.get_compiler('c')

gst_dep = dependency('gstreamer-1.0')
gstaudio_dep = dependency('gstreamer-audio-1.0')
gstfft_dep = dependency('gstreamer-fft-1.0', required: true) 
fftw_dep = dependency('fftw3', required: false)
libm_dep = cc.find_library('m', required: true)

conf = configuration_data()
conf.set('VERSION', '"@0@"'.format(meson.project_version()))
conf.set('PACKAGE', '"gst-cepstrum"')
configure_file(output: 'config.h', configuration: conf)
add_project_arguments('-DHAVE_CONFIG_H', language: 'c')

fftw_cflags = []
if fftw_dep.found()
  fftw_cflags += ['-DHAVE_LIBFFTW']
endif

shared_library('gstcepstrum', 'src/gstcepstrum.c',
  dependencies: [gst_dep, gstaudio_dep, gstfft_dep, fftw_dep, libm_dep],
  include_directories: include_directories('src'),
  c_args : fftw_cflags,
  install: true,
  install_dir: get_option('libdir') / 'gstreamer-1.0'
)
