project('pyframe3dd', 'c',
  meson_version: '>= 1.1',
  default_options: [
    'buildtype=debugoptimized',
  ],
  license : 'Apache')

cc = meson.get_compiler('c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
m_dep = cc.find_library('m', required : false)
if m_dep.found()
  add_project_link_arguments('-lm', language : 'c')
endif

is_windows = host_machine.system() == 'windows'
add_project_arguments('-DANSI=1', language : 'c')
if is_windows
  add_project_arguments('-D_WIN64=1', language : 'c')
endif

# https://mesonbuild.com/Python-module.html
# Here we differentiate from the python used by meson, py3_command, and that python target, py3_target. This is useful
# when cross compiling like on conda-forge
py3 = import('python').find_installation(pure: false)
py3_dep = py3.dependency()

message(py3.path())
message(py3.get_install_dir())

subdir('pyframe3dd')
