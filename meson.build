# Much of this is from SciPy

project(
  'pyoptsparse',
  'c',
# unnecessary metadata commented out until Meson supports PEP517 and installation with pip
#  version: 'x.x.x',
#  license: 'Apache',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
  ],
)

cc = meson.get_compiler('c')

# We need -lm for all C code (assuming it uses math functions, which is safe to
# assume for SciPy). For C++ it isn't needed, because libstdc++/libc++ is
# guaranteed to depend on it. For Fortran code, Meson already adds `-lm`.
#m_dep = cc.find_library('m', required : false)
#if m_dep.found()
#  add_project_link_arguments('-lm', language : 'c')
#endif

# https://mesonbuild.com/Python-module.html
# Here we differentiate from the python used by meson, py3_command, and that python target, py3_target. This is useful
# when cross compiling like on conda-forge
py_mod = import('python')
py3_command = py_mod.find_installation()
if get_option('python_target') != ''
    py3_target = py_mod.find_installation(get_option('python_target'))
else
    py3_target = py3_command
endif
py3_dep = py3_target.dependency()

# NumPy include directory - needed in all submodules
incdir_numpy = get_option('incdir_numpy')
if incdir_numpy == ''
    incdir_numpy = run_command(py3_target,
      [
        '-c',
        'import os; os.chdir(".."); import numpy; print(numpy.get_include())'
      ],
      check: true
    ).stdout().strip()
endif
# this creates a raw string which is useful for Windows use of '\' for paths
incdir_numpy = '''@0@'''.format(incdir_numpy)

subdir('pyframe3dd')
