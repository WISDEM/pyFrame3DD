project('pyframe3dd', 'c',
  version : '1.2',
  meson_version: '>= 0.60',
  default_options: [
    'buildtype=debugoptimized',
  ],
  license : 'Apache')

cc = meson.get_compiler('c')

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

subdir('pyframe3dd')

			    
