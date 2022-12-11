cmake_minimum_required(VERSION 3.11.0)
project(pyframe3dd)
find_package(PythonExtensions REQUIRED)

include_directories(
        ${PROJECT_SOURCE_DIR}/pyframe3dd/src
)

set(my_sources
    py_HPGmatrix.h
    HPGutil.h
    NRutil.h
    coordtrans.h
    py_eig.h
    py_frame3dd.h
    py_io.h
    py_HPGmatrix.c
    HPGutil.c
    NRutil.c
    coordtrans.c
    preframe.c
    py_eig.c
    py_frame3dd.c
    py_io.c
    py_main.c
)
list(TRANSFORM my_sources PREPEND "pyframe3dd/src/")

add_library(_pyframe3dd MODULE ${my_sources})
python_extension_module(_pyframe3dd)
install(TARGETS _pyframe3dd LIBRARY DESTINATION pyframe3dd)

