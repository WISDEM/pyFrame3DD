# setup.py
# only if building in place: ``python setup.py build_ext --inplace``
import os
import re
import shutil
import setuptools
import subprocess


def run_meson_build():
    prefix = os.path.join(os.getcwd(), staging_dir)
    purelibdir = "."

    # check if meson extra args are specified
    meson_args = ""
    if "MESON_ARGS" in os.environ:
        meson_args = os.environ["MESON_ARGS"]

    # configure
    meson_path = shutil.which("meson")
    meson_call = (
        f"{meson_path} setup {staging_dir} --prefix={prefix} "
        + f"-Dpython.purelibdir={purelibdir} -Dpython.platlibdir={purelibdir} {meson_args}"
    )
    sysargs = meson_call.split(" ")
    sysargs = [arg for arg in sysargs if arg != ""]
    p1 = subprocess.run(sysargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    setup_log = os.path.join(staging_dir, "setup.log")
    with open(setup_log, "wb") as f:
        f.write(p1.stdout)
    if p1.returncode != 0:
        with open(setup_log, "r") as f:
            print(f.read())
        raise OSError(sysargs, f"The meson setup command failed! Check the log at {setup_log} for more information.")

    # build
    meson_call = f"{meson_path} compile -vC {staging_dir}"
    sysargs = meson_call.split(" ")
    p2 = subprocess.run(sysargs, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    compile_log = os.path.join(staging_dir, "compile.log")
    with open(compile_log, "wb") as f:
        f.write(p2.stdout)
    if p2.returncode != 0:
        with open(compile_log, "r") as f:
            print(f.read())
        raise OSError(
            sysargs, f"The meson compile command failed! Check the log at {compile_log} for more information."
        )


def copy_shared_libraries():
    build_path = os.path.join(staging_dir, "pyframe3dd")
    for root, _dirs, files in os.walk(build_path):
        for file in files:
            # move pyframe3dd to just under staging_dir
            if file.endswith((".so", ".lib", ".pyd", ".pdb", ".dylib", ".dll")):
                if ".so.p" in root or ".pyd.p" in root:  # excludes intermediate object files
                    continue
                file_path = os.path.join(root, file)
                new_path = str(file_path)
                match = re.search(staging_dir, new_path)
                new_path = new_path[match.span()[1] + 1 :]
                print(f"Copying build file {file_path} -> {new_path}")
                shutil.copy(file_path, new_path)


if __name__ == "__main__":
    # This is where the meson build system will install to, it is then
    # used as the sources for setuptools
    staging_dir = "meson_build"

    # this keeps the meson build system from running more than once
    if "dist" not in str(os.path.abspath(__file__)) and not os.path.isdir(staging_dir):
        cwd = os.getcwd()
        run_meson_build()
        os.chdir(cwd)
        copy_shared_libraries()

    #docs_require = ""
    #req_txt = os.path.join("doc", "requirements.txt")
    #if os.path.isfile(req_txt):
    #    with open(req_txt) as f:
    #        docs_require = f.read().splitlines()

    init_file = os.path.join("pyframe3dd", "__init__.py")
    #__version__ = re.findall(
    #    r"""__version__ = ["']+([0-9\.]*)["']+""",
    #    open(init_file).read(),
    #)[0]

    setuptools.setup(
        name="pyFrame3DD",
        version="1.2",
        description="Python bindings to Frame3DD",
        long_description="pyFrame3DD is a Python package that provides a library-like interface to the structural analysis program, Frame3DD",
        author='NREL WISDEM Team',
        author_email='systems.engineering@nrel.gov',
        install_requires=[
            "numpy",
        ],
        extras_require={
            "testing": ["pytest"],
        },
        python_requires=">=3.7",
        packages=["pyframe3dd"],
        license='Apache License, Version 2.0',
        zip_safe=False,
    )

#os.environ['NPY_DISTUTILS_APPEND_FLAGS'] = '1'

#if os.name == 'nt':  # Windows.
#    extra_compile_args = ['/TC', '/D', 'ANSI'] # for msvs
#     # TODO: Not with Anaconda MINGW
#else:
#extra_compile_args = ''

#froot = 'pyframe3dd' + os.sep + 'src' + os.sep
#pyframeExt = Extension('pyframe3dd._pyframe3dd', sources=[froot+'py_HPGmatrix.c',
#                                               froot+'HPGutil.c',
#                                               froot+'NRutil.c',
#                                               froot+'coordtrans.c',
#                                               froot+'preframe.c',
#                                               froot+'py_eig.c',
#                                               froot+'py_frame3dd.c',
#                                               froot+'py_io.c',
#                                               froot+'py_main.c'])
