# pyFrame3DD

Python bindings to [Frame3DD](http://frame3dd.sourceforge.net)

## Documentation

Browse the [documentation for Frame3DD](http://svn.code.sourceforge.net/p/frame3dd/code/trunk/doc/Frame3DD-manual.html).  This is a Python wrapper to that code with the following modifications:

* Elimination of all input/output files in favor of direct variable passing
* Arbitrary stiffness values can be passed in (rather than only rigid or free).  ``ReactionData(node, Rx, Ry, Rz, Rxx, Ryy, Rzz, rigid=1)`` takes as input the optional parameter rigid (defaults to 1, which is what Frame3DD uses), which defines what number in the reaction inputs corresponds to a rigid connection.  If a user wants to input spring constants in Rx, Ry, etc. those will be used directly in the stiffness matrix.  The parameter ``rigid`` can then be set to anything else like ``-1``.
* Frame3DD allows inclusion of concentrated masses but they only affect the modal analysis.  In pyFrame3DD they also affect the loads.

There is example code that shows usage contained in ``examples/exB.py``.  This follows example (B) Pyramid Frame contained on the [Frame3DD home page](http://frame3dd.sourceforge.net).

## Prerequisites

pyFrame3DD requires a C compiler

## Install (as a library)

For detailed installation instructions of WISDEM modules see <https://github.com/WISDEM/WISDEM> or to install pyFrame3DD by itself do:

    $ pip install WISDEM-pyFrame3DD


## Install (from source)

If you would like to build the project locally from source for easier access to the underlying methods and tests, do:

    $ git clone https://github.com/WISDEM/pyFrame3DD.git
    $ cd pyFrame3DD
    $ pip install .

If developer/editable mode, do the same `git clone` step, but on install do:

    $ pip install --no-build-isolation -e .

The `--no-build-isolation` option is important per [Meson guidelines](https://meson-python.readthedocs.io/en/latest/how-to-guides/editable-installs.html).


## Unit Tests

    $ pytest test

For software issues please use <https://github.com/WISDEM/pyFrame3DD/issues>.  For functionality and theory related questions and comments please use the NWTC forum for [Systems Engineering Software Questions](https://wind.nrel.gov/forum/wind/viewtopic.php?f=34&t=1002).


## License

Frame3DD uses the GNU GPL so this code must also be under the same license.  The larger WISDEM code has a special dispensation to use the Apache License.

