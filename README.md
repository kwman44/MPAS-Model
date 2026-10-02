AMPS Notes
==========

Changes for AMPS:

Call mynn_bl_driver with psfc_hyd_p rather than psfc_p.  Without this
change, our high vertical resolution near the surface seems to trigger
NaNs now and then.

Use a constant snow heat capacity (older value still in the code, but
commented out in the official version) rather than the computed value.
This (combined with using config_noahmp_iopt_tksno=3) seems to
mitigate (though not eliminate) the extreme cold biases we would see
episodically over the Ross Ice Shelf.

Add a wait loop for reading lateral boundary conditions.  This is
vital to our real-time regional run, where the boundary condition
files are created on the fly from a real-time global MPAS forecast.
The wait loop allows some extra time (i.e., in case real-time
processes are delayed) for creating the LBC files, rather than
immediately stopping the forecast if LBC files are not available.

Added AMPS diagnostics.  Currently, only the ceiling product (ported
from the AMPS version of RIP4) is created.

---

For building the scotch library (adapt directories and version numbers
as appropriate):
>```
> export SCOTCH=/glade/work/ampsrt/derecho/src-8.4.2/scotch
>
> git clone https://gitlab.inria.fr/scotch/scotch.git scotch-7.0.16
>
> cd scotch-7.0.16
>
> cmake -DCMAKE_INSTALL_PREFIX=${SCOTCH} \
>    -DBUILD_SHARED_LIBS=ON \
>    -DINTSIZE:STRING=32 \
>    -DIDXSIZE:STRING=32 \
>    -DBISON_EXECUTABLE=/glade/u/apps/derecho/25.10/opt/view/bin/bison \
>    -DFLEX_EXECUTABLE=/glade/u/apps/derecho/25.10/opt/view/bin/flex
>
> make -j 4
>
> make test
>
> make install
>```

MPAS-v8.4.2
====

The Model for Prediction Across Scales (MPAS) is a collaborative project for
developing atmosphere, ocean, and other earth-system simulation components for
use in climate, regional climate, and weather studies. The primary development
partners are the climate modeling group at Los Alamos National Laboratory
(COSIM) and the National Center for Atmospheric Research. Both primary
partners are responsible for the MPAS framework, operators, and tools common to
the applications; LANL has primary responsibility for the ocean model, and NCAR
has primary responsibility for the atmospheric model.

The MPAS framework facilitates the rapid development and prototyping of models
by providing infrastructure typically required by model developers, including
high-level data types, communication routines, and I/O routines. By using MPAS,
developers can leverage pre-existing code and focus more on development of
their model.

BUILDING
========

This README is provided as a brief introduction to the MPAS framework. It does
not provide details about each specific model, nor does it provide building
instructions.

For information about building and running each core, please refer to each
core's user's guide, which can be found at the following web sites:

[MPAS-Atmosphere](http://mpas-dev.github.io/atmosphere/atmosphere_download.html)

[MPAS-Albany Land Ice](http://mpas-dev.github.io/land_ice/download.html)

[MPAS-Ocean](http://mpas-dev.github.io/ocean/releases.html)

[MPAS-Seaice](http://mpas-dev.github.io/sea_ice/releases.html)


Code Layout
----------

Within the MPAS repository, code is laid out as follows. Sub-directories are
only described below the src directory.

	MPAS-Model
	├── src
	│   ├── driver -- Main driver for MPAS in stand-alone mode (Shared)
	│   ├── external -- External software for MPAS (Shared)
	│   ├── framework -- MPAS Framework (Includes DDT Descriptions, and shared routines. Shared)
	│   ├── operators -- MPAS Opeartors (Includes Operators for MPAS meshes. Shared)
	│   ├── tools -- Empty directory for include files that Registry generates (Shared)
	│   │   ├── registry -- Code for building Registry.xml parser (Shared)
	│   │   └── input_gen -- Code for generating streams and namelist files (Shared)
	│   └── core_* -- Individual model cores.
	│       └── inc -- Empty directory for include files that Registry generates
	├── testing_and_setup -- Tools for setting up configurations and test cases (Shared)
	└── default_inputs -- Copies of default stream and namelists files (Shared)

Model cores are typically developed independently. For information about
building and running a particular core, please refer to that core's user's
guide.
