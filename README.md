![Workflow Status](https://github.com/pahdb/AmesPAHdbPythonSuite/actions/workflows/ci.yml/badge.svg) [![Coverage Status]( https://codecov.io/gh/PAHdb/AmesPAHdbPythonSuite/graph/badge.svg)](https://codecov.io/gh/PAHdb/AmesPAHdbPythonSuite) [![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](https://pahdb.github.io/AmesPAHdbPythonSuite/)

# AmesPAHdbPythonSuite

The AmesPAHdbPythonSuite is a package to work with a downloaded PAHdb
XML-file.

A Python module to work with a downloaded PAHdb XML-file.

## Requirements

This software requires:

``python3``
``astropy``
``lxml``
``matplotlib``
``scipy``

## Installation

### Using pip

The AmesPAHdbPythonSuite can be directly installed from its
[repository](https://github.com/PAHdb/AmesPAHdbPythonSuite) using pip:

``pip install git+git://github.com/PAHdb/AmesPAHdbPythonSuite.git``

### From source

Alternatively the AmesPAHdbPythonSuite can be cloned and then installed:

``git clone https://github.com/PAHdb/AmesPAHdbPythonSuite.git``

Then change directories to the new AmesPAHdbPythonSuite directory and install:

``pip install -e .``

## Examples

```python
from pkg_resources import resource_filename
from amespahdbpythonsuite.amespahdb import AmesPAHdb

# Read the database.
xml = 'resources/pahdb-theoretical_cutdown.xml'
pahdb = AmesPAHdb(filename=resource_filename('amespahdbpythonsuite', xml),
                  check=False, cache=False)

# Retrieve the transitions from the database for coronene.
transitions = pahdb.gettransitionsbyuid([18])

# Plot the emission 'stick' spectrum.
transitions.plot(show=True)

# Calculate the emission spectrum at the temperature reached
# after absorbing a 6 eV (CGS units) photon.
transitions.cascade(6 * 1.603e-12, multiprocessing=False)

# Plot the emission 'stick' spectrum at that temperature.
transitions.plot(show=True)

# Convolve the bands with a Gaussian with FWHM of 15 /cm.
convolved = transitions.convolve(fwhm=15.0, gaussian=True,
                                 multiprocessing=False)
convolved.plot(show=True)
```

More examples can be found in the
[examples](examples)-directory.

## Background

The NASA Ames PAH IR Spectroscopic Database and the
AmesPAHdbPythonSuite are being supported through a directed Work
Package at NASA Ames titled: *"Laboratory Astrophysics – The NASA Ames
PAH IR Spectroscopic Database"*.

Additional information can be found at the NASA Ames PAH IR
Spectroscopic Database website, which is located at
[www.astrochemistry.org/pahdb/](https://www.astrochemistry.org/pahdb/).

You are kindly asked to consider the following references for citation
when using the AmesPAHdbPythonSuite:

* C.W. Bauschlicher, Jr., A. Ricca, C. Boersma, and
  L.J. Allamandola, "THE NASA AMES PAH IR SPECTROSCOPIC DATABASE:
  COMPUTATIONAL VERSION 3.00 WITH UPDATED CONTENT AND THE
  INTRODUCTION OF MULTIPLE SCALING FACTORS", The Astrophysical
  Journal Supplement Series, 234, 32, 2018
  [https://doi.org/10.3847/1538-4365/aaa019](https://doi.org/10.3847/1538-4365/aaa019])

* C. Boersma, C.W. Bauschlicher, Jr., A. Ricca, A.L. Mattioda,
  J. Cami, E. Peeters, F. Sanchez de Armas, G. Puerta Saborido,
  D.M. Hudgins, and L.J. Allamandola, "THE NASA AMES PAH IR
  SPECTROSCOPIC DATABASE VERSION 2.00: UPDATED CONTENT, WEBSITE AND
  ON/OFFLINE TOOLS", The Astrophysical Journal Supplement Series,, 211, 8, 2014 [https://doi.org/10.1088/0067-0049/211/1/8](https://doi.org/10.1088/0067-0049/211/1/8)

* Mattioda, A. L., Hudgins, D. M., Boersma, C., Ricca, A.,
  Peeters, E., Cami, J., Sanchez de Armas, F., Puerta Saborido,
  G., Bauschlicher, C. W., J., and Allamandola, L. J. "THE NASA
  AMES PAH IR SPECTROSCOPIC DATABASE: THE LABORATORY SPECTRA", The
  Astrophysical Journal Supplement Series, 251, 22, 2020,
  [https://doi.org/10.3847/1538-4365/abc2c8](https://doi.org/10.3847/1538-4365/abc2c8)

## Contributing

Please read [CONTRIBUTING.md](CONTRIBUTING.md) for details on the code
of conduct, and the process for submitting pull requests.

## Versioning

For the versions available, see the [tags on this
repository](https://github.com/pahdb/amespahdbpythonsuite/tags).

## Authors

* **Christiaan Boersma** - *Initial work* - [PAHdb](https://github.com/pahdb)
* **Matthew J. Shannon** - *Initial work* - [PAHdb](https://github.com/pahdb)
* **Alexandros Maragkoudakis** - *Initial work* - [PAHdb](https://github.com/pahdb)

See also the list of [contributors](CONTRIBUTORS.md) who participated
in this project.

## License

This project is licensed under the BSD 3-Clause License - see the
[LICENSE](LICENSE) file for details

## Acknowledgments

* The NASA Ames PAH IR Spectroscopic Database Team -
  [www.astrochemistry.org/pahdb/](https://www.astrochemistry.org/pahdb/theoretical/3.00/help/about)
* The Astrophysics & Astrochemistry Laboratory at NASA Ames Research
  Center - [www.astrochemistry.org](https://www.astrochemistry.org)
