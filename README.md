[![Travis Status](https://img.shields.io/travis/com/PAHdb/AmesPAHdbPythonSuite)](https://app.travis-ci.com/github/PAHdb/AmesPAHdbPythonSuite) [![Coverage Status]( https://codecov.io/gh/PAHdb/AmesPAHdbPythonSuite/graph/badge.svg)](https://codecov.io/gh/PAHdb/AmesPAHdbPythonSuite) [![Documentation](https://img.shields.io/badge/docs-available-brightgreen.svg)](https://pahdb.github.io/AmesPAHdbPythonSuite/)

# AmesPAHdbPythonSuite

The AmesPAHdbPythonSuite is a package to work with a downloaded PAHdb
XML-file.

A Python module to work with a downloaded PAHdb XML-file.

## Requirements

This software requires:

``python``
``scipy``
``astropy``

## Installation

The AmesPAHdbPythonSuite can be directly installed from the
[repository](https://github.com/PAHdb/AmesPAHdbPythonSuite) using pip:

``pip install git+git://github.com/PAHdb/AmesPAHdbPythonSuite.git``

## Examples

```python
import pkg_resources
from amespahdbpythonsuite.xmlparser import XMLparser
import matplotlib.pyplot as plt

path = 'resources/pahdb-theoretical_cutdown.xml'
xml = pkg_resources.resource_filename('amespahdbpythonsuite', path)
parser = XMLparser(xml)
parser.verify_schema()
library = parser.to_pahdb_dict()
plt.bar([d['frequency'] for d in library['species'][18]['transitions']],
        [d['intensity'] for d in library['species'][18]['transitions']],
        20, color='red', edgecolor="none")
plt.title('stick absorption spectrum of coronene (UID=18)')
plt.xlabel('frequency [cm$^{-1}$]')
plt.ylabel('integrated cross-section [km mol$^{-1}$]')
plt.show()
```

More examples can be found in the
[examples](examples)-directory.

## Documentation

Documentation can be found at
[PAHdb.github.io](https://PAHdb.github.io).

## Background

The NASA Ames PAH IR Spectroscopic Database and the
AmesPAHdbPythonSuite are being supported through a directed Work
Package at NASA Ames titled: *"Laboratory Astrophysics – The NASA Ames
PAH IR Spectroscopic Database"*.

Additional information can be found at the NASA Ames PAH IR
Spectroscopic Database website, which is located at
[www.astrochemistry.org/pahdb](https://www.astrochemistry.org/pahdb/).

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
  [www.astrochemistry.org/pahdb](https://www.astrochemistry.org/pahdb/theoretical/3.00/help/about)
* The Astrophysics & Astrochemistry Laboratory at NASA Ames Research
  Center - [www.astrochemistry.org](https://www.astrochemistry.org)
