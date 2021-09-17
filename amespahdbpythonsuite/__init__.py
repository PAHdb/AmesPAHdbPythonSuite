"""The AmesPAHdbPythonSuite is developed by the NASA Ames PAH IR
Spectroscopic Database team, asscociated with the Astrophysics &
Astrochemistry Laboratory at NASA Ames Research Center. The NASA Ames
PAH IR Spectroscopic Database and pypahdb are being supported through
a directed Work Package at NASA Ames titled: "Laboratory Astrophysics
- The NASA Ames PAH IR Spectroscopic Database".

Additional information can be found at the NASA Ames PAH IR
Spectroscopic Database website, which is located at
https://www.astrochemistry.org/pahdb You are kindly asked to consider
the following references for citation when using the
AmesPAHdbPythonSuite:

    * C.W. Bauschlicher, Jr., A. Ricca, C. Boersma, and
      L.J. Allamandola, "THE NASA AMES PAH IR SPECTROSCOPIC DATABASE:
      COMPUTATIONAL VERSION 3.00 WITH UPDATED CONTENT AND THE
      INTRODUCTION OF MULTIPLE SCALING FACTORS", The Astrophysical
      Journal Supplement Series, 234, 32, 2018
      https://doi.org/10.3847/1538-4365/aaa019

    * C. Boersma, C.W. Bauschlicher, Jr., A. Ricca, A.L. Mattioda,
      J. Cami, E. Peeters, F. Sanchez de Armas, G. Puerta Saborido,
      D.M. Hudgins, and L.J. Allamandola, "THE NASA AMES PAH IR
      SPECTROSCOPIC DATABASE VERSION 2.00: UPDATED CONTENT, WEBSITE
      AND ON/OFFLINE TOOLS", The Astrophysical Journal Supplement
      Series, 211, 8, 2014 https://doi.org/10.1088/0067-0049/211/1/8

    * Mattioda, A. L., Hudgins, D. M., Boersma, C., Ricca, A.,
      Peeters, E., Cami, J., Sanchez de Armas, F., Puerta Saborido,
      G., Bauschlicher, C. W., J., and Allamandola, L. J. "THE NASA
      AMES PAH IR SPECTROSCOPIC DATABASE: THE LABORATORY SPECTRA", The
      Astrophysical Journal Supplement Series, 251, 22, 2020,
      https://doi.org/10.3847/1538-4365/abc2c8

"""

from ._version import get_versions
__version__ = get_versions()['version']
del get_versions
