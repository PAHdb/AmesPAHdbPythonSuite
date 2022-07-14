#!/usr/bin/env python3

from typing import Optional

import numpy as np

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.spectrum import Spectrum

message = AmesPAHdb.message


class Coadded(Spectrum):
    """
    AmesPAHdbPythonSuite coadded class

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.__set(d, **keywords)
        return None

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.Spectrum.spectrum.set` to parse keywords.

        """
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        self.weights = keywords.get("weights", dict())
        self.averaged = keywords.get("averaged", False)

        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "weights" not in keywords:
                    self.weights = d["weights"]
                if "averaged" not in keywords:
                    self.averaged = d["averaged"]

    def get(self) -> dict:
        """
        Calls class: :class:`amespahdbpythonsuite.spectrum.Spectrum.get`.
        Assigns class variables from inherited dictionary.

        """
        d = Spectrum.get(self)
        d["type"] = self.__class__.__name__
        d["weights"] = self.weights
        d["averaged"] = self.averaged

        return d

    def __repr__(self) -> str:
        """
        Class representation.

        """
        return f"{self.__class__.__name__}(" f"{self.uids=})"

    def __str__(self) -> str:
        """
        A description of the instance.
        """

        return f"AmesPAHdbPythonSuite Coadded instance.\n" f"{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the coadded spectrum to file as an IPAC-table.

        """
        import sys
        import datetime
        from astropy.io import ascii  # type: ignore
        from astropy.table import Table  # type: ignore

        if filename == "":
            filename = self.__class__.__name__ + ".tbl"

        hdr = list()

        kv = {
            "DATE": datetime.datetime.now()
            .astimezone()
            .replace(microsecond=0)
            .isoformat(),
            "ORIGIN": "NASA Ames Research Center",
            "CREATOR": f"Python {sys.version_info.major}.{sys.version_info.minor}.{sys.version_info.micro}",
            "SOFTWARE": "AmesPAHdbPythonSuite",
            "AUTHOR": "Dr. C. Boersma",
            "TYPE": self.__class__.__name__.upper(),
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(
            [
                np.array([f for _ in self.data.values() for f in self.grid])
                * self.units["abscissa"]["unit"],
                np.array([t for v in self.data.values() for t in v])
                * self.units["ordinate"]["unit"],
            ],
            names=["FREQUENCY", "INTENSITY"],
            meta={"comments": hdr},
        )

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def plot(self, **keywords) -> None:
        """
        Plot the spectrum.

        """
        import matplotlib.pyplot as plt  # type: ignore
        import matplotlib.cm as cm  # type: ignore

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which="major", right="on", top="on", direction="in")
        colors = cm.rainbow(np.linspace(0, 1, len(self.uids)))
        for y, col in zip(self.data.values(), colors):
            ax.plot(self.grid, y, color=col)

        ax.set_xlim((max(self.grid), min(self.grid)))

        ax.set_xlabel(
            self.units["abscissa"]["label"]
            + " ["
            + self.units["abscissa"]["unit"].to_string("latex_inline")
            + "]",
        )

        unit = self.units["ordinate"]["unit"]
        scale = unit.scale
        unit /= scale
        unit = unit.decompose().cgs.unit
        pre = ""

        if scale != 1.0:
            s = np.log10(scale)
            pre = r"$\times10^{" + f"{s:.0f}" + r"}$ "

        ax.set_ylabel(
            self.units["ordinate"]["label"]
            + " ["
            + pre
            + unit.to_string("latex_inline")
            + "]",
        )
