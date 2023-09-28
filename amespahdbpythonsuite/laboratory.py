#!/usr/bin/env python3

from typing import Optional

import numpy as np

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.data import Data

message = AmesPAHdb.message


class Laboratory(Data):
    """
    AmesPAHdbPythonSuite laboratory class.
    Contains methods to work with a laboratory spectrum.

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.set(d, **keywords)

    def set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.set` to parse keywords.

        """
        Data.set(self, d, **keywords)

    def get(self) -> dict:
        """
        Assigns class variables from inherited dictionary.

        """
        d = Data.get(self)
        d["type"] = self.__class__.__name__

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

        return f"AmesPAHdbPythonSuite Laboratory instance.\n" f"{self.uids=}"

    def write(self, filename: str = "") -> None:
        """
        Write the laboratory spectra to file as an IPAC-table.

        """
        import datetime
        import sys

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
            "SPECIES": str(len(self.data)),
        }

        for key, value in kv.items():
            if not value.isnumeric():
                hdr.append(f"{key:8} = '{value}'")
            else:
                hdr.append(f"{key:8} = {value}")

        tbl = Table(
            [
                [uid for uid, v in self.data.items() for _ in v["frequency"]],
                np.array([f for v in self.data.values() for f in v["frequency"]])
                * self.units["abscissa"]["unit"],
                np.array([t for v in self.data.values() for t in v["intensity"]])
                * self.units["ordinate"]["unit"],
            ],
            names=["UID", "FREQUENCY", "INTENSITY"],
            meta={"comments": hdr},
        )

        ascii.write(tbl, filename, format="ipac", overwrite=True)

        message(f"WRITTEN: {filename}")

    def plot(self, **keywords) -> None:
        """
        Plot the spectrum.

        """
        import matplotlib as mpl  # type: ignore
        import matplotlib.pyplot as plt  # type: ignore

        _, ax = plt.subplots()
        ax.minorticks_on()
        ax.tick_params(which="major", right="on", top="on", direction="in")
        colors = mpl.colormaps["rainbow"](np.linspace(0, 1, len(self.uids)))
        for d, col in zip(self.data.values(), colors):
            ax.plot(d["frequency"], d["intensity"], color=col)

        ax.set_xlabel(
            self.units["abscissa"]["label"]
            + " ["
            + self.units["abscissa"]["unit"].to_string("latex_inline")
            + "]",
        )
        ax.set_ylabel(
            self.units["ordinate"]["label"]
            + " ["
            + self.units["ordinate"]["unit"].to_string("latex_inline")
            + "]",
        )

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "laboratory"
            plt.savefig(f"{basename}.pdf")
        elif keywords.get("show", False):
            plt.show()
