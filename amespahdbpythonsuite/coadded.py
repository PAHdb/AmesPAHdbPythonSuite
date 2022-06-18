#!/usr/bin/env python3

from typing import Optional

import numpy as np

from amespahdbpythonsuite.spectrum import Spectrum


class Coadded(Spectrum):
    """
    AmesPAHdbPythonSuite coadded class

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
        super().__init__(d, **keywords)
        self.__set(d, **keywords)
        return None

    def set(self, d: Optional[dict] = None, **keywords):
        """
        Calls class: :class:`amespahdbpythonsuite.Spectrum.spectrum.set` to parse keywords.

        """
        Spectrum.set(self, d, **keywords)
        self.__set(d, **keywords)

    def __set(self, d: Optional[dict] = None, **keywords) -> None:
        """
        Populate data dictionary helper.

        """
        if isinstance(d, dict):
            if d.get("type", "") == self.__class__.__name__:
                if "weights" not in keywords:
                    self.weights = d["weights"]
                if "averaged" not in keywords:
                    self.averaged = d["averaged"]

        self.weights = keywords.get("weights", dict())
        self.averaged = keywords.get("averaged", False)

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

    def plot(self, **keywords):
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
