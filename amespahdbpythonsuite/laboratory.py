#!/usr/bin/env python3

from typing import Optional

import numpy as np

from amespahdbpythonsuite.data import Data


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
