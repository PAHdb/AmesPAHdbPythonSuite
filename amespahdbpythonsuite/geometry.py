#!/usr/bin/env python3

import copy
import numpy as np

from amespahdbpythonsuite.data import Data


class Geometry(Data):
    """
    AmesPAHdbPythonSuite geometry class

    """

    def __init__(self, d=None, **keywords):
        super().__init__(d, **keywords)
        self._atomic_mass = np.array(
            [
                0.0,
                1.007940,
                4.002602,
                6.941000,
                9.012182,
                10.811000,
                12.011000,
                14.006740,
                15.999400,
                18.998404,
                20.179701,
                22.989767,
                24.305000,
                26.981539,
                28.085501,
                30.973763,
                32.066002,
                35.452702,
                39.948002,
                39.098301,
                40.077999,
                44.955910,
                47.880001,
                50.941502,
                51.996101,
                54.938049,
                55.847000,
                58.933201,
                58.693401,
                63.546001,
                65.389999,
                69.723000,
                72.610001,
                74.921593,
                78.959999,
                79.903999,
                83.800003,
                85.467796,
                87.620003,
                88.905853,
                91.223999,
                92.906380,
                95.940002,
                98.000000,
                101.070000,
                102.905502,
                106.419998,
                107.868202,
                112.411003,
                114.820000,
                118.709999,
                121.757004,
                127.599998,
                126.904472,
                131.289993,
                132.905426,
                137.326996,
                138.905502,
                140.115005,
                140.907654,
                144.240005,
                145.000000,
                150.360001,
                151.964996,
                157.250000,
                158.925339,
                162.500000,
                164.930313,
                167.259995,
                168.934204,
                173.039993,
                174.966995,
                178.490005,
                180.947906,
                183.850006,
                186.207001,
                190.199997,
                192.220001,
                195.080002,
                196.966537,
                200.589996,
                204.383301,
                207.199997,
                208.980377,
                209.000000,
                210.000000,
                222.000000,
                223.000000,
                226.024994,
                227.028000,
                232.038101,
                231.035904,
                238.028900,
                237.048004,
                244.000000,
                243.000000,
                247.000000,
                247.000000,
                251.000000,
                252.000000,
                257.000000,
                258.000000,
                259.000000,
                262.000000,
                261.000000,
                262.000000,
                263.000000,
                262.000000,
                265.000000,
                266.000000,
            ]
        )
        self.set(d, **keywords)

    def set(self, d, **keywords) -> None:
        """
        Calls class: :class:`amespahdbpythonsuite.data.Data.set to parse keywords.

        """
        Data.set(self, d, **keywords)

    def get(self) -> dict:
        """
        Assigns class variables from inherited dictionary.

        """
        d = Data.get(self)
        d["type"] = self.__class__.__name__

        return d

    def inertia(self) -> dict:
        """
        Compute the moment of inertia.

        """
        inertia = dict()

        for uid, geometry in self.data.items():

            m = self._atomic_mass[np.array([g["type"] for g in geometry], dtype=int)]
            x = np.array([g["x"] for g in geometry])
            y = np.array([g["y"] for g in geometry])
            z = np.array([g["z"] for g in geometry])

            tensor11 = np.sum(m * (y**2 + z**2))
            tensor22 = np.sum(m * (x**2 + z**2))
            tensor33 = np.sum(m * (x**2 + y**2))
            tensor12 = -np.sum(m * x * y)
            tensor13 = -np.sum(m * x * z)
            tensor23 = -np.sum(m * y * z)

            inertia[uid] = np.array(
                [
                    [tensor11, tensor12, tensor13],
                    [tensor12, tensor22, tensor23],
                    [tensor13, tensor23, tensor33],
                ],
            )

        return inertia

    def diagonalize(self, full=False, equal=False) -> None:
        """
        Diagonalize the moment of inertia and align structure
        with the x-y plane

        """
        if not equal:
            masses = copy.copy(self._atomic_mass)
            masses[[12, 26]] = 0.0
        else:
            masses = np.ones(len(self._atomic_mass))

        for geometry in self.data.values():
            coordinates = np.array(
                [
                    [g["x"] for g in geometry],
                    [g["y"] for g in geometry],
                    [g["z"] for g in geometry],
                ]
            )

            m = np.resize(
                masses[np.array([g["type"] for g in geometry], dtype=int)],
                coordinates.shape,
            )

            coordinates -= np.resize(
                np.sum(coordinates * m, 1) / np.sum(m, 1),
                tuple(reversed(coordinates.shape)),
            ).T

            tensor = np.diag(
                [
                    np.sum(m[0, :] * coordinates[0, :] ** 2),
                    np.sum(m[1, :] * coordinates[1, :] ** 2),
                    np.sum(m[2, :] * coordinates[2, :] ** 2),
                ]
            )

            if full:
                tensor -= np.diag(
                    [
                        np.sum(m[:, 0] * coordinates[:, 0] * coordinates[:, 1]),
                        np.sum(m[:, 2] * coordinates[:, 1] * coordinates[:, 2]),
                        1.0,
                    ]
                )

                tensor += np.diag(np.diag(tensor, 1), -1)

                tensor[0, 2] = -np.sum(m[:, 1] * coordinates[:, 0] * coordinates[:, 2])

                tensor[2, 0] = [0, 2]

            v, w = np.linalg.eig(tensor)

            for i in range(len(geometry)):
                coordinates[:, i] = np.matmul(w, coordinates[:, i])

            coordinates = coordinates[np.argsort(v)[::-1], :]

            for g, i in zip(geometry, range(len(geometry))):
                g["x"] = coordinates[0, i]
                g["y"] = coordinates[1, i]
                g["z"] = coordinates[2, i]

    def mass(self) -> dict:
        """
        Computes molecular mass.

        """

        mass = dict()

        for uid, geometry in self.data.items():
            mass[uid] = np.sum(
                self._atomic_mass[np.array([g["type"] for g in geometry], dtype=int)]
            )

        return mass

    def rings(self) -> list:
        """
        Computes the number of rings per type.

        """

        rings = dict()

        for uid, geometry in self.data.items():

            num = {"three": 0, "four": 0, "five": 0, "six": 0, "seven": 0, "eight": 0}

            x = np.array([g["x"] for g in geometry])

            y = np.array([g["y"] for g in geometry])

            z = np.array([g["z"] for g in geometry])

            numn = np.zeros(x.shape, dtype=int)

            nlist = np.zeros((x.shape[0], 6), dtype=int)

            for i in range(x.shape[0]):
                dd = np.sqrt((x - x[i]) ** 2 + (y - y[i]) ** 2 + (z - z[i]) ** 2)

                bounds = np.where((dd > 0.0) & (dd < 1.7))[0]

                numn[i] = len(bounds)

                nlist[i, 0: numn[i]] = bounds

            iring = np.zeros(9, dtype=int)

            for i in range(x.shape[0]):

                iring[0] = i

                for j in range(numn[i]):

                    i2 = nlist[i, j]

                    if i2 < i:
                        continue

                    iring[1] = i2

                    for k in range(numn[i2]):

                        i3 = nlist[i2, k]

                        if i3 < i:
                            continue

                        iring[2] = i3

                        if i3 != i:

                            for r in range(numn[i3]):

                                i4 = nlist[i3, r]

                                if i4 < i:
                                    continue

                                iring[3] = i4

                                if i4 == i2:
                                    continue

                                if i4 != i or iring[1] < iring[2]:

                                    for m in range(numn[i4]):

                                        i5 = nlist[i4, m]

                                        if i5 < i:
                                            continue

                                        iring[4] = i5

                                        if i5 == i2 or i5 == i3:
                                            continue

                                        if i5 != i or iring[1] < iring[3]:

                                            for n in range(numn[i5]):

                                                i6 = nlist[i5, n]

                                                if i6 < i:
                                                    continue

                                                iring[5] = i6

                                                if i6 == i2 or i6 == i3 or i6 == i4:
                                                    continue

                                                if i6 != i or iring[1] < iring[4]:

                                                    for o in range(numn[i6]):

                                                        i7 = nlist[i6, o]

                                                        if i7 < i:
                                                            continue

                                                        iring[6] = i7

                                                        if (
                                                            i7 == i2
                                                            or i7 == i3
                                                            or i7 == i4
                                                            or i7 == i5
                                                        ):
                                                            continue

                                                        if (
                                                            i7 != i
                                                            or iring[1] < iring[5]
                                                        ):

                                                            for p in range(numn[i7]):

                                                                i8 = nlist[i7, p]

                                                                if i8 < i:
                                                                    continue

                                                                iring[7] = i8

                                                                if (
                                                                    i8 == i2
                                                                    or i8 == i3
                                                                    or i8 == i4
                                                                    or i8 == i5
                                                                    or i8 == i6
                                                                ):
                                                                    continue

                                                                if (
                                                                    i8 != i
                                                                    or iring[1]
                                                                    < iring[6]
                                                                ):

                                                                    for q in range(
                                                                        numn[i8]
                                                                    ):

                                                                        i9 = nlist[
                                                                            i8, q
                                                                        ]

                                                                        if i9 < i:
                                                                            continue

                                                                        iring[8] = i9

                                                                        if (
                                                                            i9 == i2
                                                                            or i9 == i3
                                                                            or i9 == i4
                                                                            or i9 == i5
                                                                            or i9 == i6
                                                                            or i9 == i7
                                                                        ):
                                                                            continue

                                                                        if (
                                                                            i9 != i
                                                                            or iring[1]
                                                                            < iring[7]
                                                                        ):
                                                                            continue

                                                                        else:
                                                                            num[
                                                                                "eight"
                                                                            ] += 1
                                                                else:
                                                                    num["seven"] += 1

                                                        else:
                                                            num["six"] += 1

                                                else:
                                                    num["five"] += 1

                                        else:
                                            num["four"] += 1

                                else:
                                    num["three"] += 1

            rings[uid] = num

        return rings

    def area(self) -> dict:
        """
        Computes the area.

        """

        per_ring = {
            "three": 0.848,
            "four": 1.96,
            "five": 3.37,
            "six": 5.09,
            "seven": 7.12,
            "eight": 9.46,
        }

        area = dict()

        rings = self.rings()
        for uid, nr in rings.items():
            a = 0.0
            for r, n in nr.items():
                a += per_ring[r] * n
            area[uid] = a

        return area
