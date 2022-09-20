#!/usr/bin/env python3
import copy
from typing import Optional

import numpy as np
from scipy.spatial import KDTree  # type: ignore
from scipy.spatial.distance import euclidean as edist  # type: ignore

from amespahdbpythonsuite.amespahdb import AmesPAHdb
from amespahdbpythonsuite.data import Data

message = AmesPAHdb.message


class Geometry(Data):
    """
    AmesPAHdbPythonSuite geometry class

    """

    def __init__(self, d: Optional[dict] = None, **keywords) -> None:
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

        return f"AmesPAHdbPythonSuite Geometry instance.\n" f"{self.uids=}"

    def plot(self, uid: int, **keywords) -> None:
        """
        Plot the structure

        """

        import matplotlib.pyplot as plt  # type: ignore

        # atom_names = ["H", "C", "N", "O", "Mg", "Si", "Fe"]

        atom_numbers = [1, 6, 7, 8, 12, 14, 26]

        atom_colors = ["grey", "black", "blue", "red", "green", "pink", "orange"]

        atom_symsize = np.array([1, 2, 3, 3, 3.5, 4, 4]) * keywords.get("scale", 100.0)

        g = self.data[uid]
        ng = len(g)

        numn = np.zeros(ng, dtype=int)
        nlist = np.full((ng, 6), -1, dtype=int)

        px = np.array([g["x"] for g in self.data[uid]])
        py = np.array([g["y"] for g in self.data[uid]])
        pz = np.array([g["z"] for g in self.data[uid]])

        m = np.max([np.max(px), np.max(py)])

        for x, y, z, i in zip(px, py, pz, range(ng)):
            dd = np.sqrt((px - x) ** 2 + (py - y) ** 2 + (pz - z) ** 2)
            sel = np.where((dd < 1.6))
            nsel = len(sel[0])
            numn[i] = nsel
            nlist[i, 0:nsel] = sel[0]

        _, ax = plt.subplots()
        ax.set_xlim(-m, m)
        ax.set_ylim(-m, m)

        for x, y, i in zip(px, py, range(ng)):
            for j in range(numn[i]):
                ax.plot([x, px[nlist[i, j]]], [y, py[nlist[i, j]]], c="black", lw=5)

                numn[nlist[i, j]] -= 1
                if numn[nlist[i, j]] > 0:
                    nlist[nlist[i, j], np.where(nlist[nlist[i, j], :] == i)[0]] = -1
                    nlist[nlist[i, j], :] = nlist[
                        nlist[i, j], np.argsort(nlist[nlist[i, j], :])[::-1]
                    ]

        pt = np.array([g["type"] for g in self.data[uid]])

        for i in range(len(atom_numbers)):
            ii = np.where(pt == atom_numbers[i])[0]
            if len(ii) > 0:
                ax.scatter(
                    px[ii],
                    py[ii],
                    marker="o",
                    s=atom_symsize[i],
                    zorder=2,
                    c=atom_colors[i],
                )

        ratio = 1.0
        x_left, x_right = ax.get_xlim()
        y_low, y_high = ax.get_ylim()
        ax.set_aspect(abs((x_right - x_left) / (y_low - y_high)) * ratio)
        plt.axis("off")

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "laboratory"
            plt.savefig(f"{basename}.pdf")
        elif keywords.get("show", False):
            plt.show()

    def structure(self, uid: int, **keywords) -> None:
        """
        Render the 3D structure.

        """

        import vtkmodules.vtkInteractionStyle  # type:  ignore
        import vtkmodules.vtkRenderingOpenGL2  # type:  ignore # noqa:  F401
        from vtk import vtkTransform  # type:  ignore
        from vtkmodules.vtkFiltersSources import (  # type:  ignore
            vtkCylinderSource,
            vtkSphereSource,
        )
        from vtkmodules.vtkRenderingCore import (  # type:  ignore
            vtkRenderWindow,
            vtkActor,
            vtkPolyDataMapper,
            vtkRenderer,
            vtkRenderWindowInteractor,
            vtkWindowToImageFilter,
        )  # type:  ignore
        from vtkmodules.vtkIOImage import vtkPNGWriter  # type:  ignore

        atom_colors = {
            1: [0.78, 0.78, 0.78],
            6: [0.11, 0.11, 0.11],
            7: [0.19, 0.31, 0.97],
            8: [1.0, 0.05, 0.05],
            12: [0.54, 1.0, 0.0],
            14: [0.94, 0.78, 0.63],
            26: [0.89, 0.40, 0.20],
        }

        atom_radii = {1: 0.25, 6: 0.5, 7: 0.75, 8: 0.75, 12: 0.875, 14: 1, 26: 1}

        scale = keywords.get("scale", 1.0)
        for number in atom_radii:
            atom_radii[number] *= scale

        g = self.data[uid]
        ng = len(g)

        numn = np.zeros(ng, dtype=int)
        nlist = np.full((ng, 6), -1, dtype=int)

        px = np.array([g["x"] for g in self.data[uid]])
        py = np.array([g["y"] for g in self.data[uid]])
        pz = np.array([g["z"] for g in self.data[uid]])
        pt = np.array([g["type"] for g in self.data[uid]])

        for x, y, z, i in zip(px, py, pz, range(ng)):
            dd = np.sqrt((px - x) ** 2 + (py - y) ** 2 + (pz - z) ** 2)
            sel = np.where((dd < 1.6))
            nsel = len(sel[0])
            numn[i] = nsel
            if nsel == 0:
                continue
            if nsel > 6:
                ii = np.argsort(dd[sel])
                sel = tuple([sel[0][ii[0:5]]])
                nsel = 6

            nlist[i, 0:nsel] = sel[0]

        actors = list()

        for x, y, z, t, i in zip(px, py, pz, pt, range(ng)):
            for j in range(numn[i]):
                numn[nlist[i, j]] -= 1
                if numn[nlist[i, j]] > 0:
                    nlist[nlist[i, j], np.where(nlist[nlist[i, j], :] == i)[0]] = -1
                    nlist[nlist[i, j], :] = nlist[
                        nlist[i, j], np.argsort(nlist[nlist[i, j], :])[::-1]
                    ]

                if nlist[i, j] < 0:
                    continue

                vec = np.array(
                    (x - px[nlist[i, j]], y - py[nlist[i, j]], z - pz[nlist[i, j]])
                )

                norm = np.linalg.norm(vec)

                cylinder = vtkCylinderSource()
                cylinder.SetResolution(32)
                cylinder.SetRadius(0.1)
                cylinder.SetHeight(norm)
                cylinderMapper = vtkPolyDataMapper()
                cylinderMapper.SetInputConnection(cylinder.GetOutputPort())
                cylinderActor = vtkActor()
                cylinderActor.SetMapper(cylinderMapper)
                if t == 1 or pt[nlist[i, j]] == 1:
                    cylinderActor.GetProperty().SetColor(0.78, 0.78, 0.78)
                else:
                    cylinderActor.GetProperty().SetColor(0.11, 0.11, 0.11)
                cylinderTransform = vtkTransform()
                cylinderTransform.Identity()
                cylinderTransform.PostMultiply()
                cylinderTransform.Translate(0.0, norm / 2.0, 0.0)
                cylinderTransform.RotateX(90.0)
                angle = 0.0
                if norm != 0.0:
                    angle = 180.0 * np.arccos(-vec[2] / norm) / np.pi
                    if angle == 180.0:
                        cylinderTransform.Translate(0.0, 0.0, -norm)
                cylinderTransform.RotateWXYZ(angle, vec[1], -vec[0], 0.0)
                cylinderTransform.Translate(x, y, z)
                cylinderActor.SetUserTransform(cylinderTransform)
                actors.append(cylinderActor)

        if not keywords.get("frame", False):
            for x, y, z, t in zip(px, py, pz, pt):
                sphere = vtkSphereSource()
                sphere.SetThetaResolution(32)
                sphere.SetPhiResolution(32)
                sphere.SetCenter(x, y, z)
                sphere.SetRadius(atom_radii[t])
                sphereMapper = vtkPolyDataMapper()
                sphereMapper.SetInputConnection(sphere.GetOutputPort())
                sphereActor = vtkActor()
                sphereActor.SetMapper(sphereMapper)
                sphereActor.GetProperty().SetColor(atom_colors[t])

                sphereActor.GetProperty().SetSpecular(0.25)

                actors.append(sphereActor)

        renderer = vtkRenderer()
        renderer.SetBackground(0.24, 0.24, 0.24)

        for actor in actors:
            renderer.AddActor(actor)

        renderWindow = vtkRenderWindow()
        if not keywords.get("show", False):
            renderWindow.SetOffScreenRendering(True)
            if keywords.get("transparent", False):
                renderWindow.AlphaBitPlanesOn()
        else:
            renderWindow.SetWindowName(self.__class__.__name__)
            interactor = vtkRenderWindowInteractor()
            interactor.SetRenderWindow(renderWindow)

        renderWindow.AddRenderer(renderer)
        renderWindow.SetSize(1024, 1024)

        renderWindow.Render()
        if renderWindow.GetAlphaBitPlanes() == 0 and keywords.get("transparent", False):
            message("TRANSPARENCY NOT SUPPORTED")

        basename = keywords.get("save")
        if basename:
            if not isinstance(basename, str):
                basename = "structure"
            windowToImageFilter = vtkWindowToImageFilter()
            windowToImageFilter.SetInput(renderWindow)
            if keywords.get("transparent", False):
                windowToImageFilter.SetInputBufferTypeToRGBA()
            windowToImageFilter.SetScale(1)
            windowToImageFilter.Update()
            writer = vtkPNGWriter()
            writer.SetInputConnection(windowToImageFilter.GetOutputPort())
            writer.SetFileName(f"{basename}.png")
            writer.Write()

        if keywords.get("show", False):
            interactor.Initialize()
            interactor.Start()

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

    def diagonalize(self, full: bool = False, equal: bool = False) -> None:
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

    def rings(self) -> dict:
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

    def bec(self) -> dict:
        """Convert PAH carbon and hydrogen positions to a boundary-edge
        code. Only makes sense to use this for regular PAHs. Return of
        a tuple with a list of boundary carbons traversed in order and
        the boundary-edge code. By Dr. Joseph E. Roser
        <Joseph.E.Roser@nasa.gov>"""

        becodes = dict()
        for uid, geometry in self.data.items():
            becode = self.__becode(
                [(g["x"], g["y"], g["z"]) for g in geometry if g["type"] == 6],
                [(g["x"], g["y"], g["z"]) for g in geometry if g["type"] == 1],
            )
            becodes[uid] = becode[1]

        return becodes

    def __becode(
        self,
        allcarbons: list[tuple[float, float, float]],
        allhydrogens: list[tuple[float, float, float]],
    ) -> tuple:
        """The star of the show here is the function PAHbecode, which answers
        the challenge of converting a list of PAH carbon atom and
        hydrogen atom positions into a boundary-edge code of a PAH
        molecule. Due to the possibilities of ring distortion,
        non-planarity, ambiguity of orientation in space, and
        topological challenges of bizarrely shaped PAH molecules,
        we're going to do the best we can here. Please check the
        results.

        The function indicator simply calculates the sign of a scalar
        triple product of three three-vectors. The function
        pcone_to_be converts a PC-1 code to a boundary-edge
        code. These are relatively simple helper functions, but
        useable if you need them.

        Author: Dr. Joseph E. Roser <Joseph.E.Roser@nasa.gov>

        Given a list of 3D carbon atom positions and a list of 3D
        hydrogen atom positions, this function will do its best to
        return (a) the boundary carbon atoms traversed in order and
        (b) the boundary-edge code string

        """

        # Step 1: Find any boundary edge.
        carbontree = KDTree(allcarbons)

        # Compute the center of mass of the carbon atom distribution
        carbon_center = [
            sum([Catom[i] for Catom in allcarbons]) / len(allcarbons)
            for i in range(0, 3)
        ]

        # Rank distances of all carbon atoms from the center point
        _, indices = carbontree.query(x=carbon_center, k=len(allcarbons))

        # Initialize boundary_carbons list
        boundary_carbons = [allcarbons[indices[-1]]]
        _, nindex = carbontree.query(x=boundary_carbons[0], k=2)
        boundary_carbons.append(allcarbons[nindex[-1]])

        # Step 2: We need to specify a vector that is roughly normal
        # to the surface of the PAH molecule.
        v0 = np.subtract(boundary_carbons[0], carbon_center)
        v1 = np.subtract(boundary_carbons[1], carbon_center)
        normal = np.cross(v1, v0)

        # Step 3: A length scale for bounding carbon-carbon nearest
        # neighbor searches
        dscale = 1.366025404 * edist(boundary_carbons[0], boundary_carbons[1])

        # Step 4: Which carbon atom neighbor of a given hydrogen atom
        # is really the one that it is bound to?
        hydrogenated_carbons = [
            allcarbons[index]
            for index in map(lambda y: carbontree.query(x=y)[1], allhydrogens)
        ]

        # Step 5: Traverse the boundary step by step and add in the
        # boundary carbon atoms one by one.
        for _ in range(0, 2 * len(allhydrogens) - 8):
            # Step 6: Look for carbon atoms one "hex vertex" away from
            # the current end point of the boundary traversal.
            _, indexlist = carbontree.query(
                x=boundary_carbons[-1], k=4, distance_upper_bound=dscale
            )

            # Reject duplicate boundary points and invalid indicies in
            # indexlist
            trial_carbons = list(
                filter(
                    lambda x: x not in boundary_carbons,
                    [allcarbons[item] for item in indexlist if item < len(allcarbons)],
                )
            )

            # Step 7: Update the boundary carbon atoms list.
            if len(trial_carbons) == 0:
                message("BEC ERROR: SEARCH FOUND TOO FEW NEAREST-NEIGHBOR POINTS")
                break
            elif len(trial_carbons) == 1:
                boundary_carbons.extend(trial_carbons)
            elif len(trial_carbons) == 2:
                # Compute an edge vector for the most recently added boundary edge
                edgevector = np.subtract(boundary_carbons[-1], boundary_carbons[-2])

                # If the boundary end point is hydrogenated, we move
                # "clockwise", otherwise "counter-clockwise"
                for Catom in trial_carbons:
                    trialvector = np.subtract(Catom, boundary_carbons[-1])
                    ivalue = self.__indicator(trialvector, edgevector, normal)
                    if not np.logical_xor(
                        ivalue == 1.0, boundary_carbons[-1] in hydrogenated_carbons
                    ):
                        boundary_carbons.append(Catom)
                        break
            else:
                message("BEC ERROR: SEARCH FOUND TOO MANY NEAREST-NEIGHBOR POINTS")
                break

        # Step 8: We can then compute the boundary-edge code
        if len(allcarbons) == 6:
            message("BEC NOTICE: INPUT SUGGEST THE TRIVIAL PAH BENZENE WITH BEC '6'")
            becode = "6"
        else:
            # Compute a PC-1 code from boundary_carbons
            pcone_code = [
                "0" if Catom in hydrogenated_carbons else "1"
                for Catom in boundary_carbons
            ]

            # Compute the boundary-edge code and its reversed
            # equivalent
            becode = self.__pcone_to_be(pcone_code)
            pcone_code.reverse()
            reverse_becode = self.__pcone_to_be(pcone_code)

            # Convert the boundary-edge code to its lexicographically
            # maximum equalivalent.
            if len(becode) > 2:
                forward_best = max(
                    [becode[i:] + becode[:i] for i in range(0, len(becode))]
                )
                reverse_best = max(
                    [
                        reverse_becode[i:] + reverse_becode[:i]
                        for i in range(0, len(becode))
                    ]
                )
                becode = max(forward_best, reverse_best)

        # Step 9: Return a tuple containing the boundary-edge list
        # (with an extra point for boundary closure) and the computed
        # boundary edge code
        boundary_carbons.append(boundary_carbons[0])
        return (boundary_carbons, becode)

    def __indicator(
        self,
        v1: np.ndarray,
        v2: np.ndarray,
        normal: np.ndarray,
    ) -> float:
        """Returns the sign of normal dot (v1 cross v2) assuming that these
        are 3-element sequences of some kind. By Dr. Joseph E. Roser
        <Joseph.E.Roser@nasa.gov

        """
        value = 0.0
        value += normal[0] * (v1[1] * v2[2] - v1[2] * v2[1])
        value += normal[1] * (v1[2] * v2[0] - v1[0] * v2[2])
        value += normal[2] * (v1[0] * v2[1] - v1[1] * v2[0])
        return np.sign(value)

    def __pcone_to_be(self, pcone_code: list[str]) -> str:
        """Converts the PC-1 code of a PAH to its boundary-edge code.  By
        Dr. Joseph E. Roser <Joseph.E.Roser@nasa.gov

        """
        becode = ""
        csum = 0
        x = pcone_code.index("1")
        for item in pcone_code[x + 1:] + pcone_code[: x + 1]:
            if item == "0":
                csum += 1
            else:
                becode += str(csum + 1)
                csum = 0
        return becode
