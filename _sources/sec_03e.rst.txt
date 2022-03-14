.. sectnum::
   :start: 5
   :prefix: 3.

=====================================
Working with Molecular Geometric Data
=====================================

The 'AmesPAHdbPythonSuite_Geometry'-object exposes the molecular
geometric data (*work-in-progress*).

.. code:: python

   uids = pahdb.search("c<=20 neutral n=2 neutral")
   geometry = pahdb.getgeometrybyuid(uids)

The 'geometry'-instance provides both the 'plot' and 'structure'-methods to output the chemical structure of the PAH with provided UID. The first method uses a 2D representation, while the latter a 3D one (*work-in-progress*).

.. code:: python

   # UID 18 = coronene (C24H12)

   geometry.plot(18)
   img = geometry.structure(18)

The 'plot'-method will accept the 'scale', 'thick', and 'rgb'-
keywords to control the scale at which the atoms are to be drawn,
the thickness of the bonds, and if colors should be outputted as
decomposed RGB, respectively. In addition, the 'noerase' and
'position'-keywords control display erasing prior plotting and the
plot position in normalized coordinates, respectively. When the
'resolution'-keyword is set, the plot will be rendered to 
Z-buffer device at given resolution and outputted as an image and
returned to the caller. When the 'save'-keyword is also set, a
PNG-image will be generated as well, which will have the UID of
the PAH molecule under consideration embedded in the filename.
Lastly, the 'angle'-keyword controls rotation of the structure
around the *z*-axis (*work-in-progress*).

.. code:: python

   # UID 18 = coronene (C24H12)

   img = geometry.plot(18, scale=2.0, thick=5.0, resolution=(1200, 1200), save=True)

The 'structure'-method will accept the 'background' and
'frame'-keywords to control background color and display of a
wireframe. The 'resolution'-keyword controls the dimension of the
3D-rendered output image. The 'save'-keyword can be specified to
save the generated image as a PNG-image as well, which will have
the UID of the PAH molecule under consideration embedded in the
filename. If in addition the 'transparent'-keyword is specified as
a color-triplet array, that color will be set to transparent in
the generated PNG-file. The 'axis' and 'angle'-keywords can
control axis and angle of rotation, respectively. With the
'view'-keyword set the generated structure will be rendered on screen.
Lastly, when the 'obj'-keyword is set, the
generated 3D model is returned to the caller (*work-in-progress*).

.. code:: python

   # UID 18 = coronene (C24H12)

img = geometry.structure(18, background=(0, 255, 0), resolution=(600, 600), save=True, transparent=(0, 255, 0), view=True)

In addition, the 'geometry'-instance provides the
'mass', 'rings', and 'area'-methods that return the calculated
mass based on atomic masses, the number of 3-8 membered rings, and
total surface area of the PAH molecules under consideration (*work-in-progress*).

.. code:: python

   masses = geometry.mass()

   rings = geometry.rings()

   areas = geometry.area()

Lastly, the 'inertia'-method provides the moment of inertia
matrices, which are diagonalized with the 'diagonalize'-method,
for the PAH molecules under consideration. The latter method can
be particular useful to ensure proper alignment of the structure
with the view before calling the 'plot' or 'structure'-methods (*work-in-progress*).

.. code:: python

   matrix = geometry.inertia()

   geometry.diagonalize()
