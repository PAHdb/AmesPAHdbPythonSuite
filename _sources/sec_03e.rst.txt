.. sectnum::
   :start: 5
   :prefix: 3.

=====================================
Working with Molecular Geometric Data
=====================================

The 'AmesPAHdbIDLSuite_Geometry'-object exposes the molecular
geometric data.

.. code:: idl

   geometry = pahdb->getGeometryByUID( $
              pahdb->Search("c<=20 neutral n=2 neutral"))

The 'AmesPAHdbIDLSuite_Geometry'-object provides both the 'Plot'
and 'Structure'-methods to output the chemical structure of the
PAH with provided UID. The first method uses IDL's PLOT-procedures
to create the output, while the latter IDL object graphics.

.. code:: idl

   ; UID 18 = coronene (C24H12)

   void = geometry->Plot(18)

   img = geometry->Structure(18)

The 'Plot'-method will accept the 'Scale', 'Thick', and 'RGB'-
keywords to control the scale at which the atoms are to be drawn,
the thickness of the bonds, and if colors should be outputted as
decomposed RGB, respectively. In addition, the 'NoErase' and
'Position'-keywords control display erasing prior plotting and the
plot position in normalized coordinates, respectively. When the
'Resolution'-keyword is set, the plot will be rendered to IDL's
Z-buffer device at given resolution and outputted as an image and
returned to the caller. When the 'Save'-keyword is also set, a
PNG-image will be generated as well, which will have the UID of
the PAH molecule under consideration embedded in the filename.
Lastly, the 'Angle'-keyword controls rotation of the structure
around the *z*-axis.

.. code:: idl

   ; UID 18 = coronene (C24H12)

   img = geometry->Plot(18, Scale=2.0, Thick=5.0, $
                           Resolution=[1200, 1200], /Save)

The 'Structure'-method will accept the 'Background' and
'Frame'-keywords to control background color and display of a
wireframe. The 'Resolution'-keyword controls the dimension of the
3D-rendered output image. The 'Save'-keyword can be specified to
save the generated image as a PNG-image as well, which will have
the UID of the PAH molecule under consideration embedded in the
filename. If in addition the 'Transparent'-keyword is specified as
a color-triplet array, that color will be set to transparent in
the generated PNG-file. The 'Axis' and 'Angle'-keywords can
control axis and angle of rotation, respectively. With the
'View'-keyword set the generated structure will be displayed using
IDL's XOBJVIEW routine. Lastly, when the 'Obj'-keyword is set, the
generated IDLgrModel-object is returned to the caller.

.. code:: idl

   ; UID 18 = coronene (C24H12)

   img = geometry->Structure(18, Background=[0, 255, 0], $
                                 Resolution=[600, 600], $
                                 /Save, Transparent=[0, 255, 0], $
                                 /View)

In addition, the 'AmesPAHdbIDLSuite_Geometry'-object provides the
'Mass', 'Rings', and 'Area'-methods that return the calculated
mass based on atomic masses, the number of 3-8 membered rings, and
total surface area of the PAH molecules under consideration.

.. code:: idl

   masses = geometry->Mass()

   rings = geometry->Rings()

   areas = geometry->Area()

Lastly, the 'Inertia'-method provides the moment of inertia
matrices, which are diagonalized with the 'Diagonalize'-method,
for the PAH molecules under consideration. The latter method can
be particular useful to ensure proper alignment of the structure
with the view before calling the 'Plot' or 'Structure'-methods.

.. code:: idl

   matrix = geometry->Inertia()

   geometry->Diagonalize