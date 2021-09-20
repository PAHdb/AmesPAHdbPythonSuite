.. sectnum::
   :start: 10
   :prefix: 3.

=========================
 Graphical User Interface
=========================
         
Utilizing IDL's widget programming capabilities, a simple
graphical user interface (GUI) is available to explore the
contents of a database XML-file. The GUI provides an overview of
the embedded PAH species and the fundamental vibrational
transitions - both tabulated and graphically, and 3D-rotatable
chemical structure of the selected species. A database XML-file
can be loaded through the 'File→Open' menu-option.

.. code:: idl

   amespahdbidlsuite_gui

.. figure:: figures/gui@2x.png
   :align: center

   Screenshot of the AmesPAHdbIDLsuite GUI.

**NB** The GUI generates all necessary assets on-the-fly, which
results in a considerable start-up-time, notably when loading a
large database.
