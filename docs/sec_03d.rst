.. sectnum::
   :start: 4
   :prefix: 3.

=================================================
Working with fundamental vibrational transitions
=================================================

The 'transitions'-instance exposes the fundamental
vibrational transitions.

.. code:: python

   transitions = pahdb.gettransitionsbyuid(pahdb.search("c<=20 neutral n=2 neutral"))

The 'transitions'-instance's 'print'-method will
print out the associated fundamental vibrational transitions for
each PAH species.

.. code:: python

   transitions.print()

Optionally, the 'str'-keyword can be given to the 'print'-method,
which will return the associated fundamental vibrational
transitions for each PAH species as a single, concatenated string.

.. code:: python
   
   myStr = transitions.print(str=True)

The 'transitions'-instance's 'plot'-method will
display the fundamental vibrational transitions in a 'stick'-plot.
The transitions of each PAH species will be presented in a
different color. 
The 'show' keyword will display the plot on screen, while the 'outfile' keyword will save figure to file.

.. code:: python

   transitions.plot(show=True)
   transitions.plot(outfile='myPlot')

Optionally, the 'wavelength', 'stick', 'oplot', 'legend', and
'color'-keywords can be given to the 'Plot'-method to control the
abscissa, stick representation, overplotting, legend and color,
respectively (*work-in-progress*). 

.. code:: python

   transitions.plot(wavelength=True, color='blue', show=True)

The 'transitions'-instance's 'write'-method will write the fundamental vibrational transitions to file. The transitions of each PAH species will be written to a IPAC table (.tbl). Optionally, a filename can be provided (*work-in-progress*).

.. code:: python

   transitions.write(myFile)
