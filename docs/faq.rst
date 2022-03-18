.. sectnum::
   :start: 5

===
FAQ
===

**What Python version is required to run the package?**

*The AmesPAHdbPythonSuite requires a Python > 3.6 installation.*


**I am getting schema validation errors, is this an issue and how can I get rid of them?**

*Schema validation errors typically do not cause a problem. The
errors you see could be related to the move from
www.astrochem.org to www.astrochemistry.org for the base URL,
affecting the used XML namespaces. Re-downloading your database
XML-file should solve the problem. As a last resort, you can
disable schema validation by setting the 'check'-keyword to
'check=False' when creating the AmesPAHdbPythonSuite-object.*
 
**I am getting an unsupported protocol in URL error, what can I
do?**
 
*With the recent move to using SSL (https) on our domain
(www.astrochemistry.org) the XML-parser fails to download the
linked schema for checking the structure of the database file.
By setting the 'check'-keyword to 'False' in the AmesPAHdb-call will
mitigate the problem. Note that you might need to manually
remove the corrupted cache file.*

**Is multiprocessing supported?**

*Multiprocessing is currently supported for the `cascade` emission model and the `convolve` method of the `transitions`-instance.*

.. code:: python

    transitions.cascade(6 * 1.603e-12, multiprocessing=True)
    transitions.convolve(grid, fwhm=15.0, gaussian=True, multiprocessing=True)