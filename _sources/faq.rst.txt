.. sectnum::
   :start: 5

===
FAQ
===

**Why doesn't the code in this guide work with my installation?**

*You are probably using a version of the AmesPAHdbIDLSuite*
dated prior *February 14, 2015. Please update to a newer
version.*

**Does the AmesPAHdbIDLSuite work with GDL?**

*The AmesPAHdbIDLSuite might work with* `GDL <http://gnudatalanguage.sourceforge.net>`__ *, the open source alternative to IDL.*

**How does the *AmesPAHdbIDLSuite* deal with color output?**

*The AmesPAHdbIDLSuite 'Plot'-methods create and load a color
table containing 15 colors. When outputting to screen, the
device will be set into discrete graphics mode. Upon
completion, both the color table and the display mode will be
restored to their initial settings.*

**Where have the C-routines gone?**

*The C-code has been removed from the AmesPAHdbIDLSuite for
three reasons. First, The NNLS implementation in the suite now
makes use of IDL’s LA_LEAST_SQUARES-routine, which is
implemented using LAPACK and therefore very fast. Second, it
appears that IDL’s READS-routine has seen some significant
speed improvements in recent IDL releases such that the
CONVERT_ASCII C-implementation has become obsolete. Third and
last, making the C-routines compile and reliable across
multiple platforms has proven to be difficult.*

**I am getting schema validation errors, is this an issue and how can I get rid of them?**

*Schema validation errors typically do not cause a problem. The
errors you see could be related to the move from
www.astrochem.org to www.astrochemistry.org for the base URL,
affecting the used XML namespaces. Re-downloading your database
XML-file should solve the problem. As a last resort, you can
disable schema validation by setting the 'Check'-keyword to
'Check=0' when creating the AmesPAHdbIDLSuite-object.*
 
**I am getting an unsupported protocol in URL error, what can I
do?**
 
*With the recent move to using SSL (https) on our domain
(www.astrochemistry.org) IDL's XML-parser fails to download the
linked schema for checking the structure of the database file.
By setting the 'CHECK'-keyword to '0' in the OBJ_NEW-call will
mitigate the problem. Note that you might need to manually
remove the corrupted cache file.*
 
**I am receiving an invalid multi-byte sequence error, what can I
do?**

*It has been reported that on older versions of some operating
systems IDL's XML-parser is unable to properly handle the utf-8
characters embedded in the database XML-file. One solution is
to move to a more modern version of your operation system.
Another possible solution is converting the file to use a
different encoding with the following command.*

.. code:: bash

    xmllint --encode latin1 --output converted.xml file_with_UTF-8.xml

Note that characters that cannot be converted are simply
removed.
