.. sectnum::
   :start: 7
   :prefix: 3.

======================================================
Working Outside the *AmesPAHdbIDLSuite* Object Classes
======================================================

The *AmesPAHdbPythonSuite* object classes all provide the
'get'-method that allows extraction of the object's internal data
into a dictionary.

.. code:: python

   transitions_dict = transitions.get()

This allows, for example, more control over the presentation of
the data.

.. code:: python

   import matplotlib.pyplot as plt

   plt.plot(transitions_dict['data']['frequency'], transitions_dict['data']['intensity'])
   plt.xlabel(transitions_dict['units']['x']['str'])
   plt.ylabel('integrated intensity [cm!U-2!N mol!U-1!N]')
   plt.show()

Subsequently this dictionary can be manipulated and even
be set to the object, which will try altering its internal state
to reflect that of the manipulated data.

.. code:: python

   transitions.set(transitions_dict)

In addition, the 'set'-method accepts keywords to alter its
internal state.

.. code:: python

   transitions.set(data=transitions_dict['data'])

Lastly, it is also possible to create new intances initialized with an appropriate data representation or through keywords.

.. code:: python

   transitions = transitions(transitions_dict)

   transitions = trantions(type=transitions_dict['type'], \
                           version=transitions_dict['version'], \
                           pahdb=pahdb, \
                           data=transitions_dict['data'], \
                           uids=transitions_dict['uids'], \
                           model=transitions_dict['model'], \
                           units=transitions_dict['units'])