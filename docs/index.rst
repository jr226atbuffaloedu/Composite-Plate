.. Composite Plate documentation master file, created by
   sphinx-quickstart on Wed Nov 18 11:19:19 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Composite Plate's documentation!
===========================================

Contents:

.. toctree::
   :maxdepth: 2

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

How to Use Composite Plate
==========================

The goal of this module is to provide you with a way to determine how a 
composite plate will respond to either stresses or strains.  In order to
do this, you must define the composite plate in the following way.

#. Define a ply material using :py:class:`composite_plate.Ply`
#. Use the ply and a given orientation to define a laminae using :py:class:`composite_plate.Laminae`
#. One or more of the laminae in a list to define a laminate using :py:class:`composite_plate.Laminate`

Example: ::
    
    >>> E1      = 133.44;   # Modulus in fiber direction        [GPa]
    >>> E2      = 8.78;     # Modulus in transverse direction   [GPa]
    >>> G12     = 3.254;    # Shear Modulus                     [GPa]
    >>> nu12    = 0.26;     # Poissons Ratio (principal)        [   ]
    >>> h       = 6;        # Ply thickness                     [mm]

    >>> ply = composite_plate.Ply(E1=E1,E2=E2,G12=G12,nu12=nu12,h=h)
    
    >>> theta_deg = 45;     # Orientation (+ from X)            [deg]
    >>> theta_rad = theta_deg*np.pi/180.0;    #                 [rad]
    
    >>> laminae1 = composite_plate.Laminae(ply, theta_rad);            
    >>> laminae2 = composite_plate.Laminae(ply,-theta_rad);
    
    >>> laminate = composite_plate.Laminate([laminae1, laminae2])
    
After the laminate is created, you can then apply strains and see the 
resulting stress using :py:func:`composite_plate.Laminate.applied_strain`.

Example: ::
    
    >>> Epsilon = numpy.matrix([]);
    >>> Kappa = numpy.matrix([]);
    >>> force_and_moment_dict = laminate.applied_strain(Epsilon, Kappa)
    >>> N = force_and_moment_dict('N')
    >>> print("The forces are Nxx = {0}, Nyy = {1}, Nxy = {2}".format(N[1],N[2],N[3]))
    
    >>> M = force_and_moment_dict('M')
    >>> print("The moments are Mxx = {0}, Myy = {1}, Mxy = {2}".format(M[1],M[2],M[3]))


Modules
=======

.. automodule:: composite_plate
    :members:
    :special-members:
