Semi-Lagrangian Scheme for advection equation
------

``WENO`` class
------

``WENO`` class contains six solvers for 1D/2D/3D rectangle/cubic domain with periodical or non-periodical boundary condition.

The functions are:

``PWENO1D``
``NPWENO1D``
``PWENO2D``
``NPWENO2D``
``PWENO3D``
``NPWENO3D``

For now, constant velocity for 1D has been implemented for CFL in [-.5,.5]. Shifting is not considered right now.

``Interpolation`` class
------
TO BE ADDED.



TEST RESULT
------

``PWENO1D``

===== ========= =============== =================== =============== ============= ================
N     steps     delta_t         ORDER3              ORDER5          ORDER7        ORDER9       
===== ========= =============== =================== =============== ============= ================
  10      100           0.01          0.0754906       0.00620214     0.000515127	   4.4025e-05  
  20      200          0.005          0.0103752      0.000207935     4.39807e-06    9.61453e-08 
  40      400         0.0025         0.00132904	     6.65147e-06     3.53654e-08    1.94407e-10 
  80      800        0.00125        0.000167547	     2.09736e-07     2.79181e-10    3.84144e-13 
 160     1600       0.000625        2.10187e-05	     6.5793e-09	     2.18989e-12	   8.32236e-16 
===== ========= =============== =================== =============== ============= ================
