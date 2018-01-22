How to calculate DNA elasticity?
================================

DNA elasticity governs the DNA conformational fluctations. Therefore, its study gives insight into the DNA dynamics and
stiffness. There are two level at which elastic properties can be determined.

1) Global Elastic Properties : It can be used for determining the stiffness of longer DNA segment.
2) Local Elastic Properties: It can be used to study stiffness of a base-step (two adjacent nucleotides).


Global Elastic Properties
-------------------------

Global elastic properties include the motions governing bending, stretching and twisting motions of the DNA. To determine
these properties, these motions for DNA segment need to be quantified from the simulations.

1) Bending motion : In dnaMD, bending motion is decomposed into two orthogonal planes. These motions can be quantified
   using bending angles. Bending angle is calculated between the tangent of global helical axis present at the end of
   the DNA segment.

2) Stretching motion : It is sum of the helical-rise calculated from 3DNA/do_x3dna package.

3) Twisting motion : It is sum of the helical-twist calculated from 3DNA/do_x3dna package.

These quantities are calculated from the simulation trajectories, and subsequently, elastic matrix is calculated using
covariance matrix.

.. graphviz::

   digraph {
        a -> b-> c -> d -> e
        g -> h -> d -> f

        a [label = "Free DNA Simulation Trajectories"];
        b [label = "Bending angles, Contour length and Twist angle"];
        c [label = "Averages and Covariance matrix"];
        d [label = "Elastic Constant Matrix"];
        e [label = "Modulus Matrix"];
        f [label = "Deformation free energy of bound DNA"];

        g [label = "Bound DNA Simulation Trajectories"]
        h [label = "Bending angles, Contour length and Twist angle"];
   }


**Covariance matrix**

.. math::
    \mathbf{C} = \begin{bmatrix}
        \langle (\theta^{x}_{i} - \theta^{x}_0) (\theta^{x}_{i} - \theta^{x}_0) \rangle   &
        \langle (\theta^{x}_{i} - \theta^{x}_0) (\theta^{y}_{i} - \theta^{y}_0) \rangle   &
        \langle (\theta^{x}_{i} - \theta^{x}_0) (L_i - L_0) \rangle                       &
        \langle (\theta^{x}_{i} - \theta^{x}_0) (\phi_i - \phi_0) \rangle                 \\
        \langle (\theta^{x}_{i} - \theta^{x}_0) (\theta^{y}_{i} - \theta^{y}_0) \rangle   &
        \langle (\theta^{y}_{i} - \theta^{y}_0) (\theta^{y}_{i} - \theta^{y}_0) \rangle   &
        \langle (\theta^{y}_{i} - \theta^{y}_0) (L_i - L_0) \rangle                       &
        \langle (\theta^{y}_{i} - \theta^{y}_0) (\phi_i - \phi_0) \rangle                 \\
        \langle (\theta^{x}_{i} - \theta^{x}_0) (L_i - L_0) \rangle                       &
        \langle (\theta^{y}_{i} - \theta^{y}_0) (L_i - L_0) \rangle                       &
        \langle (L_i - L_0) (L_i - L_0) \rangle                                           &
        \langle (L_i - L_0) (\phi_i - \phi_0) \rangle                                     \\
        \langle (\theta^{x}_{i} - \theta^{x}_0) (\phi_i - \phi_0) \rangle                 &
        \langle (\theta^{y}_{i} - \theta^{y}_0) (\phi_i - \phi_0) \rangle                 &
        \langle (L_i - L_0) (\phi_i - \phi_0) \rangle                                     &
        \langle (\phi_i - \phi_0) (\phi_i - \phi_0) \rangle
    \end{bmatrix}


where:
    * :math:`\theta^{x}_{i}` - Bending angle in first plane at :math:`i` th frame
    * :math:`\theta^{x}_0` - Average Bending angle in first plane from all frames
    * :math:`\theta^{y}_{i}` - Bending angle in second plane at :math:`i` th frame
    * :math:`\theta^{y}_0` - Average Bending angle in second plane from all frames
    * :math:`L_i` - Contour length (sum of helical-rise) at :math:`i` th frame
    * :math:`L_0` - Average Contour length from all frames
    * :math:`\phi_i` - Twist angle (sum of helical-twist) at :math:`i` th frame
    * :math:`\phi_0` - Average twist angle from all frames

**Elastic Constant Matrix:**

.. math::
    \text{Elastic matrix} = \mathbf{K} = k_BT\mathbf{C^{-1}} = \begin{bmatrix}
        K_{Bx}       & K_{Bx,By} & K_{Bx,S} & K_{Bx,T} \\
        K_{Bx,By}    & K_{By}    & K_{By,S} & K_{By,T} \\
        K_{Bx,S}     & K_{By,S}  & K_{S}    & K_{S,T} \\
        K_{Bx,T}     & K_{Bx,T}  & K_{S,T}  & K_{T}
    \end{bmatrix}

where:
    * :math:`k_B` - Boltzmann constant
    * :math:`T` - Temperature
    * :math:`\mathbf{C}` - Covariance matrix
    * :math:`Bx` - Bending motion in first plane
    * :math:`By` - Bending motion in second plane
    * :math:`S` - Stretching motion
    * :math:`T` - Twisting motion

**Modulus Matrix**

.. math::
    \text{modulus matrix} =
    \begin{bmatrix}
    M_{Bx}       & M_{Bx,By} & M_{Bx,S} & M_{Bx,T} \\
    M_{Bx,By}    & M_{By}    & M_{By,S} & M_{By,T} \\
    M_{Bx,S}     & M_{By,S}  & M_{S}    & M_{S,T} \\
    M_{Bx,T}     & M_{Bx,T}  & M_{S,T}  & M_{T}
    \end{bmatrix}
    = 4.1419464 \times \mathbf{K}  \times L_0

where:
    * :math:`M_{Bx}` - Bending-1 stiffness in one plane
    * :math:`M_{By}` - Bending-2 stiffness in another orthogonal plane
    * :math:`M_{S}` - Stretch Modulus
    * :math:`M_{T}` - Twist rigidity
    * :math:`M_{Bx,By}` - Bending-1 and Bending-2 coupling
    * :math:`M_{By,S}` - Bending-2 and stretching coupling
    * :math:`M_{S,T}` - Stretching Twsiting coupling
    * :math:`M_{Bx,S}` - Bending-1 Stretching coupling
    * :math:`M_{By,T}` - Bending-2 Twisting coupling
    * :math:`M_{Bx,T}` - Bending-1 Twisting coupling


Global Deformation Free Energy
------------------------------

Deformation free energy can be calculated using following equation with either elastic constant matrix or modulus matrix.

.. math::
    G = \frac{1}{2}\mathbf{xKx^T}

**or**

.. math::
    G = \frac{1}{2L_0}\mathbf{xMx^T}

where:
    * :math:`\mathbf{x} =  \begin{bmatrix} (\theta^{x}_{i} - \theta^{x}_0)    & (\theta^{y}_{i} - \theta^{y}_0) & (L_i - L_0) & (\phi_i - \phi_0) \end{bmatrix}`


In dnaMD, deformation free energy is directly calculated from :math:`\mathbf{K}` instead of :math:`\mathbf{M}`.


Local Elastic Properties
------------------------

For a single base-step or few base-steps, local elastic properties can be calculated using either local base-step parameters
or helical base-step parameters. These parameters quantify six degree of freedom along which a base-step can fluctuates.
Therefore, elasticity along these six degrees of freedom and their coupling can be calculated under the harmonic approximation.

In case of :ref:`base-step-image`: Shift (:math:`Dx`), Slide (:math:`Dy`), Rise (:math:`Dz`),
Tilt (:math:`\tau`), Roll (:math:`\rho`) and Twist (:math:`\omega`), following elastic matrix is calculated.

.. math::

    \mathbf{K}_{base-step} = \begin{bmatrix}
    K_{Dx}        & K_{Dx,Dy}      & K_{Dx,Dz}      & K_{Dx,\tau}      & K_{Dx,\rho}      & K_{Dx,\omega} \\
    K_{Dx,Dy}     & K_{Dy}         & K_{Dy,Dz}      & K_{Dy,\tau}      & K_{Dy,\rho}      & K_{Dy,\omega} \\
    K_{Dx,Dz}     & K_{Dy,Dz}      & K_{Dz}         & K_{Dz,\tau}      & K_{Dz,\rho}      & K_{Dz,\omega} \\
    K_{Dx,\tau}   & K_{Dy,\tau}    & K_{Dz,\tau}    & K_{\tau}         & K_{\tau, \rho}   & K_{\tau,\omega} \\
    K_{Dx,\rho}   & K_{Dy,\rho}    & K_{Dz,\rho}    & K_{\tau, \rho}   & K_{\rho}         & K_{\rho,\omega} \\
    K_{Dx,\omega} & K_{Dy,\omega}  & K_{Dz,\omega}  & K_{\tau, \omega} & K_{\rho, \omega} & K_{\omega} \\
    \end{bmatrix}


In case of :ref:`helical-base-step-image`: x-displacement (:math:`dx`), y-displacement (:math:`dy`), h-rise (:math:`h`),
inclination (:math:`\eta`), tip (:math:`\theta`) and twist (:math:`\Omega`), following elastic matrix is calculated.

.. math::

    \mathbf{K}_{helical-base-step} = \begin{bmatrix}
    K_{dx}        & K_{dx,dy}      & K_{dx,h}      & K_{dx,\eta}      & K_{dx,\theta}      & K_{dx,\Omega} \\
    K_{dx,dy}     & K_{dy}         & K_{dy,h}      & K_{dy,\eta}      & K_{dy,\theta}      & K_{dy,\Omega} \\
    K_{dx,h}      & K_{dy,h}       & K_{h}         & K_{h,\eta}       & K_{h,\theta}       & K_{h,\Omega} \\
    K_{dx,\eta}   & K_{dy,\eta}    & K_{h,\eta}    & K_{\eta}         & K_{\eta, \theta}   & K_{\eta,\Omega} \\
    K_{dx,\theta} & K_{dy,\theta}  & K_{h,\theta}  & K_{\eta, \theta} & K_{\theta}         & K_{\theta,\Omega} \\
    K_{dx,\Omega} & K_{dy,\Omega}  & K_{h,\Omega}  & K_{\eta, \Omega} & K_{\theta, \Omega} & K_{\Omega} \\
    \end{bmatrix}


Local Deformation Energy
------------------------

Using the elastic matrix, local deformation energy of base-steps can be calculated under the assumption that the energy
landscape is harmonic. Therefore, deformation energy is given as follows.

.. math::

    G = \frac{1}{2L_0}\mathbf{xKx^T}


**In case of local base-steps parameters:**

.. math::
    \mathbf{K} = \mathbf{K}_{base-step}

.. math::

    \mathbf{x} =  \begin{bmatrix}
                      (Dx_{i}-Dx_0)  &  (Dy_i - Dy_0) & (Dz_i - Dz_0) & (\tau_i - \tau_0) &
                      (\rho_i - \rho_0) & (\omega_i - \omega_0)
                  \end{bmatrix}


**In case of local helical base-steps parameters:**

.. math::
    \mathbf{K} = \mathbf{K}_{helical-base-step}

.. math::
    \mathbf{x} =  \begin{bmatrix}
                      (dx_{i}-dx_0)  &  (dy_i - dy_0) & (h_i - h_0) & (\eta_i - \eta_0) &
                      (\theta_i - \theta_0) & (\Omega_i - \Omega_0)
                  \end{bmatrix}
