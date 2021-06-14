
FORMULATION 1: u-p formulation with Euler time integration

Governing equation are
1. Mixture moementum balance
2. Fluid momentum balance

Implement
1. non-incremental projection mehtod
2. incremental projection method

Stress and Strain are computed by using the B matrix computed at particle level (no enhancement such as BBar is used)

Damping is included, but needs improvements

Undrained boundary can be applied by making only the velicity in the direction of the ourward normal to the boundary to zero.

Lumped mass matrix is used in the formulation
