The paper "A simple and effective method for simulating nested exchangeable correlated binary data for longitudinal cluster randomised trials" by Rhys A. Bowden, Jessica Kasza and Andrew B. Forbes accompanies this package and explains the underlying mathematical theory. However, some of the notation in the paper differs from that in the arguments of the package. This file aims to explain those differences.

Concept                       Paper notation       rNestBin argument name
------------------------------------------------------------------------------
Desired prevalences for       pi_1, ..., pi_T      means
each period               

Within-period correlation     rho_C                rhoC

Between-period correlation    rho_{CT}             rhoCT

Observations in period j      N_j                  n (can be scalar or vector)

Number of clusters to sample  C                    C

Number of subclusters per     T                    length of means vector
cluster