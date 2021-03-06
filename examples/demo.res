PSI Design and Analysis
Version: 0.0.1 
Design Codes: All Codes

Input File: demo.inp 

*** Preprocessing, initializing analysis...
*** Switching to base units.
*** Assembling system stiffness and force matrices.
*** Solving system equations for displacements.
*** Post processing elements...
*** Calculating support reactions and internal forces.
*** Writing loadcase results data.
*** Element code checking.
*** Code checking complete.
*** Writing code checking results data.
*** Analysis complete!



PSI 0.0.1, DATE: 2020-07-26 TIME: 08:16 PM
JOB NAME: None
LICENSED TO: PSI Community
MOVEMENTS : Displacement Report


MULTIPLE LOAD CASES


LOAD CASE DEFINITIONS

L1 (SUS) W[1] + P[1]
L2 (SUS) W[1] + P[1]
L3 (SUS) L1 + L2



                         TRANSLATIONS (inch)             ROTATIONS (deg)        

NODE      LOAD CASE         DX        DY        DZ        RX        RY        RZ

10         l1 (sus)      0.000    -0.000     0.000    0.0000    0.0000   -0.0000
           l2 (sus)      0.000    -0.000     0.000    0.0000    0.0000   -0.0000
           l3 (sus)      0.000     0.000     0.000    0.0000    0.0000    0.0000
                MAX      0.000     0.000     0.000    0.0000    0.0000    0.0000

20         l1 (sus)      0.000    -0.018     0.000    0.0000   -0.0000   -0.0002
           l2 (sus)      0.000    -0.018     0.000    0.0000   -0.0000   -0.0002
           l3 (sus)      0.000     0.026     0.000    0.0000    0.0000    0.0003
                MAX      0.000     0.026     0.000    0.0000    0.0000    0.0003





PSI 0.0.1, DATE: 2020-07-26 TIME: 08:16 PM
JOB NAME: None
LICENSED TO: PSI Community
FORCES : Element Forces Report


MULTIPLE LOAD CASES


LOAD CASE DEFINITIONS

L1 (SUS) W[1] + P[1]
L2 (SUS) W[1] + P[1]
L3 (SUS) L1 + L2



                             FORCES (lbf)               MOMENTS (foot*lbf)      

NODE      LOAD CASE         FX        FY        FZ        MX        MY        MZ

10         l1 (sus)        0.0     202.2       0.0       0.0       0.0    1685.0
           l2 (sus)        0.0     202.2       0.0       0.0       0.0    1685.0
           l3 (sus)        0.0     286.0       0.0       0.0       0.0    2383.0
                MAX        0.0     286.0       0.0       0.0       0.0    2383.0

20         l1 (sus)        0.0    -202.2       0.0       0.0       0.0     337.0
           l2 (sus)        0.0    -202.2       0.0       0.0       0.0     337.0
           l3 (sus)        0.0     286.0       0.0       0.0       0.0     476.6
                MAX        0.0     286.0       0.0       0.0       0.0     476.6





PSI 0.0.1, DATE: 2020-07-26 TIME: 08:16 PM
JOB NAME: None
LICENSED TO: PSI Community
REACTIONS : Support Reactions Report


MULTIPLE LOAD CASES


LOAD CASE DEFINITIONS

L1 (SUS) W[1] + P[1]
L2 (SUS) W[1] + P[1]
L3 (SUS) L1 + L2



                             FORCES (lbf)               MOMENTS (foot*lbf)      

NODE      LOAD CASE         FX        FY        FZ        MX        MY        MZ

10         l1 (sus)        0.0     404.4       0.0       0.0       0.0    2022.0
           l2 (sus)        0.0     404.4       0.0       0.0       0.0    2022.0
           l3 (sus)        0.0     571.9       0.0       0.0       0.0    2859.6
                MAX        0.0     571.9       0.0       0.0       0.0    2859.6

20         l1 (sus)        0.0       0.0       0.0       0.0       0.0       0.0
           l2 (sus)        0.0       0.0       0.0       0.0       0.0       0.0
           l3 (sus)        0.0       0.0       0.0       0.0       0.0       0.0
                MAX        0.0       0.0       0.0       0.0       0.0       0.0





PSI 0.0.1, DATE: 2020-07-26 TIME: 08:16 PM
JOB NAME: None
LICENSED TO: PSI Community
STRESSES : Element Stresses Report


MULTIPLE LOAD CASES


LOAD CASE DEFINITIONS

L1 (SUS) W[1] + P[1]
L2 (SUS) W[1] + P[1]



                    --------------- STRESSES (psi) -------------------

NODE      LOAD CASE       HOOP     AXIAL   TORSION  PRESSURE   BENDING      CODE

10         l1 (sus)     3310.9       0.0       0.0    1655.5     676.2     B31.1
           l2 (sus)     3310.9       0.0       0.0    1655.5     676.2
                MAX     3310.9       0.0       0.0    1655.5     676.2          

20         l1 (sus)     3310.9       0.0       0.0    1655.5     135.2     B31.1
           l2 (sus)     3310.9       0.0       0.0    1655.5     135.2
                MAX     3310.9       0.0       0.0    1655.5     135.2          





PSI 0.0.1, DATE: 2020-07-26 TIME: 08:16 PM
JOB NAME: None
LICENSED TO: PSI Community
CODES : Codes Report


MULTIPLE LOAD CASES


LOAD CASE DEFINITIONS

L1 (SUS) W[1] + P[1]
L2 (SUS) W[1] + P[1]



NODE      LOAD CASE     STRESS     ALLOW      SIFi      SIFo     RATIO      CODE
                           psi       psi

10         l1 (sus)     2331.6   30000.0      1.00      1.00      0.08     B31.1
           l2 (sus)     2331.6   30000.0      1.00      1.00      0.08
                MAX     2331.6   30000.0      1.00      1.00      0.08

20         l1 (sus)     1790.7   30000.0      1.00      1.00      0.06     B31.1
           l2 (sus)     1790.7   30000.0      1.00      1.00      0.06
                MAX     1790.7   30000.0      1.00      1.00      0.06





