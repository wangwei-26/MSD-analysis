# MSD analysis
A script designed to analyze sets of MSD data to determine the concentration of analyte in unknown samples. 
Using raw signals from standards, it creats a standard curve with 4-parameter logistic (FourPL), and then uses
the standard curve to calculate concentrations for unknows and controls.

# How to run:
1. files to prepare (please keep formats the same as in the example files)
    a. raw signals (RP2_raw.txt)
    b. plate design/layout (SP1_RP2_layout.txt)
    c. concentrations of standards (Standards_conc_K15659U.txt)
    d. a list of analytes in each well (Analytes_K15659U.txt)
    c. sample list for output (samplelist_SP2_RP2.txt)
2. run MSD_4PL.py
