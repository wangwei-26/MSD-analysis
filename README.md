## Overview
A script designed to analyze sets of MSD data to determine the concentration of analyte in unknown samples in a 96-well format. Using raw signals from standards, it creats a standard curve with 4-parameter logistic (FourPL), and then uses
the standard curve to calculate concentrations for unknows and controls.

## Demo:
1. Prepare correct input files as below (please keep formats the same as in the [Example Files](example_files)):
    * raw signals (RP2_raw.txt)
    * plate design/layout (SP1_RP2_layout.txt)
    * concentrations of standards (Standards_conc_K15659U.txt)
    * a list of analytes in each well (Analytes_K15659U.txt)
    * sample list for output (samplelist_SP2_RP2.txt)
2. run MSD_4PL.py in [Code](Code)

3. When the provided example plate is run, it should produce three files, as seen in the [Expected Results](Expected_Results) folder.
