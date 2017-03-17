# BASS 0.2.2
- added print and summary methods for bass objects

# BASS 0.2.1
- added "small" parameter to allow for smaller memory footprint
- removed uses of timestamp to avoid a bug in Rstudio on Windows
- corrected example to use temperatures rather than inverse temperatures, and to match vignette

# BASS 0.2.0
- vignette added
- argument changes
    - temp.ladder now refers to temperatures rather than inverse temperatures
- Bug fixes:
    - Sobol index rounding for functional case
    - Handling responses in a dataframe
    - Error handling for w1 or w1 equal to 0
- Documentation updates
    - Changed input names a.beta.prec and b.beta.prec from the g-prior to a.tau and b.tau

# BASS 0.1.0
- First release.