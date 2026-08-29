Write a subroutine to compare max refl and OLR from WRF vs HRRR analysis at a given time that satisfies the following requirements:

    1. Use the same Lambert map projection when plotting (use cartopy or whatever that is the best)
    2. The WRF domain is smaller than that of HRRR. Make sure the plots are the same size and cover the same area for better comparison. The corners in the WRF plot will have not data but that is okay. 
    3. Use pcolormesh and allow me to adjust the contour intervals 
    4. put on state borders and long/lat lines
    