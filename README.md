# Python Script for Plotting WRF Precipitation

### Script Overview

This Python script downloads meteorological data from the **WRF model** via the **THREDDS server of MeteoGalicia**, processes it, and generates visual representations of sea-level pressure and precipitation for the current day. 

### Key Features

1. **Data Download and Loading**:  
   The script automatically downloads a NetCDF file from the **MeteoGalicia THREDDS server** containing WRF model output for the current date. If the file already exists locally, it opens it instead of downloading it again.

2. **Georeferencing**:  
   The function `georeferencio()` computes the latitude and longitude grid using parameters from the NetCDF file, converting them into a map projection with the `Basemap` library.

3. **Daily Rain Calculation**:  
   The `dailyrain()` function calculates the accumulated rain over the last three hours by subtracting rain data at the current and previous time steps.

4. **Sea Level Pressure (SLP) Calculation**:  
   The `pressure()` function computes sea-level pressure by adjusting the surface pressure data with elevation and temperature corrections.

5. **Plotting Meteorological Data**:  
   For each time step in the NetCDF file, the script generates a plot with the following:
   - **Sea-Level Pressure Contours**: These are plotted as isobars on a map.
   - **Highs and Lows**: The script identifies and labels areas of high and low pressure on the map.
   - **Optional Rain Visualization**: You can uncomment the corresponding section to display rain data on the map.

6. **Saving Plots**:  
   Each plot is saved as a PNG image file. After all plots are created, the script compiles them into an animated GIF showing the progression of the weather over time.

7. **Cleaning Up**:  
   The script deletes the raw NetCDF and optional PNG, keeping only the final GIF as the output.

