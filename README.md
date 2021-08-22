# Downscaling-code
Downscaling code for TRMM and AMSR-E/2

The Tropical Rainfall Measuring Mission (TRMM) is a joint space mission between NASA and Japan's National Space Development Agency designed to monitor and study tropical and subtropical precipitation and the associated release of energy.

The Advanced Microwave Scanning Radiometer for EOS (AMSR-E) is a twelve-channel, six-frequency, total power passive-microwave radiometer system. It measures brightness temperatures at 6.925, 10.65, 18.7, 23.8, 36.5, and 89.0 GHz. Vertically and horizontally polarized measurements are taken at all channels. Measures precipitation rate, cloud water, water vapor, sea surface winds, sea surface temperature, ice, snow, and soil moisture.

Both the datasets are at 10km and 25km resolution. But for carrying out any scientific study on these datasets for a smaller region which may be a village or a point location it becomes difficult. Also when we use these datasets along any other datasets like NDVI from MODIS at 500m resolution or DEM at 30m resolution we need to resample these datasets which gives wrong results.

Hence to find a way around to use this datasets at a smaller resolution the given downscaling code was created.

The downsaling code for TRMM uses SRTM DEM and MODIS NDVI data. Whereas downscaling code for AMSR-E/2 uses MODIS NDVI data and MODIS LST datasets.

Using the Geographically Weighted Regression method both datasets were downscaled toa resolution of 500m.
