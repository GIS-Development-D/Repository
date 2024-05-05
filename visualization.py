from osgeo import gdal
import matplotlib.pyplot as plt

file_path1 = 'SOLWEIG_Analysis_Initial.tif'
file_path2 = 'SOLWEIG_Analysis_Improved.tif'

ds1 = gdal.Open(file_path1)
if ds1 is None:
    print("Error: The data is not exist")
    exit()

band1 = ds1.GetRasterBand(1)
data1 = band1.ReadAsArray()

ds2 = gdal.Open(file_path2)
if ds2 is None:
    print("Error: The data is not exist")
    exit()

band2 = ds2.GetRasterBand(1)
data2 = band2.ReadAsArray()

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

cax1 = ax1.imshow(data1, cmap='gray')
fig.colorbar(cax1, ax=ax1, orientation='vertical')
ax1.set_title('Before planting trees')

cax2 = ax2.imshow(data2, cmap='gray')
fig.colorbar(cax2, ax=ax2, orientation='vertical')
ax2.set_title('After planting trees')

plt.show()
