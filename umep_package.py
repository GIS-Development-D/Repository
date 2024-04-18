import sys
import os
from qgis.core import QgsApplication

# ----------------- Initialization of the QGIS related dependencies -----------------
# Initiating a QGIS application
home = 'D:/GIS-Development'
QgsApplication.setPrefixPath(home, True)
app = QgsApplication([], False)
app.initQgis()

# This line import the default plugins of the QGIS
sys.path.append(r'D:\QGIS 3.36.1\apps\qgis\python\plugins')
import processing
from processing.core.Processing import Processing

# It is important to initialize a Processing before using it
Processing.initialize()
# Remember to change it as your own local path of the QGIS plugins(This is for UMEP tools)
sys.path.append(r'C:\Users\ASUS\AppData\Roaming\QGIS\QGIS3\profiles\default\python\plugins')
from processing_umep.processing_umep_provider import ProcessingUMEPProvider
from qgis.core import QgsCoordinateReferenceSystem

umep_provider = ProcessingUMEPProvider()
QgsApplication.processingRegistry().addProvider(umep_provider)


# ----------------- UMEP Function Definition -----------------
def CropDatasetWithMask(function_input_data, function_input_mask, function_input_crs, function_output_dir):
    return processing.run("gdal:cliprasterbymasklayer",
                          {'INPUT': function_input_data,
                           'MASK': function_input_mask,
                           'SOURCE_CRS': function_input_crs,
                           'TARGET_CRS': function_input_crs,
                           'NODATA': None, 'ALPHA_BAND': False, 'CROP_TO_CUTLINE': True,
                           'KEEP_RESOLUTION': True, 'SET_RESOLUTION': False, 'X_RESOLUTION': None,
                           'Y_RESOLUTION': None, 'MULTITHREADING': False, 'OPTIONS': '',
                           'DATA_TYPE': 0, 'EXTRA': '', 'OUTPUT': function_output_dir})


def SkyViewFactor(function_input_cdsm, function_input_dsm, function_output_dir, function_output_file):
    return processing.run("umep:Urban Geometry: Sky View Factor",
                          {'ANISO': True,
                           'INPUT_CDSM': function_input_cdsm,
                           'INPUT_DSM': function_input_dsm,
                           'INPUT_TDSM': None, 'INPUT_THEIGHT': 25,
                           'OUTPUT_DIR': function_output_dir,
                           'OUTPUT_FILE': function_output_file,
                           'TRANS_VEG': 3})


def WallHeightAspect(function_input_cdsm, function_output_height_file, function_output_aspect_file):
    return processing.run("umep:Urban Geometry: Wall Height and Aspect",
                          {'INPUT': function_input_cdsm,
                           'INPUT_LIMIT': 3,
                           'OUTPUT_HEIGHT': function_output_height_file,
                           'OUTPUT_ASPECT': function_output_aspect_file})


def TreePlanter(function_input_SOLWEIG, function_input_polygonlayer, function_output_dir):
    return processing.run("umep:Outdoor Thermal Comfort: TreePlanter",
                          {'SOLWEIG_DIR': function_input_SOLWEIG,
                           'INPUT_POLYGONLAYER': function_input_polygonlayer,
                           'START_HOUR': 13,
                           'END_HOUR': 15,
                           'TTYPE': 0,
                           'HEIGHT': 10,
                           'DIA': 5,
                           'TRUNK': 3,
                           'TRANS_VEG': 3,
                           'NTREE': 3,
                           'OUTPUT_CDSM': True, 'OUTPUT_POINTFILE': True, 'OUTPUT_OCCURRENCE': False,
                           'OUTPUT_DIR': function_output_dir, 'ITERATIONS': 2000,
                           'INCLUDE_OUTSIDE': True, 'RANDOM_STARTING': False, 'GREEDY_ALGORITHM': False})


# ----------------- SOLWEIG Function Example -----------------
def SOLWEIG_Example_Goteborg():
    # Input files definition
    input_directory = "Data/Goteborg_SWEREF99_1200"
    input_mask = "mask_layer.geojson"
    input_cdsm = 'CDSM_KRbig.asc'
    input_dsm = 'DSM_KRbig.tif'
    input_dem = 'DEM_KRbig.tif'
    input_landcover = 'landcover.tif'
    input_meteo = 'gbg19970606_2015a.txt'
    cdsm_epsg = QgsCoordinateReferenceSystem('EPSG:3007')
    input_cdsm_filename = input_cdsm.split(".")[0]

    # Defines an output directory where will be stored your outputs (and intermediate results)
    output_dir = "Output_temp"
    output_SOLWEIG_dir = "Output_SOLWEIG"

    crop_cdsm = CropDatasetWithMask(os.path.join(input_directory, input_cdsm),
                                    os.path.join(input_directory, input_mask),
                                    cdsm_epsg, os.path.join(output_dir,
                                                            "Crop_" + \
                                                            input_cdsm_filename + ".tif"))

    input_dsm_filename = input_dsm.split(".")[0]
    crop_dsm = CropDatasetWithMask(os.path.join(input_directory, input_dsm), os.path.join(input_directory, input_mask),
                                   None, os.path.join(output_dir,
                                                      "Crop_" + \
                                                      input_dsm_filename + ".tif"))

    input_dem_filename = input_dem.split(".")[0]
    crop_dem = CropDatasetWithMask(os.path.join(input_directory, input_dem), os.path.join(input_directory, input_mask),
                                   None, os.path.join(output_dir,
                                                      "Crop_" + \
                                                      input_dem_filename + ".tif"))

    input_landcover_filename = input_landcover.split(".")[0]
    crop_landcover = CropDatasetWithMask(os.path.join(input_directory, input_landcover),
                                         os.path.join(input_directory, input_mask), None, os.path.join(output_dir,
                                                                                                       "Crop_" + \
                                                                                                       input_landcover_filename + ".tif"))

    # Calculates SVF from cropped data
    svf_outputs = SkyViewFactor(crop_cdsm["OUTPUT"], crop_dsm["OUTPUT"], output_dir,
                                os.path.join(output_dir, 'SkyViewFactor.tif'))

    # Calculates wall height and wall aspect from cropped data
    wallHeightRatioOutputs = WallHeightAspect(crop_dsm["OUTPUT"], os.path.join(output_dir, 'wallHeight.tif'),
                                              os.path.join(output_dir, 'WallAspect.tif'))

    processing.run("umep:Outdoor Thermal Comfort: SOLWEIG",
                   {'INPUT_DSM': crop_dsm["OUTPUT"],
                    'INPUT_SVF': os.path.join(svf_outputs['OUTPUT_DIR'], 'svfs.zip'),
                    'INPUT_HEIGHT': wallHeightRatioOutputs['OUTPUT_HEIGHT'],
                    'INPUT_ASPECT': wallHeightRatioOutputs['OUTPUT_ASPECT'],
                    'INPUT_CDSM': crop_cdsm["OUTPUT"],
                    'TRANS_VEG': 3, 'INPUT_TDSM': None, 'INPUT_THEIGHT': 25,
                    'INPUT_LC': crop_landcover["OUTPUT"],
                    'USE_LC_BUILD': False,
                    'INPUT_DEM': crop_dsm["OUTPUT"],
                    'SAVE_BUILD': False,
                    'INPUT_ANISO': os.path.join(svf_outputs['OUTPUT_DIR'], 'shadowmats.npz'),
                    'ALBEDO_WALLS': 0.2, 'ALBEDO_GROUND': 0.15, 'EMIS_WALLS': 0.9, 'EMIS_GROUND': 0.95,
                    'ABS_S': 0.7, 'ABS_L': 0.95, 'POSTURE': 0, 'CYL': True,
                    'INPUTMET': os.path.join(input_directory, input_meteo),
                    'ONLYGLOBAL': False, 'UTC': 0, 'POI_FILE': None, 'POI_FIELD': '', 'AGE': 35,
                    'ACTIVITY': 80, 'CLO': 0.9, 'WEIGHT': 75, 'HEIGHT': 180, 'SEX': 0, 'SENSOR_HEIGHT': 10,
                    'OUTPUT_TMRT': True, 'OUTPUT_KDOWN': False, 'OUTPUT_KUP': False, 'OUTPUT_LDOWN': False,
                    'OUTPUT_LUP': False, 'OUTPUT_SH': False, 'OUTPUT_TREEPLANTER': False,
                    'OUTPUT_DIR': output_SOLWEIG_dir})


def TreePlanter_Example():
    # Input files definition
    input_directory = "Data/TreePlanterTestData"
    input_cdsm = os.path.join(input_directory, 'CDSM.tif')
    input_dsm = os.path.join(input_directory, 'DSM.tif')
    input_dem = os.path.join(input_directory, 'DEM.tif')
    input_meteo = os.path.join(input_directory, 'metfile.txt')
    input_polygonlayer = os.path.join(input_directory, 'planting_area.shp')

    output_dir = "Output_temp"
    output_TreePlanter_SOLWEIG_dir = "Output_TreePlanter"
    output_TreePlanter_position = os.path.join(output_TreePlanter_SOLWEIG_dir, 'Output_TreePlanter_Position')

    # Calculates SVF from cropped data
    svf_outputs = SkyViewFactor(input_cdsm, input_dsm, output_dir,
                                os.path.join(output_dir, 'SkyViewFactor.tif'))

    # Calculates wall height and wall aspect from cropped data
    wallHeightRatioOutputs = WallHeightAspect(input_dsm, os.path.join(output_dir, 'wallHeight.tif'),
                                              os.path.join(output_dir, 'WallAspect.tif'))

    processing.run("umep:Outdoor Thermal Comfort: SOLWEIG",
                                     {'INPUT_DSM': input_dsm,
                                      'INPUT_SVF': os.path.join(svf_outputs['OUTPUT_DIR'], 'svfs.zip'),
                                      'INPUT_HEIGHT': wallHeightRatioOutputs['OUTPUT_HEIGHT'],
                                      'INPUT_ASPECT': wallHeightRatioOutputs['OUTPUT_ASPECT'],
                                      'INPUT_CDSM': input_cdsm,
                                      'TRANS_VEG': 3, 'INPUT_TDSM': None, 'INPUT_THEIGHT': 25,
                                      'INPUT_LC': None,
                                      'USE_LC_BUILD': False,
                                      'INPUT_DEM': input_dem,
                                      'SAVE_BUILD': False,
                                      'INPUT_ANISO': os.path.join(svf_outputs['OUTPUT_DIR'], 'shadowmats.npz'),
                                      'ALBEDO_WALLS': 0.2, 'ALBEDO_GROUND': 0.15, 'EMIS_WALLS': 0.9,
                                      'EMIS_GROUND': 0.95,
                                      'ABS_S': 0.7, 'ABS_L': 0.95, 'POSTURE': 0, 'CYL': True,
                                      'INPUTMET': input_meteo,
                                      'ONLYGLOBAL': False, 'UTC': 0, 'POI_FILE': None, 'POI_FIELD': '', 'AGE': 35,
                                      'ACTIVITY': 80, 'CLO': 0.9, 'WEIGHT': 75, 'HEIGHT': 180, 'SEX': 0,
                                      'SENSOR_HEIGHT': 10,
                                      'OUTPUT_TMRT': True, 'OUTPUT_KDOWN': False, 'OUTPUT_KUP': False,
                                      'OUTPUT_LDOWN': False,
                                      'OUTPUT_LUP': False, 'OUTPUT_SH': False, 'OUTPUT_TREEPLANTER': True,
                                      'OUTPUT_DIR': output_TreePlanter_SOLWEIG_dir})

    processing.run("umep:Outdoor Thermal Comfort: TreePlanter",
                                  {'SOLWEIG_DIR': output_TreePlanter_SOLWEIG_dir,
                                   'INPUT_POLYGONLAYER': input_polygonlayer,
                                   'START_HOUR': 13,
                                   'END_HOUR': 15,
                                   'TTYPE': 0,
                                   'HEIGHT': 10,
                                   'DIA': 5,
                                   'TRUNK': 3,
                                   'TRANS_VEG': 3,
                                   'NTREE': 3,
                                   'OUTPUT_CDSM': True, 'OUTPUT_POINTFILE': True, 'OUTPUT_OCCURRENCE': False,
                                   'OUTPUT_DIR': output_TreePlanter_position, 'ITERATIONS': 2000,
                                   'INCLUDE_OUTSIDE': True, 'RANDOM_STARTING': False, 'GREEDY_ALGORITHM': False})


# SOLWEIG_Example_Goteborg()
TreePlanter_Example()
