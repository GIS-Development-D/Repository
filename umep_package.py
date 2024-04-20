import sys
import os
from qgis.core import QgsApplication
from qgis.core import QgsVectorLayer
from qgis.core import QgsField
from qgis.core import QgsVectorFileWriter
from qgis.PyQt.QtCore import QVariant

# ----------------- Initialization of the QGIS related dependencies -----------------
# Initiating a QGIS application
QgsApplication.setPrefixPath("D:/GIS-Development", True)
qgs = QgsApplication([], False)
qgs.initQgis()

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


# ----------------- UMEP Related Function Definition -----------------
def cropDatasetWithMask(function_input_data, function_input_mask, function_input_crs, function_output_dir):
    return processing.run("gdal:cliprasterbymasklayer",
                          {'INPUT': function_input_data,
                           'MASK': function_input_mask,
                           'SOURCE_CRS': function_input_crs,
                           'TARGET_CRS': function_input_crs,
                           'NODATA': None, 'ALPHA_BAND': False, 'CROP_TO_CUTLINE': True,
                           'KEEP_RESOLUTION': True, 'SET_RESOLUTION': False, 'X_RESOLUTION': None,
                           'Y_RESOLUTION': None, 'MULTITHREADING': False, 'OPTIONS': '',
                           'DATA_TYPE': 0, 'EXTRA': '', 'OUTPUT': function_output_dir})


def skyViewFactor(function_input_cdsm, function_input_dsm, function_output_dir, function_output_file):
    return processing.run("umep:Urban Geometry: Sky View Factor",
                          {'ANISO': True,
                           'INPUT_CDSM': function_input_cdsm,
                           'INPUT_DSM': function_input_dsm,
                           'INPUT_TDSM': None, 'INPUT_THEIGHT': 25,
                           'OUTPUT_DIR': function_output_dir,
                           'OUTPUT_FILE': function_output_file,
                           'TRANS_VEG': 3})


def wallHeightAspect(function_input_cdsm, function_output_height_file, function_output_aspect_file):
    return processing.run("umep:Urban Geometry: Wall Height and Aspect",
                          {'INPUT': function_input_cdsm,
                           'INPUT_LIMIT': 3,
                           'OUTPUT_HEIGHT': function_output_height_file,
                           'OUTPUT_ASPECT': function_output_aspect_file})


def run_solweig(input_dsm, input_cdsm, input_dem, input_meteo, svf_path, aniso_path, wallHeightRatioInputs, output_dir,
                trans_veg=3, t_height=25, albedo_walls=0.2, albedo_ground=0.15, emis_walls=0.9, emis_ground=0.95,
                abs_s=0.7, abs_l=0.95, posture=0, cyl=True, only_global=False, utc=0, poi_file=None,
                poi_field='', age=35, activity=80, clo=0.9, weight=75, height=180, sex=0, sensor_height=10):
    """
    Run the SOLWEIG model for outdoor thermal comfort.

    Parameters:
    - input_dsm: Path to the Digital Surface Model (DSM) input file.
    - input_cdsm: Path to the Canopy DSM input file.
    - input_dem: Path to the Digital Elevation Model (DEM) input file.
    - input_meteo: Path to the meteorological input file.
    - svf_path: Dictionary containing the paths to SVF-related inputs.
    - aniso_path: Dictionary containing the paths to aniso-related inputs.
    - wallHeightRatioOutputs: Dictionary containing the paths to wall height ratio outputs.
    - output_dir: Path to the output directory.
    - Additional parameters correspond to SOLWEIG model options.

    Returns:
    - Output of SOLWEIG processing run.
    """
    solweig_params = {
        'INPUT_DSM': input_dsm,
        'INPUT_SVF': svf_path,
        'INPUT_HEIGHT': wallHeightRatioInputs['OUTPUT_HEIGHT'],
        'INPUT_ASPECT': wallHeightRatioInputs['OUTPUT_ASPECT'],
        'INPUT_CDSM': input_cdsm,
        'TRANS_VEG': trans_veg,
        'INPUT_TDSM': None,
        'INPUT_THEIGHT': t_height,
        'INPUT_LC': None,
        'USE_LC_BUILD': False,
        'INPUT_DEM': input_dem,
        'SAVE_BUILD': False,
        'INPUT_ANISO': aniso_path,
        'ALBEDO_WALLS': albedo_walls,
        'ALBEDO_GROUND': albedo_ground,
        'EMIS_WALLS': emis_walls,
        'EMIS_GROUND': emis_ground,
        'ABS_S': abs_s,
        'ABS_L': abs_l,
        'POSTURE': posture,
        'CYL': cyl,
        'INPUTMET': input_meteo,
        'ONLYGLOBAL': only_global,
        'UTC': utc,
        'POI_FILE': poi_file,
        'POI_FIELD': poi_field,
        'AGE': age,
        'ACTIVITY': activity,
        'CLO': clo,
        'WEIGHT': weight,
        'HEIGHT': height,
        'SEX': sex,
        'SENSOR_HEIGHT': sensor_height,
        'OUTPUT_TMRT': True,
        'OUTPUT_KDOWN': False,
        'OUTPUT_KUP': False,
        'OUTPUT_LDOWN': False,
        'OUTPUT_LUP': False,
        'OUTPUT_SH': False,
        'OUTPUT_TREEPLANTER': True,
        'OUTPUT_DIR': output_dir
    }

    return processing.run("umep:Outdoor Thermal Comfort: SOLWEIG", solweig_params)


def treePlanter(function_input_SOLWEIG, function_input_polygonlayer, function_output_dir):
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


def treeGenerator(function_input_pointlayer, function_input_dsm, function_input_dem, function_input_cdsm,
                  function_output_cdsm, function_output_tdsm):
    processing.run("umep:Spatial Data: Tree Generator",
                   {
                       'INPUT_POINTLAYER': function_input_pointlayer,
                       'TREE_TYPE': 'ttype', 'TOT_HEIGHT': 'height', 'TRUNK_HEIGHT': 'trunk', 'DIA': 'diameter',
                       'INPUT_BUILD': None,
                       'INPUT_DSM': function_input_dsm,
                       'INPUT_DEM': function_input_dem,
                       'INPUT_CDSM': function_input_cdsm, 'INPUT_TDSM': None,
                       'CDSM_GRID_OUT': function_output_cdsm, 'TDSM_GRID_OUT': function_output_tdsm})


# This function is right now used to modify the attribute table of the generated trees
def attributeTableMod(path, field_name, field_value, field_type=QVariant.Int):
    layer = QgsVectorLayer(path)

    if not layer.isValid():
        print("Layer failed to load!")
        return

    print("Layer was loaded successfully!")

    # Check if the field already exists
    fields = layer.fields()
    if field_name in fields.names():
        print(f"Field '{field_name}' already exists. Skipping field addition and update.")
        return

    # Define the field type based on the input parameter
    field = QgsField(field_name, field_type)

    # Add the new field
    layer.dataProvider().addAttributes([field])
    layer.updateFields()

    # Start editing to update field values
    layer.startEditing()
    for feature in layer.getFeatures():
        feature.setAttribute(field_name, field_value)
        layer.updateFeature(feature)
    layer.commitChanges()
    print(f"Field '{field_name}' added and all features updated with the value '{field_value}'.")


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

    crop_cdsm = cropDatasetWithMask(os.path.join(input_directory, input_cdsm),
                                    os.path.join(input_directory, input_mask),
                                    cdsm_epsg, os.path.join(output_dir,
                                                            "Crop_" + \
                                                            input_cdsm_filename + ".tif"))

    input_dsm_filename = input_dsm.split(".")[0]
    crop_dsm = cropDatasetWithMask(os.path.join(input_directory, input_dsm), os.path.join(input_directory, input_mask),
                                   None, os.path.join(output_dir,
                                                      "Crop_" + \
                                                      input_dsm_filename + ".tif"))

    input_dem_filename = input_dem.split(".")[0]
    crop_dem = cropDatasetWithMask(os.path.join(input_directory, input_dem), os.path.join(input_directory, input_mask),
                                   None, os.path.join(output_dir,
                                                      "Crop_" + \
                                                      input_dem_filename + ".tif"))

    input_landcover_filename = input_landcover.split(".")[0]
    crop_landcover = cropDatasetWithMask(os.path.join(input_directory, input_landcover),
                                         os.path.join(input_directory, input_mask),
                                         None, os.path.join(output_dir,
                                                            "Crop_" + input_landcover_filename + ".tif"))

    # Calculates SVF from cropped data
    svf_outputs = skyViewFactor(crop_cdsm["OUTPUT"], crop_dsm["OUTPUT"], output_dir,
                                os.path.join(output_dir, 'SkyViewFactor.tif'))

    # Calculates wall height and wall aspect from cropped data
    wallHeightRatioOutputs = wallHeightAspect(crop_dsm["OUTPUT"], os.path.join(output_dir, 'wallHeight.tif'),
                                              os.path.join(output_dir, 'WallAspect.tif'))

    run_solweig(input_dsm=crop_dsm["OUTPUT"],
                input_cdsm=crop_cdsm["OUTPUT"],
                input_dem=crop_dsm["OUTPUT"],
                input_meteo=os.path.join(input_directory, input_meteo),
                svf_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'svfs.zip'),
                aniso_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'shadowmats.npz'),
                wallHeightRatioInputs={'OUTPUT_HEIGHT': wallHeightRatioOutputs['OUTPUT_HEIGHT'],
                                                 'OUTPUT_ASPECT': wallHeightRatioOutputs['OUTPUT_ASPECT']},
                output_dir=output_SOLWEIG_dir)


# ----------------- TreePlanter Function Example -----------------
def treePlanter_Example(input_cdsm="Data/TreePlanterTestData/CDSM.tif", tree_planter_flag=True):
    # Input files definition
    input_directory = "Data/TreePlanterTestData"
    input_dsm = os.path.join(input_directory, 'DSM.tif')
    input_dem = os.path.join(input_directory, 'DEM.tif')
    input_meteo = os.path.join(input_directory, 'metfile.txt')
    input_polygonlayer = os.path.join(input_directory, 'planting_area.shp')

    output_dir = "Output_temp"
    output_TreePlanter_SOLWEIG_dir = "Output_TreePlanter"
    output_TreePlanter_position = 'Output_TreePlanter_Position'

    # Calculates SVF from cropped data
    svf_outputs = skyViewFactor(input_cdsm, input_dsm, output_dir,
                                os.path.join(output_dir, 'SkyViewFactor.tif'))

    # Calculates wall height and wall aspect from cropped data
    wallHeightRatioOutputs = wallHeightAspect(input_dsm, os.path.join(output_dir, 'wallHeight.tif'),
                                              os.path.join(output_dir, 'WallAspect.tif'))

    run_solweig(input_dsm=input_dsm,
                input_cdsm=input_cdsm,
                input_dem=input_dem,
                input_meteo=input_meteo,
                svf_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'svfs.zip'),
                aniso_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'shadowmats.npz'),
                wallHeightRatioInputs={'OUTPUT_HEIGHT': wallHeightRatioOutputs['OUTPUT_HEIGHT'],
                                                 'OUTPUT_ASPECT': wallHeightRatioOutputs['OUTPUT_ASPECT']},
                output_dir=output_TreePlanter_SOLWEIG_dir)

    if tree_planter_flag:
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


# ----------------- Update the cdsm with New Trees Example -----------------
def regenerate_cdsm():
    input_directory = "Data/TreePlanterTestData"
    input_cdsm = os.path.join(input_directory, 'CDSM.tif')
    input_dsm = os.path.join(input_directory, 'DSM.tif')
    input_dem = os.path.join(input_directory, 'DEM.tif')

    output_TreePlanter_position = 'Output_TreePlanter_Position'
    input_TreePlanter_layer = os.path.join(output_TreePlanter_position, 'treePlanterPoint.shp')
    output_cdsm = "Output_TreeGenerator/cdsm.tif"
    output_tdsm = "Output_TreeGenerator/tdsm.tif"

    attributeTableMod(input_TreePlanter_layer, "ttype", 1)
    treeGenerator(input_TreePlanter_layer, input_dsm, input_dem, input_cdsm, output_cdsm, output_tdsm)


# SOLWEIG_Example_Goteborg()

# treePlanter_Example()
# regenerate_cdsm()
# treePlanter_Example(input_cdsm="Output_TreeGenerator/cdsm.tif", tree_planter_flag=False)
