import sys
import os
from tqdm import tqdm
from qgis.core import QgsApplication
from qgis.core import QgsCoordinateTransform
from qgis.core import QgsVectorLayer
from qgis.core import QgsRasterLayer
from qgis.core import QgsField
from qgis.core import QgsProject
from qgis.core import QgsRasterFileWriter
from qgis.core import QgsRasterPipe
from qgis.core import QgsRasterProjector
from osgeo import ogr, osr

from PyQt5.QtWidgets import QApplication
from qgis.PyQt.QtCore import QVariant

# ----------------- Initialization of the QGIS related dependencies -----------------
# Initiating a QGIS application
QgsApplication.setPrefixPath("D:/GIS-Development", True)
qgs = QgsApplication([], False)
qgs.initQgis()

# Create the application instance for the GUI
app = QApplication([])

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


def run_SOLWEIG(input_dsm, input_cdsm, input_dem, input_meteo, svf_path, aniso_path, wallHeightRatioInputs, output_dir,
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
    output_dir = "Output_temp_SOLWEIG"
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

    run_SOLWEIG(input_dsm=crop_dsm["OUTPUT"],
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

    output_dir = "Output_temp_TreePlanter"
    output_TreePlanter_SOLWEIG_dir = "Output_TreePlanter_SOLWEIG"
    output_TreePlanter_position = 'Output_TreePlanter_Position'

    pbar = tqdm(total=4)

    # Calculates SVF from cropped data
    svf_outputs = skyViewFactor(input_cdsm, input_dsm, output_dir,
                                os.path.join(output_dir, 'SkyViewFactor.tif'))
    pbar.update(1)

    # Calculates wall height and wall aspect from cropped data
    wallHeightRatioOutputs = wallHeightAspect(input_dsm, os.path.join(output_dir, 'wallHeight.tif'),
                                              os.path.join(output_dir, 'WallAspect.tif'))
    pbar.update(1)

    run_SOLWEIG(input_dsm=input_dsm,
                input_cdsm=input_cdsm,
                input_dem=input_dem,
                input_meteo=input_meteo,
                svf_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'svfs.zip'),
                aniso_path=os.path.join(svf_outputs['OUTPUT_DIR'], 'shadowmats.npz'),
                wallHeightRatioInputs={'OUTPUT_HEIGHT': wallHeightRatioOutputs['OUTPUT_HEIGHT'],
                                       'OUTPUT_ASPECT': wallHeightRatioOutputs['OUTPUT_ASPECT']},
                output_dir=output_TreePlanter_SOLWEIG_dir)
    pbar.update(1)

    if tree_planter_flag:
        processing.run("umep:Outdoor Thermal Comfort: TreePlanter",
                       {'SOLWEIG_DIR': output_TreePlanter_SOLWEIG_dir,
                        'INPUT_POLYGONLAYER': input_polygonlayer,
                        'START_HOUR': 8,
                        'END_HOUR': 18,
                        'TTYPE': 0,
                        'HEIGHT': 25,
                        'DIA': 1.5,
                        'TRUNK': 15,
                        'TRANS_VEG': 30,
                        'NTREE': 50,
                        # Problem existed in the OUTPUT_CDSM, we use regenerate_cdsm instead
                        'OUTPUT_CDSM': True, 'OUTPUT_POINTFILE': True, 'OUTPUT_OCCURRENCE': False,
                        'OUTPUT_DIR': output_TreePlanter_position, 'ITERATIONS': 2000,
                        'INCLUDE_OUTSIDE': True, 'RANDOM_STARTING': False, 'GREEDY_ALGORITHM': True})
        pbar.update(1)
    pbar.close()


# ----------------- URock preparation Example -----------------
def URock_prepare(input_dsm, input_dem, input_cdsm, input_building, output_building, output_vegetation):
    processing.run("umep:Urban Wind Field: URock Prepare",
                   {'INPUT_BUILD_FOOTPRINT': input_building,
                    'INPUT_BUILD_DSM': input_dsm,
                    'INPUT_BUILD_DEM': input_dem,
                    'INPUT_VEG_CDSM': input_cdsm,
                    'INPUT_VEG_POINTS': None,
                    'HEIGHT_VEG_FIELD': '',
                    'RADIUS_VEG_FIELD': '',
                    'VEGETATION_ASPECT': '0.75',
                    'BUILDINGS_WITH_HEIGHT': output_building,
                    'OUTPUT_BUILD_HEIGHT_FIELD': 'ROOF_HEIGHT',
                    'VEGETATION_WITH_HEIGHT': output_vegetation,
                    'OUTPUT_VEG_HEIGHT_FIELD': 'VEG_HEIGHT'})

    return output_building, output_vegetation


# ----------------- URock analysis Example -----------------
def URock(input_buildings, input_building_height, input_vegetation, input_vegetation_height, output_urock,
          output_filename, wind_height=1.5):
    processing.run("umep:Urban Wind Field: URock", {'BUILDINGS': input_buildings,
                                                    'HEIGHT_FIELD_BUILD': input_building_height,
                                                    'VEGETATION': input_vegetation,
                                                    'VEGETATION_CROWN_TOP_HEIGHT': input_vegetation_height,
                                                    'VEGETATION_CROWN_BASE_HEIGHT': '', 'ATTENUATION_FIELD': '',
                                                    'INPUT_PROFILE_FILE': '', 'INPUT_PROFILE_TYPE': 0,
                                                    'INPUT_WIND_HEIGHT': 10, 'INPUT_WIND_SPEED': 2,
                                                    'INPUT_WIND_DIRECTION': 45, 'RASTER_OUTPUT': None,
                                                    'HORIZONTAL_RESOLUTION': 2, 'VERTICAL_RESOLUTION': 2,
                                                    'WIND_HEIGHT': wind_height,
                                                    'UROCK_OUTPUT': output_urock,
                                                    'OUTPUT_FILENAME': output_filename, 'SAVE_RASTER': True,
                                                    'SAVE_VECTOR': True, 'SAVE_NETCDF': False, 'LOAD_OUTPUT': False})


def transform_coordinate_vector_layer(input_file, output_file, target_crs=3006):
    target_proj = osr.SpatialReference()
    target_proj.ImportFromEPSG(target_crs)

    ds = ogr.Open(input_file, 0)  # 0 means read-only mode
    layer = ds.GetLayer(0)

    driver = ogr.GetDriverByName('GPKG')
    out_ds = driver.CreateDataSource(output_file)  # 指定新文件的名字
    out_layer = out_ds.CreateLayer(layer.GetName(), target_proj, layer.GetGeomType())

    layer_defn = layer.GetLayerDefn()
    for i in range(layer_defn.GetFieldCount()):
        field_defn = layer_defn.GetFieldDefn(i)
        out_layer.CreateField(field_defn)

    transform = osr.CoordinateTransformation(layer.GetSpatialRef(), target_proj)
    for feature in layer:
        geom = feature.GetGeometryRef()
        geom.Transform(transform)

        new_feature = ogr.Feature(out_layer.GetLayerDefn())
        new_feature.SetGeometry(geom)
        for i in range(layer_defn.GetFieldCount()):
            new_feature.SetField(layer_defn.GetFieldDefn(i).GetNameRef(), feature.GetField(i))
        out_layer.CreateFeature(new_feature)

    del ds
    del out_ds


def transform_coordinate_raster_layer(input_file, output_file, target_crs='EPSG:3006'):

    raster_layer = QgsRasterLayer(input_file, 'Original Raster')
    if not raster_layer.isValid():
        print("Loaded layer that would be transformed failed!")
    else:
        print("Loaded layer that would be transformed successfully!")
        source_crs = raster_layer.crs()
        target_crs = QgsCoordinateReferenceSystem(target_crs)

        transform_context = QgsProject.instance().transformContext()
        coord_transform = QgsCoordinateTransform(source_crs, target_crs, transform_context)

        raster_projector = QgsRasterProjector()
        raster_projector.setCrs(source_crs, target_crs, transform_context)

        pipe = QgsRasterPipe()
        pipe.set(raster_layer.dataProvider().clone())
        pipe.insert(2, raster_projector)

        writer = QgsRasterFileWriter(output_file)
        error = writer.writeRaster(
            pipe,
            raster_layer.width(),
            raster_layer.height(),
            raster_layer.extent(),
            target_crs
        )

        if error == QgsRasterFileWriter.NoError:
            print("Raster file reprojected and saved successfully.")
        else:
            print("Failed to reproject raster file: Error code", error)


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

    # Modify the attribute table
    attributeTableMod(input_TreePlanter_layer, "ttype", 1)

    # treeGenerator is used to generate the new cdsm
    treeGenerator(input_TreePlanter_layer, input_dsm, input_dem, input_cdsm, output_cdsm, output_tdsm)


def SOLWEIG_Analysis(solweig_dir="Output_TreePlanter_SOLWEIG", stat_out="SOLWEIG_Analysis.tif",
                     tmrt_stat_out="SOLWEIG_Analysis_Trmt.tif", stat_type=0):
    return processing.run("umep:Outdoor Thermal Comfort: SOLWEIG Analyzer",
                          {'SOLWEIG_DIR': solweig_dir, 'BUILDINGS': None, 'VARIA_IN': 0,
                           'STAT_TYPE': stat_type, 'THRES_TYPE': 1, 'TMRT_THRES_NUM': 50,
                           'STAT_OUT': stat_out, 'TMRT_STAT_OUT': tmrt_stat_out})


# SOLWEIG_Example_Goteborg()

# ----------------- The workflow of the green infrastructure solution -----------------
# treePlanter_Example()
# SOLWEIG_Analysis(stat_out="SOLWEIG_Analysis_Initial.tif",
#                  tmrt_stat_out="SOLWEIG_Analysis_Initial_tmrt.tif")
# regenerate_cdsm()
# treePlanter_Example(input_cdsm="Output_TreeGenerator/cdsm.tif", tree_planter_flag=False)
# SOLWEIG_Analysis(stat_out="SOLWEIG_Analysis_Improved.tif",
#                  tmrt_stat_out="SOLWEIG_Analysis_Improved_tmrt.tif")

# ----------------- The workflow of the URock(wind speed) analysis -----------------
# prepared_building, prepared_vegetation = URock_prepare(input_dsm='URock/Annedal_EPSG3006/dsm.tif',
#                                                        input_dem='URock/Annedal_EPSG3006/dem.tif',
#                                                        input_cdsm='URock/Annedal_EPSG3006/cdsm.tif',
#                                                        input_building='URock/Annedal_EPSG3006/buildings.shp',
#                                                        output_building='URock/Output/building_urock.gpkg',
#                                                        output_vegetation='URock/Output/veg_urock.gpkg')
# transform_coordinate_vector_layer(input_file=prepared_building, output_file="URock/Output/building_urock_3006.gpkg")
# transform_coordinate_vector_layer(input_file=prepared_vegetation, output_file="URock/Output/veg_urock_3006.gpkg")
#
# URock(input_buildings="URock/Output/building_urock_3006.gpkg",
#       input_building_height='ROOF_HEIGHT',
#       input_vegetation="URock/Output/veg_urock_3006.gpkg",
#       input_vegetation_height='VEG_HEIGHT',
#       output_urock='URock/Output',
#       output_filename='urock_output',
#       wind_height=1.5)
#
# transform_coordinate_raster_layer(input_file="URock/Output/z1_5/urock_outputWS.Gtiff",
#                                   output_file="URock/Output/z1_5/urock_outputWS_3006.Gtiff",
#                                   target_crs='EPSG:3006')
