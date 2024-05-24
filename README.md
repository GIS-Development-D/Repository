# UMEP Package for Improving Green Infrastructure

## Description

This package provides tools for urban microclimate modeling using the Urban Multi-scale Environmental Predictor (UMEP) within QGIS. It includes functions for preprocessing geospatial data, running urban climate models, and analyzing outputs. Key components integrated include SOLWEIG, Tree Planter, and URock, making it suitable for various urban climate and thermal comfort studies.

## Installation

### Prerequisites

1. **QGIS**: Install QGIS 3.36.1 on your system.
2. **Python**: Install Python 3.8.1 on your system.
3. **Virtual Environment**: Create a virtual environment to manage dependencies.

### Steps

1. **Install Dependencies**:
    ```bash
    pip install -r requirements.txt
    ```
2. **Update `config.json`**: Ensure UMEP(Our version 4.0.4) and UMEP Processing(Our version 2.0.13) plugins are installed in QGIS. Update the following paths:
    ```json
    {
        "qgis_prefix_path": "your_project_directory",
        "plugin_paths": [
            "your_qgis_path/apps/qgis/python/plugins",
            "C:/Users/your_user_name/AppData/Roaming/QGIS/QGIS3/profiles/default/python/plugins"
        ],
    }
    ```

## Functions

1. **cropDatasetWithMask(input_data, input_mask, input_crs, output_dir)**
    - Crops a dataset using a mask layer.
    - **Parameters**:
        - `input_data`: Path to the input raster dataset.
        - `input_mask`: Path to the mask layer.
        - `input_crs`: CRS of the input data.
        - `output_dir`: Directory to save the cropped output.
    - **Returns**: Clipping process result.

2. **skyViewFactor(input_cdsm, input_dsm, output_dir, output_file)**
    - Calculates the Sky View Factor (SVF).
    - **Parameters**:
        - `input_cdsm`: Path to the Canopy Digital Surface Model (CDSM).
        - `input_dsm`: Path to the Digital Surface Model (DSM).
        - `output_dir`: Directory to save the SVF output.
        - `output_file`: Name of the SVF output file.
    - **Returns**: SVF calculation result.

3. **wallHeightAspect(input_cdsm, output_height_file, output_aspect_file)**
    - Calculates wall height and aspect.
    - **Parameters**:
        - `input_cdsm`: Path to the CDSM.
        - `output_height_file`: Path to save wall height output.
        - `output_aspect_file`: Path to save wall aspect output.
    - **Returns**: Calculation result.

4. **run_SOLWEIG(input_dsm, input_cdsm, input_dem, input_meteo, svf_path, aniso_path, wallHeightRatioInputs, output_dir, \*\*kwargs)**
    - Runs the SOLWEIG model.
    - **Parameters**:
        - `input_dsm`, `input_cdsm`, `input_dem`, `input_meteo`: Paths to input files.
        - `svf_path`, `aniso_path`, `wallHeightRatioInputs`: Dictionaries for related inputs.
        - `output_dir`: Path to the output directory.
        - Additional SOLWEIG options as keyword arguments.
    - **Returns**: SOLWEIG processing output.

5. **treePlanter(input_SOLWEIG, input_polygonlayer, output_dir)**
    - Runs the Tree Planter model.
    - **Parameters**:
        - `input_SOLWEIG`: Path to SOLWEIG output directory.
        - `input_polygonlayer`: Path to polygon layer for planting areas.
        - `output_dir`: Directory to save Tree Planter output.
    - **Returns**: Tree Planter processing result.

6. **treeGenerator(input_pointlayer, input_dsm, input_dem, input_cdsm, output_cdsm, output_tdsm)**
    - Generates new trees in the CDSM.
    - **Parameters**:
        - Paths to input point layer, DSM, DEM, and CDSM.
        - Paths to save updated CDSM and TDSM.
    - **Returns**: None.

7. **transform_coordinate_vector_layer(input_file, output_file, target_crs)**
    - Transforms the coordinate system of a vector layer.
    - **Parameters**:
        - `input_file`: Path to the input vector file.
        - `output_file`: Path to the output vector file.
        - `target_crs`: EPSG code of the target CRS.
    - **Returns**: None.

8. **transform_coordinate_raster_layer(input_file, output_file, target_crs)**
    - Transforms the coordinate system of a raster layer.
    - **Parameters**:
        - `input_file`: Path to the input raster file.
        - `output_file`: Path to the output raster file.
        - `target_crs`: EPSG code of the target CRS.
    - **Returns**: None.

9. **URock_prepare(input_dsm, input_dem, input_cdsm, input_building, output_building, output_vegetation)**
    - Prepares building and vegetation data for URock.
    - **Parameters**:
        - `input_dsm`: Path to the Digital Surface Model (DSM).
        - `input_dem`: Path to the Digital Elevation Model (DEM).
        - `input_cdsm`: Path to the Canopy Digital Surface Model (CDSM).
        - `input_building`: Path to the building footprints.
        - `output_building`: Path to save the building height output.
        - `output_vegetation`: Path to save the vegetation height output.
    - **Returns**: Paths to the prepared building and vegetation data.

10. **URock(input_buildings, input_building_height, input_vegetation, input_vegetation_height, output_urock, output_filename, wind_height=1.5)**
    - Runs the URock model for urban wind field analysis.
    - **Parameters**:
        - `input_buildings`: Path to the building data.
        - `input_building_height`: Building height field name.
        - `input_vegetation`: Path to the vegetation data.
        - `input_vegetation_height`: Vegetation height field name.
        - `output_urock`: Directory to save URock outputs.
        - `output_filename`: Base name for URock output files.
        - `wind_height`: Wind height for analysis (default is 1.5).
    - **Returns**: None.

11. **regenerate_cdsm(input_directory="Data/TreePlanterTestData", output_TreePlanter_position='Output_TreePlanter_Position', input_TreePlanter_layer=os.path.join(output_TreePlanter_position, 'treePlanterPoint.shp'))**
    - Regenerates the Canopy Digital Surface Model (CDSM) using Tree Planter output.
    - **Parameters**:
        - `input_directory`: Directory containing input data.
        - `output_TreePlanter_position`: Directory containing Tree Planter output.
        - `input_TreePlanter_layer`: Path to the Tree Planter point layer.
    - **Returns**: None.

12. **SOLWEIG_Analysis(solweig_dir="Output_TreePlanter_SOLWEIG", stat_out="SOLWEIG_Analysis.tif", tmrt_stat_out="SOLWEIG_Analysis_Trmt.tif", stat_type=0)**
    - Analyzes SOLWEIG model output.
    - **Parameters**:
        - `solweig_dir`: Directory containing SOLWEIG output.
        - `stat_out`: Path to save the analysis output.
        - `tmrt_stat_out`: Path to save the mean radiant temperature (TMRT) analysis output.
        - `stat_type`: Type of statistic to calculate (default is 0).
    - **Returns**: Analysis result.

## Example Functions and Workflow

### SOLWEIG Example

1. **SOLWEIG_Example_Helsinki()**
    - Runs the SOLWEIG model for Helsinki data.

2. **SOLWEIG_Example_Goteborg()**
    - Runs the SOLWEIG model for Goteborg data.

### Green Infrastructure Solution Workflow

1. **treePlanter_Example()**

2. **SOLWEIG_Analysis()**
    - Runs initial SOLWEIG analysis and generates output files.

3. **regenerate_cdsm()**
    - Regenerates CDSM.

4. **treePlanter_Example()**
    - Runs Tree Planter with updated CDSM.

5. **SOLWEIG_Analysis()**
    - Runs improved SOLWEIG analysis and generates output files.

### URock Workflow

1. **URock_prepare()**
    - Prepares building and vegetation data for URock.

2. **transform_coordinate_vector_layer()**
    - Transforms coordinate systems for URock-prepared data.

3. **URock()**
    - Runs URock model.

4. **transform_coordinate_raster_layer()**
    - Transforms coordinate system for URock output raster.
  

## Web Service
