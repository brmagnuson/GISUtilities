import sys
import ogr
import os
import osr
import pandas


def copyLayerFieldDefinitions(fromLayer, toLayer):
    """
    Copy field definitions (not data) from one layer to another
    :param fromLayer: layer object that contains the field definitions to copy
    :param toLayer: layer object to copy the field definitions into
    :return
    """
    featureDefinition = fromLayer.GetLayerDefn()
    for i in range(featureDefinition.GetFieldCount()):
        fieldDefinition = featureDefinition.GetFieldDefn(i)
        toLayer.CreateField(fieldDefinition)

    return


def copyFeatureFieldValues(fromFeature, toFeature, isShapefile=True):
    """
    Copy field values from one feature to another. Note that this assumes that the features have the same fields!
    :param fromFeature: feature object that contains the data to copy
    :param toFeature: feature object that the data is to be copied into
    :param isShapefile: True if toFeature is from a shapefile, False otherwise
    :return
    """
    for i in range(fromFeature.GetFieldCount()):
        fieldName = fromFeature.GetFieldDefnRef(i).GetName()

        # Set field as long as its value is not None
        if fromFeature.GetField(fieldName) is not None:

            # Deal with shapefile field name length limitation of 10 characters
            if isShapefile:
                if len(fieldName) > 10:
                    toFeature.SetField(fieldName[:10], fromFeature.GetField(fieldName))
                else:
                    toFeature.SetField(fieldName, fromFeature.GetField(fieldName))
            else:
                toFeature.SetField(fieldName, fromFeature.GetField(fieldName))

    return


def reprojectSpatialFile(inputFilePath, outputFilePath, driverName, outputEPSG):
    """
    Reprojects a spatial file from one projected coordinate system to another based on EPSG number.
    :param inputFilePath: String. Path to existing spatial file
    :param outputFilePath: String. Path to new spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param outputEPSG: Int. EPSG code for new projection
    :return:
    """

    # Get geometry from input file
    inputLayer, inputDataSource = getLayer(inputFilePath, driverName)
    geomType = inputLayer.GetLayerDefn().GetGeomType()

    # Create empty out file based on input file
    outputLayer, outputDataSource = createLayer(outputFilePath, driverName, geomType)
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Get layer definition for the output file
    outputLayerDefinition = outputLayer.GetLayerDefn()

    # Create Coordinate Transform using EPSG codes
    inputSpatialRef = inputLayer.GetSpatialRef()
    outputSpatialRef = osr.SpatialReference()
    outputSpatialRef.ImportFromEPSG(outputEPSG)
    coordinateTransformation = osr.CoordinateTransformation(inputSpatialRef, outputSpatialRef)

    # Loop through the features in input file
    for i in range(inputLayer.GetFeatureCount()):

        # Get input feature
        inputFeature = inputLayer.GetFeature(i)

        # Get feature geometry and reproject
        geometry = inputFeature.GetGeometryRef()
        geometry.Transform(coordinateTransformation)

        # Create the output feature
        outputFeature = ogr.Feature(outputLayerDefinition)

        # Set geometry and field values
        outputFeature.SetGeometry(geometry)
        copyFeatureFieldValues(inputFeature, outputFeature)

        # Write out to output shapefile
        outputLayer.CreateFeature(outputFeature)

    # Generate .prj file
    createPrjFile(outputSpatialRef, outputFilePath)
    print 'Reprojected File: ' + inputFilePath

    return


def mergeSpatialFiles(inputFilePaths, outputFilePath, driverName):
    """
    This function stitches together a list of spatial files into one large file
    :param inputFilePaths: List of files to merge. (Files must all be in same projection)
    :param outputFilePath: New file containing merged features
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :return:
    """

    # Determine if we're working with a shapefile
    if driverName == 'ESRI Shapefile':
        isShapefile = True
    else:
        isShapefile = False

    # Get geometry from arbitrary input file
    inputLayer, inputDataSource = getLayer(inputFilePaths[0], driverName)
    geomType = inputLayer.GetLayerDefn().GetGeomType()

    # Create empty out file based on arbitrary input file
    outputLayer, outputDataSource = createLayer(outputFilePath, driverName, geomType)
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Set spatial reference for output file using arbitrary input file
    inputSpatialReference = inputLayer.GetSpatialRef()
    createPrjFile(inputSpatialReference, outputFilePath)

    # Loop through input files
    outputLayerDefinition = outputLayer.GetLayerDefn()
    for inputFilePath in inputFilePaths:

        # Read input file
        inputLayer, inputDataSource = getLayer(inputFilePath, driverName)

        # Loop through features in input file
        for i in range(inputLayer.GetFeatureCount()):

            # Get the input feature
            inputFeature = inputLayer.GetFeature(i)

            # Create the output feature
            outputFeature = ogr.Feature(outputLayerDefinition)

            # Set geometry and field values
            geometry = inputFeature.GetGeometryRef()
            outputFeature.SetGeometry(geometry)
            copyFeatureFieldValues(inputFeature, outputFeature, isShapefile)

            # Write out to output shapefile
            outputLayer.CreateFeature(outputFeature)

        print 'Processed File: ' + inputFilePath

    return


def addNewFieldToSpatialFile(filePath, driverName, fieldName, fieldType, fieldValues):
    """
    Adds a new field with values to a spatial file
    WARNING: Should not use on file already open in write mode!
    :param filePath: String. Path to existing spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param fieldName: Name to give new field
    :param fieldType: 'Integer', 'Float', or 'String'
    :param fieldValues: Must be a dictionary of form FID:Value
    :return
    """
    # Open file
    layer, dataSource = getLayer(filePath, driverName, mode=1)

    # Make sure fieldName is 10 characters or less because of shapefile limitations
    if len(fieldName) > 10:
        raise Exception('Field name must be 10 characters or less.')

    # Translate fieldType to OGR
    if fieldType == 'Integer':
        ogrFieldType = ogr.OFTInteger
    elif fieldType == 'Float':
        ogrFieldType = ogr.OFTReal
    elif fieldType == 'String':
        ogrFieldType = ogr.OFTString
    else:
        raise Exception('Field type not recognized: ' + fieldType)

    # Check for the field and create if necessary, overwrite otherwise
    fieldIndex = layer.GetLayerDefn().GetFieldIndex(fieldName)
    if fieldIndex != -1:
        layer.DeleteField(fieldIndex)
    fieldDefinition = ogr.FieldDefn(fieldName, ogrFieldType)
    layer.CreateField(fieldDefinition)

    # Loop through features to add values
    for FID, fieldValue in fieldValues.items():

        if fieldValue is None:
            print 'FID ' + str(FID) + ' had None value.'
            continue

        # Get the feature
        feature = layer.GetFeature(FID)

        # Set value of feature
        feature.SetField(fieldName, fieldValue)
        layer.SetFeature(feature)

    return


def idToFID(shapefilePath, idFieldName):
    """
    Maps explicit IDs named in one of the shapefile's fields to that feature's FID
    :param shapefilePath: String. Path to existing shapefile
    :param idFieldName: each value in this field should be a unique identifier
    :return: Dictionary in form {idFieldValue: FID}
    """
    # Open shapefile
    layer, dataSource = getLayer(shapefilePath, 'ESRI Shapefile')

    # Loop through all the features to build dictionary of id & FID
    idToFIDDictionary = {}
    feature = layer.GetNextFeature()
    while feature:
        fid = feature.GetFID()
        explicitID = feature.GetField(idFieldName)
        idToFIDDictionary[explicitID] = fid
        feature = layer.GetNextFeature()

    return idToFIDDictionary


def fieldDifference(filePath, driverName, fieldName1, fieldName2, newFieldName):
    """
    Finds the difference between two fields and adds that to the spatial file
    :param filePath: String. Path to existing spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param fieldName1: String. Name of first field
    :param fieldName2: String. Name of second field
    :param newFieldName: String. Name of newly generated field
    :return:
    """
    # Open spatial file
    layer, dataSource = getLayer(filePath, driverName)

    # Loop through layer to calculate difference between fields {FID:newFieldValue}
    newFieldValues = {}
    for feature in layer:
        fid = feature.GetFID()
        fieldValue1 = feature.GetField(fieldName1)
        fieldValue2 = feature.GetField(fieldName2)
        if fieldValue1 is None or fieldValue2 is None:
            continue
        newFieldValue = fieldValue1 - fieldValue2
        newFieldValues[fid] = newFieldValue

    # Feed that new dictionary into addNewFieldToSpatialFile
    addNewFieldToSpatialFile(filePath, driverName, newFieldName, 'Float', newFieldValues)
    return


def scaleField(filePath, driverName, existingFieldName, newFieldName):
    """
    Scales field from 0 to 1 and adds that information to the spatial file
    :param filePath: String. Path to existing spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param existingFieldName: String. Name of field of interest
    :param newFieldName: String. Name of new field
    :return:
    """
    # Open spatial file
    layer, dataSource = getLayer(filePath, driverName)

    # Loop through features to get minimum and maximum value
    minValue = float('inf')
    maxValue = float('-inf')
    for feature in layer:
        fieldValue = feature.GetField(existingFieldName)
        if fieldValue is None:
            continue
        if fieldValue < minValue:
            minValue = fieldValue
        if fieldValue > maxValue:
            maxValue = fieldValue

    # Loop through layer to calculate scaled version of the field {FID:newFieldValue}
    layer.ResetReading()
    scaledFieldValues = {}
    for feature in layer:

        fid = feature.GetFID()
        existingValue = feature.GetField(existingFieldName)

        # Pass None (NULL) values through
        if existingValue is None:
            continue

        # Discard 'nan' or 0 values, as these are likely data problems
        # TODO: Want to keep 0 in some cases, so maybe replacing 0 with NULL for certain fields should be own function
        # if pandas.isnull(existingValue) or existingValue == 0:
        if pandas.isnull(existingValue):
            scaledFieldValues[fid] = None
        else:
            scaledValue = ((existingValue - minValue) * 1.0) / (maxValue - minValue)
            scaledFieldValues[fid] = scaledValue

    # Feed that new dictionary into addNewFieldToSpatialFile
    addNewFieldToSpatialFile(filePath, driverName, newFieldName, 'Float', scaledFieldValues)
    return


def deleteField(shapefilePath, fieldName):
    """
    Deletes a field from a spatialFile
    WARNING: Should not run this function on a shapefile already open in write mode!
    :param shapefilePath: String. Path to existing shapefile
    :param fieldName: String. Name of field of interest
    :return:
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath, 1)
    layer = shapefile.GetLayer()

    # Delete field
    fieldIndex = layer.GetLayerDefn().GetFieldIndex(fieldName)
    if fieldIndex == -1:
        print('Field does not exist:', fieldName)
    else:
        layer.DeleteField(fieldIndex)

    return


def conditionallyDeleteFeature(shapefilePath, fieldName, fieldValue):
    """
    Deletes a feature who have a specific value for a certain field
    WARNING: Should not run this function on a shapefile already open in write mode!
    :param shapefilePath: String. Path to existing shapefile
    :param fieldName: String. Name of field of interest
    :param fieldValue: Value in field that triggers deletion of that feature
    :return:
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath, 1)
    layer = shapefile.GetLayer()

    # Check each feature for fieldValue
    for fid in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(fid)

        # Delete feature if value matches undesired value
        if feature.GetField(fieldName) == fieldValue:
            print 'Deleting field ' + str(fid)
            layer.DeleteFeature(fid)

    return


def osmToShapefile(inputOsmPath, layerToUse, outputShapefilePath):
    """
    This function converts an OSM file to a shapefile
    :param inputOsmPath: String. Can be .osm or .pbf file
    :param layerToUse: 0 for points, 1 for lines, 2 for multilines, 3 for multipolygons, 4 for other relations
    :param outputShapefilePath: String. Path to output location and file
    :return:
    """
    # Get geometry from input file
    inputLayer, inputDataSource = getLayer(inputOsmPath, 'OSM', osmLayer=layerToUse)
    geomType = inputLayer.GetLayerDefn().GetGeomType()

    # Create empty out file based on input file
    outputLayer, outputDataSource = createLayer(outputShapefilePath, 'ESRI Shapefile', geomType)
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Get layer definition for the output file
    outputLayerDefinition = outputLayer.GetLayerDefn()

    # Set spatial reference for output shapefile using arbitrary input file
    inputSpatialReference = inputLayer.GetSpatialRef()
    outputSpatialReference = osr.SpatialReference()
    outputSpatialReference.ImportFromWkt(inputSpatialReference.ExportToWkt())

    # Generate .prj file
    createPrjFile(outputSpatialReference, outputShapefilePath)

    # Loop through the features in input file
    for inputFeature in inputLayer:

        # Create the output feature
        outputFeature = ogr.Feature(outputLayerDefinition)

        # Set geometry and field values
        geometry = inputFeature.GetGeometryRef()
        outputFeature.SetGeometry(geometry)
        copyFeatureFieldValues(inputFeature, outputFeature)

        # Write out to output shapefile
        outputLayer.CreateFeature(outputFeature)

    print 'Created Shapefile from ' + inputOsmPath
    return


def pointCSVToShapefile(csvPath, x, y, shapefilePath, EPSG=4326):
    """
    This function creates a shapefile from a .csv file containing points and their attributes
    :param csvPath: String. Path to .csv containing data
    :param x: String. Name of field representing X (longitude)
    :param y: String. Name of field representing Y (latitude)
    :param shapefilePath: String. Path of shapefile that should be created
    :param EPSG: Int. Code to represent projection
    :return:
    """

    # Read in point data from CSV and get field names
    pointDataFrame = pandas.read_csv(csvPath)
    fieldNames = pointDataFrame.columns.values.tolist()
    fieldNames.remove(x)
    fieldNames.remove(y)

    # Create empty output file
    layer, dataSource = createLayer(shapefilePath, 'ESRI Shapefile', ogr.wkbPoint)

    # Add fields to the output shapefile
    for fieldName in fieldNames:
        fieldDefinition = ogr.FieldDefn()
        fieldDefinition.SetName(fieldName)
        layer.CreateField(fieldDefinition)

    # Add features to the shapefile based on feature definition from above
    featureDefinition = layer.GetLayerDefn()
    feature = ogr.Feature(featureDefinition)
    point = ogr.Geometry(ogr.wkbPoint)
    for index, values in pointDataFrame.iterrows():

        # Add geometry
        point.AddPoint(values[x], values[y])
        feature.SetGeometry(point)

        # Add field values
        for fieldName in fieldNames:
            feature.SetField(fieldName, values[fieldName])

        # Add feature to shapefile
        layer.CreateFeature(feature)

    # Define projection
    createPrjFile(EPSG, shapefilePath)
    return


def createPrjFile(projectionInfo, shapefilePath):
    """
    This function creates a .prj file for an existing shapefile according to a specified projection
    :param projectionInfo: EPSG code (as integer) or osr.SpatialReference object representing shapefile's projection
    :param shapefilePath: String. Path to shapefile we're specifying the projection for
    :return: String representing .prj file's path
    """

    # Get spatial reference information depending on how the projection info was passed in
    if isinstance(projectionInfo, int):
        spatialReference = osr.SpatialReference()
        spatialReference.ImportFromEPSG(projectionInfo)
    elif isinstance(projectionInfo, osr.SpatialReference):
        spatialReference = projectionInfo
    else:
        sys.exit('Projection info not recognized.')

    # Convert to ESRI's format
    spatialReference.MorphToESRI()

    # Write out to .prj file
    prjPath = shapefilePath[:-3] + 'prj'
    with open(prjPath, 'w') as prjFile:
        prjFile.write(spatialReference.ExportToWkt())

    return prjPath


def calculateDistance(place1, place2, method='nearestToNearest'):
    """
    This function calculates shortest distance (Euclidean) from place1 to place2
    :param place1: Starting Feature/Geometry
    :param place2: Ending Feature/Geometry
    :param method: String. Specifies how to calculate distance. Options are 'nearestToNearest' (nearest edge to nearest
    edge), 'nearestToCentroid' (nearest edge to centroid), 'centroidToNearest' (centroid to nearest edge), and
    'centroidToCentroid'. Each option is in order 'place1ToPlace2'
    :return: Float
    """

    # Get geometry object for each place, since Feature doesn't have a geometry method
    if isinstance(place1, ogr.Geometry):
        geometry1 = place1
    elif isinstance(place1, ogr.Feature):
        geometry1 = place1.GetGeometryRef()
    else:
        sys.exit('Place 1\'s geometry not derived from recognizable type.')

    if isinstance(place2, ogr.Geometry):
        geometry2 = place2
    elif isinstance(place2, ogr.Feature):
        geometry2 = place2.GetGeometryRef()
    else:
        sys.exit('Place 2\'s geometry not derived from recognizable type.')

    # Calculate distance according to specified method and return that value
    if method == 'nearestToNearest':
        return geometry1.Distance(geometry2)
    elif method == 'centroidToCentroid':
        return geometry1.Centroid().Distance(geometry2.Centroid())
    elif method == 'nearestToCentroid':
        return geometry1.Distance(geometry2.Centroid())
    elif method == 'centroidToNearest':
        return geometry1.Centroid().Distance(geometry2)


def getLayer(filePath, driverName, mode=0, osmLayer=None):
    """
    This function reads in an existing spatial file and returns its layer and data source (to keep it in scope)
    :param filePath: String. Path to spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param mode: 0 to read, 1 to write
    :param osmLayer: 0 for points, 1 for lines, 2 for multilines, 3 for multipolygons, 4 for other relations
    :return: Layer, DataSource
    """

    # Open the data source (with specific driver if given)
    driver = ogr.GetDriverByName(driverName)
    dataSource = driver.Open(filePath, mode)

    # Make sure data source was properly opened before getting the layer
    if dataSource is None:
        raise Exception('Could not open', dataSource)

    # Get the layer
    if osmLayer is not None:
        layer = dataSource.GetLayer(osmLayer)
    else:
        layer = dataSource.GetLayer()

    # Must return both layer and dataSource to previous scope to prevent layer from becoming unusable & causing crash
    return layer, dataSource


def createLayer(filePath, driverName, geometryType):
    """
    This function creates a new spatial file and returns its layer and data source (to keep it in scope)
    :param filePath: String. Location of new spatial file
    :param driverName: String. Type of spatial file (eg, 'ESRI Shapefile', 'OSM', 'GeoJSON', 'KML', 'SQLite')
    :param geometryType: Geometry. Ex: ogr.wkbPoint, ogr.wkbPolygon
    :return: Layer, DataSource
    """

    # Create appropriate driver
    driver = ogr.GetDriverByName(driverName)

    # Create empty file
    if os.path.exists(filePath):
        driver.DeleteDataSource(filePath)
    dataSource = driver.CreateDataSource(filePath)
    layerName = os.path.splitext(os.path.basename(filePath))[0]
    layer = dataSource.CreateLayer(layerName, geom_type=geometryType)

    # Must return both layer and dataSource to previous scope to prevent layer from becoming unusable & causing crash
    return layer, dataSource


def getFeatureList(layer):
    """
    Takes a layer and returns a list of all its feature objects
    :param layer: Layer object containing the set of features we want in list form
    :return: List of Features
    """
    features = []
    for feature in layer:
        features.append(feature)
    return features


def findNearest(searchFeature, targetLayer, distanceMethod='nearestToNearest',
                startBuffer=10000, bufferStep=200000, stopBuffer=1000000):
    """
    This finds the nearest feature to the "search Feature" in the targetLayer
    :param searchFeature: The feature or geometry we are starting our distance calculations from
    :param targetLayer: Layer object containing the set of features we are comparing the searchFeature to
    :param distanceMethod: Chosen method to calculate distance (see calculateDistance for options)
    :param startBuffer: Starting size of the search buffer in layer's units
    :param bufferStep: Increment of the search buffer in layer's units
    :param stopBuffer: Distance at which you are no longer interested in finding the nearest (the search ceiling)
    :return: nearest Feature
    """

    # Get geometry of search feature so we can buffer it
    searchGeometry = searchFeature.GetGeometryRef()

    # Set initial parameters for search
    currentBuffer = startBuffer
    nearest = None
    minDistance = float('inf')

    # Search the area around searchFeature for a targetFeature (avoids calculating distance for every targetFeature)
    while nearest is None and currentBuffer < stopBuffer:

        # Reduce the target layer to only those features within the buffered search area
        searchArea = searchGeometry.Buffer(currentBuffer)
        targetLayer.SetSpatialFilter(searchArea)

        # If any target features are within the search area, find which one is closest to the search feature
        if targetLayer.GetFeatureCount() > 0:

            for targetFeature in targetLayer:

                # If distance is less than the current minDistance, this is our new closest feature
                distance = calculateDistance(searchFeature, targetFeature, distanceMethod)
                if distance < minDistance:
                    nearest = targetFeature
                    minDistance = distance

        # Increment buffer to try the next largest search area
        currentBuffer += bufferStep

    # Clear the spatial filter
    targetLayer.SetSpatialFilter(None)

    return nearest


if __name__ == '__main__':

    # Set parameters
    skiShapefile = '../../Data/SkiAreas/skiAreas.shp'
    driverType = 'ESRI Shapefile'
    albersSkiShapefile = os.path.splitext(skiShapefile)[0] + '_AlbersEqualAreaConic.shp'

    # Read in layers
    skiAreaLayer, skiAreaDataSource = getLayer(albersSkiShapefile, driverType)
    subLayer, subDataSource = getLayer('../../Data/USMergedCountySubs/USMergedCountySubs_AlbersEqualAreaConic.shp',
                                       driverType)

    # Find the nearest ski area to feature with FID 100
    testFeature = subLayer.GetFeature(100)
    nearestSkiArea = findNearest(testFeature, skiAreaLayer)
    print nearestSkiArea.GetField('name')
