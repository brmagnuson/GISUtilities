import sys
import ogr
import os
import osr
import pandas


def copyLayerFieldDefinitions(fromLayer, toLayer):
    """Copy field definitions (not data) from one layer to another.
    :param fromLayer: layer object that contains the field definitions to copy
    :param toLayer: layer object to copy the field definitions into
    """
    featureDefinition = fromLayer.GetLayerDefn()
    for i in range(featureDefinition.GetFieldCount()):
        fieldDefinition = featureDefinition.GetFieldDefn(i)
        toLayer.CreateField(fieldDefinition)

    return


def copyFeatureFieldValues(fromFeature, toFeature, isShapefile=True):
    """Copy field values from one feature to another.
    This assumes that the features have the same fields!
    :param fromFeature: feature object that contains the data to copy
    :param toFeature: feature object that the data is to be copied into
    :param isShapefile: True if toFeature is from a shapefile, False otherwise.
    """
    for i in range(fromFeature.GetFieldCount()):
        fieldName = fromFeature.GetFieldDefnRef(i).GetName()

        # Set field as long as its value is not None
        if fromFeature.GetField(fieldName) != None:

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
    :param inputFilePaths: List of files to merge. (Files must all be in same projection)
    :param outputFilePath: New file containing merged features
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
    :param filePath:
    :param fieldName: Name to give new field
    :param fieldType: 'Integer', 'Float', or 'String'
    :param fieldValues: Must be a dictionary of form FID:Value
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

        if fieldValue == None:
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
    :param shapefilePath:
    :param idFieldName: each value in this field should be a unique identifier
    :return: dictionary in form {idFieldValue: FID}
    """
    # Open shapefile
    layer, dataSource = getLayer(shapefilePath, 'ESRI Shapefile')

    # Loop through all the features to build dictionary of id & FID
    idToFIDDictionary = {}
    feature = layer.GetNextFeature()
    while feature:
        fid = feature.GetFID()
        id = feature.GetField(idFieldName)
        idToFIDDictionary[id] = fid
        feature = layer.GetNextFeature()

    return idToFIDDictionary


def fieldDifference(filePath, driverName, fieldName1, fieldName2, newFieldName):
    """
    Finds the difference between two fields and adds that two the spatial file
    :param filePath:
    :param driverName:
    :param existingFieldName:
    :param newFieldName:
    :return:
    """
    # Open shapefile
    layer, dataSource = getLayer(filePath, driverName)

    # Loop through layer to calculate difference between fields {FID:newFieldValue}
    newFieldValues = {}
    for feature in layer:
        fid = feature.GetFID()
        fieldValue1 = feature.GetField(fieldName1)
        fieldValue2 = feature.GetField(fieldName2)
        if fieldValue1 == None or fieldValue2 == None:
            continue
        newFieldValue = fieldValue1 - fieldValue2
        newFieldValues[fid] = newFieldValue

    # Feed that new dictionary into addNewFieldToShapefile
    addNewFieldToSpatialFile(filePath, driverName, newFieldName, 'Float', newFieldValues)
    return


def scaleField(filePath, driverName, existingFieldName, newFieldName):
    """
    Scales field from 0 to 1 and adds that to the shapefile
    :param filePath:
    :param existingFieldName:
    :param newFieldName:
    :return:
    """
    # Open shapefile
    layer, dataSource = getLayer(filePath, driverName)

    # Loop through features to get minimum and maximum value
    minValue = float('inf')
    maxValue = float('-inf')
    for feature in layer:
        fieldValue = feature.GetField(existingFieldName)
        if fieldValue == None:
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
        if existingValue == None:
            continue

        # Discard 'nan' or 0 values, as these are likely data problems
        if pandas.isnull(existingValue) or existingValue == 0:
            scaledFieldValues[fid] = None
        else:
            scaledValue = ((existingValue - minValue) * 1.0) / (maxValue - minValue)
            scaledFieldValues[fid] = scaledValue

    # Feed that new dictionary into addNewFieldToShapefile
    addNewFieldToSpatialFile(filePath, driverName, newFieldName, 'Float', scaledFieldValues)
    return


def deleteField(shapefilePath, fieldName):
    """
    Should not run this function on a shapefile already open in write mode.
    :param shapefilePath:
    :param fieldName:
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
    Should not run this function on a shapefile already open in write mode.
    :param shapefilePath:
    :param fieldName:
    :param fieldValue:
    :return:
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath, 1)
    layer = shapefile.GetLayer()

    # Check each feature for fieldValue
    for fid in range(layer.GetFeatureCount()):
        feature = layer.GetFeature(fid)
        if feature.GetField(fieldName) == fieldValue:
            print 'Deleting field ' + fid
            layer.DeleteFeature(fid)

    return


def osmToShapefile(inputOsmPath, layerToUse, outputShapefilePath):
    """

    :param inputOsmPath: can be .osm or .pbf file
    :param layerToUse: 0 for points, 1 for lines, 2 for multilines, 3 for multipolygons, 4 for other relations
    :param outputShapefilePath:
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

    if isinstance(projectionInfo, int):
        spatialReference = osr.SpatialReference()
        spatialReference.ImportFromEPSG(projectionInfo)
    elif isinstance(projectionInfo, osr.SpatialReference):
        spatialReference = projectionInfo
    else:
        sys.exit('Projection info not recognized.')

    spatialReference.MorphToESRI()

    prjPath = shapefilePath[:-3] + 'prj'
    with open(prjPath, 'w') as prjFile:
        prjFile.write(spatialReference.ExportToWkt())

    return prjPath


def calculateDistance(place1, place2, method='nearestToNearest'):

    # Get geometry
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

    # Calculate distance according to specified method
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

    :param filePath:
    :param driverName:
    :param mode: 0 to read, 1 to write
    :param osmLayer: 0 for points, 1 for lines, 2 for multilines, 3 for multipolygons, 4 for other relations
    :return:
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

    # Must return both layer and dataSource to previous scope to prevent layer from becoming unusable & causing crash.
    return layer, dataSource


def createLayer(filePath, driverName, geometryType):

    driver = ogr.GetDriverByName(driverName)

    # Create empty file
    if os.path.exists(filePath):
        driver.DeleteDataSource(filePath)
    dataSource = driver.CreateDataSource(filePath)
    layerName = os.path.splitext(os.path.basename(filePath))[0]
    layer = dataSource.CreateLayer(layerName, geom_type=geometryType)

    # Must return both layer and dataSource to previous scope to prevent layer from becoming unusable & causing crash.
    return layer, dataSource


def getFeatureList(layer):
    """
    Takes a layer and returns a list of all its feature objects
    :param layer:
    :return:
    """
    features = []
    for feature in layer:
        features.append(feature)
    return features


def findNearest(searchFeature, targetLayer, distanceMethod='nearestToNearest',
                startBuffer=10000, bufferStep=200000, stopBuffer=1000000):
    """
    This finds the nearest feature to the "search Feature" in the targetLayer.
    :param searchFeature:
    :param targetLayer:
    :return:
    """

    # Get geometry of search feature so we can buffer it
    searchGeometry = searchFeature.GetGeometryRef()

    currentBuffer = startBuffer
    nearest = None
    minDistance = float('inf')

    while nearest is None and currentBuffer < stopBuffer:

        # Reduce the target layer to only those features within the buffered search area
        searchArea = searchGeometry.Buffer(currentBuffer)
        targetLayer.SetSpatialFilter(searchArea)

        # If any target features are within the search area, find which one is closest to the search feature.
        if targetLayer.GetFeatureCount() > 0:
            for targetFeature in targetLayer:
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
    skiShapefile = '../../Data/SkiAreas/skiAreas.shp'
    driverName = 'ESRI Shapefile'
    albersSkiShapefile = os.path.splitext(skiShapefile)[0] + '_AlbersEqualAreaConic.shp'

    skiAreaLayer, skiAreaDataSource = getLayer(albersSkiShapefile, driverName)
    subLayer, subDataSource = getLayer('../../Data/USMergedCountySubs/USMergedCountySubs_AlbersEqualAreaConic.shp', driverName)

    testFeature = subLayer.GetFeature(100)
    nearestSkiArea = findNearest(testFeature, skiAreaLayer)
    print nearestSkiArea



