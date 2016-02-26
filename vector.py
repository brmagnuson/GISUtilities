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


def copyFeatureFieldValues(fromFeature, toFeature, toDriver):
    """Copy field values from one feature to another.
    This assumes that the features have the same fields!
    :param fromFeature: feature object that contains the data to copy
    :param toFeature: feature object that the data is to be copied into
    """
    for i in range(fromFeature.GetFieldCount()):
        fieldName = fromFeature.GetFieldDefnRef(i).GetName()

        # Set field as long as its value is not None
        if fromFeature.GetField(fieldName) != None:

            # Deal with shapefile field name length limitation of 10 characters
            if toDriver.GetName == 'ESRI Shapefile':
                if len(fieldName) > 10:
                    toFeature.SetField(fieldName[:10], fromFeature.GetField(fieldName))
                else:
                    toFeature.SetField(fieldName, fromFeature.GetField(fieldName))
            else:
                toFeature.SetField(fieldName, fromFeature.GetField(fieldName))


def reprojectShapefile(inputShapefilePath, outputShapefilePath, outputEPSG):
    """
    Reprojects a shapefile from one projected coordinate system to another based on EPSG number.
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # Create empty out file
    if os.path.exists(outputShapefilePath):
        driver.DeleteDataSource(outputShapefilePath)
    outputDataSource = driver.CreateDataSource(outputShapefilePath)
    outputLayerName = os.path.splitext(os.path.basename(outputShapefilePath))[0]

    # Set geometry and field definitions from input file
    inputDataSource = driver.Open(inputShapefilePath)
    if inputDataSource is None:
        raise Exception('Could not open', inputDataSource)
    inputLayer = inputDataSource.GetLayer()
    inputLayerDefinition = inputLayer.GetLayerDefn()
    outputLayer = outputDataSource.CreateLayer(outputLayerName, geom_type=inputLayerDefinition.GetGeomType())
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Get layer definition for the output file
    outputLayerDefinition = outputLayer.GetLayerDefn()

    # Create Coordinate Transform using EPSG codes
    inputSpatialRef = inputLayer.GetSpatialRef()
    outputSpatialRef = osr.SpatialReference()
    outputSpatialRef.ImportFromEPSG(outputEPSG)
    coordinateTransformation = osr.CoordinateTransformation(inputSpatialRef, outputSpatialRef)

    # Loop through the features in input file
    for i in range(0, inputLayer.GetFeatureCount()):

        # Get input feature
        inputFeature = inputLayer.GetFeature(i)

        # Get feature geometry and reproject
        geometry = inputFeature.GetGeometryRef()
        geometry.Transform(coordinateTransformation)

        # Create the output feature
        outputFeature = ogr.Feature(outputLayerDefinition)

        # Set geometry and field values
        outputFeature.SetGeometry(geometry)
        copyFeatureFieldValues(inputFeature, outputFeature, driver)

        # Write out to output shapefile
        outputLayer.CreateFeature(outputFeature)
        outputFeature.Destroy()

    # Generate .prj file
    outputSpatialRef.MorphToESRI()
    with open(outputShapefilePath[:-3] + 'prj', 'w') as file:
        file.write(outputSpatialRef.ExportToWkt())

    print 'Reprojected File: ' + inputShapefilePath

    # Close files
    inputDataSource.Destroy()
    outputDataSource.Destroy()


def mergeShapefiles(inputShapefilePaths, outputShapefilePath):
    """
    :param inputShapefilesPaths: List of shapefiles to merge. (Files must all be in same projection)
    :param outputShapefilePath: New shapefile containing merged features
    :return:
    """
    driver = ogr.GetDriverByName('ESRI Shapefile')

    # Create empty output file
    if os.path.exists(outputShapefilePath):
        driver.DeleteDataSource(outputShapefilePath)
    outputDataSource = driver.CreateDataSource(outputShapefilePath)
    outputLayerName = os.path.splitext(os.path.basename(outputShapefilePath))[0]

    # Set geometry and field definitions from arbitrary input file
    inputDataSource = driver.Open(inputShapefilePaths[0])
    inputLayer = inputDataSource.GetLayer()
    inputLayerDefinition = inputLayer.GetLayerDefn()
    outputLayer = outputDataSource.CreateLayer(outputLayerName, geom_type=inputLayerDefinition.GetGeomType())
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Set spatial reference for output shapefile using arbitrary input file
    inputSpatialReference = inputLayer.GetSpatialRef()
    outputSpatialReference = osr.SpatialReference()
    outputSpatialReference.ImportFromWkt(inputSpatialReference.ExportToWkt())
    outputSpatialReference.MorphToESRI()
    with open(outputShapefilePath[:-3] + 'prj', 'w') as file:
        file.write(outputSpatialReference.ExportToWkt())

    inputDataSource.Destroy()

    # Loop through input files
    outputLayerDefinition = outputLayer.GetLayerDefn()
    for inputShapefilePath in inputShapefilePaths:

        # Read input file
        inputDataSource = driver.Open(inputShapefilePath)
        inputLayer = inputDataSource.GetLayer()

        # Loop through features in input file
        for i in range(0, inputLayer.GetFeatureCount()):

            # Get the input feature
            inputFeature = inputLayer.GetFeature(i)

            # Create the output feature
            outputFeature = ogr.Feature(outputLayerDefinition)

            # Set geometry and field values
            geometry = inputFeature.GetGeometryRef()
            outputFeature.SetGeometry(geometry)
            copyFeatureFieldValues(inputFeature, outputFeature, driver)

            # Write out to output shapefile
            outputLayer.CreateFeature(outputFeature)

            outputFeature.Destroy()

        # Close input file
        inputDataSource.Destroy()

        print 'Processed File: ' + inputShapefilePath

    # Close output file
    outputDataSource.Destroy()


def addNewFieldToShapefile(shapefilePath, fieldName, fieldType, fieldValues):
    """
    :param shapefilePath:
    :param fieldName: Name to give new field
    :param fieldType: 'Integer', 'Float', or 'String'
    :param fieldValues: Must be a dictionary of form FID:Value
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath, 1)
    layer = shapefile.GetLayer()

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

        if FID % 25 == 0:
            print 'Processed FID: ' + str(FID)


def idToFID(shapefilePath, idFieldName):
    """
    Maps explicit IDs named in one of the shapefile's fields to that feature's FID
    :param shapefilePath:
    :param idFieldName: each value in this field should be a unique identifier
    :return: dictionary in form {idFieldValue: FID}
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath)
    layer = shapefile.GetLayer()

    # Loop through all the features to build dictionary of id & FID
    idToFIDDictionary = {}
    feature = layer.GetNextFeature()
    while feature:
        fid = feature.GetFID()
        id = feature.GetField(idFieldName)
        idToFIDDictionary[id] = fid
        feature = layer.GetNextFeature()

    return idToFIDDictionary


def fieldDifference(shapefilePath, fieldName1, fieldName2, newFieldName):
    """

    :param shapefilePath:
    :param existingFieldName:
    :param newFieldName:
    :return:
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath)
    layer = shapefile.GetLayer()

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
    addNewFieldToShapefile(shapefilePath, newFieldName, 'Float', newFieldValues)


def scaleField(shapefilePath, existingFieldName, newFieldName):
    """

    :param shapefilePath:
    :param existingFieldName:
    :param newFieldName:
    :return:
    """
    # Open shapefile
    driver = ogr.GetDriverByName('ESRI Shapefile')
    shapefile = driver.Open(shapefilePath)
    layer = shapefile.GetLayer()

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
    addNewFieldToShapefile(shapefilePath, newFieldName, 'Float', scaledFieldValues)


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


def osmToShapefile(inputOsmPath, layerToUse, outputShapefilePath):
    """

    :param inputOsmPath: can be .osm or .pbf file
    :param layerToUse: 0 for points, 1 for lines, 2 for multilines, 3 for multipolygons, 4 for other relations
    :param outputShapefilePath:
    :return:
    """
    # Open input file
    osmDriver = ogr.GetDriverByName('OSM')
    inputDataSource = osmDriver.Open(inputOsmPath)
    if inputDataSource is None:
        raise Exception('Could not open', inputDataSource)
    inputLayer = inputDataSource.GetLayer(layerToUse)

    # Create empty output file
    shapefileDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(outputShapefilePath):
        shapefileDriver.DeleteDataSource(outputShapefilePath)
    outputDataSource = shapefileDriver.CreateDataSource(outputShapefilePath)
    outputLayerName = os.path.splitext(os.path.basename(outputShapefilePath))[0]

    # Set geometry and field definitions from input file
    inputLayerDefinition = inputLayer.GetLayerDefn()
    outputLayer = outputDataSource.CreateLayer(outputLayerName, geom_type=inputLayerDefinition.GetGeomType())
    copyLayerFieldDefinitions(inputLayer, outputLayer)

    # Get layer definition for the output file
    outputLayerDefinition = outputLayer.GetLayerDefn()

    # Set spatial reference for output shapefile using arbitrary input file
    inputSpatialReference = inputLayer.GetSpatialRef()
    outputSpatialReference = osr.SpatialReference()
    outputSpatialReference.ImportFromWkt(inputSpatialReference.ExportToWkt())

    # Generate .prj file
    outputSpatialReference.MorphToESRI()
    with open(outputShapefilePath[:-3] + 'prj', 'w') as file:
        file.write(outputSpatialReference.ExportToWkt())

    # Loop through the features in input file
    for inputFeature in inputLayer:

        # Create the output feature
        outputFeature = ogr.Feature(outputLayerDefinition)

        # Set geometry and field values
        geometry = inputFeature.GetGeometryRef()
        outputFeature.SetGeometry(geometry)
        copyFeatureFieldValues(inputFeature, outputFeature, shapefileDriver)

        # Write out to output shapefile
        outputLayer.CreateFeature(outputFeature)
        outputFeature.Destroy()

    print 'Created Shapefile from ' + inputOsmPath

    # Close files
    inputDataSource.Destroy()
    outputDataSource.Destroy()


def pointCSVToShapefile(csvPath, x, y, shapefilePath, EPSG=4326):

    # Read in point data from CSV and get field names
    pointDataFrame = pandas.read_csv(csvPath)
    fieldNames = pointDataFrame.columns.values.tolist()
    fieldNames.remove(x)
    fieldNames.remove(y)

    # Create empty output file
    shapefileDriver = ogr.GetDriverByName('ESRI Shapefile')
    if os.path.exists(shapefilePath):
        shapefileDriver.DeleteDataSource(shapefilePath)
    dataSource = shapefileDriver.CreateDataSource(shapefilePath)
    layerName = os.path.splitext(os.path.basename(shapefilePath))[0]
    layer = dataSource.CreateLayer(layerName)

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

    dataSource.Destroy()
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
    with open(prjPath, 'w') as file:
        file.write(spatialReference.ExportToWkt())

    return prjPath


