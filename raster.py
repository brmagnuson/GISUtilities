import sys
import numpy
import gdal
import ogr
import osr
import vector


def zonal_stats(feat, layer, raster, minAllowableValue, maxAllowableValue):

    # Get raster georeference info
    transform = raster.GetGeoTransform()
    xOrigin = transform[0]
    yOrigin = transform[3]
    pixelWidth = transform[1]
    pixelHeight = transform[5]

    # Get extent of feat
    geom = feat.GetGeometryRef()
    if (geom.GetGeometryName() == 'MULTIPOLYGON'):
        count = 0
        pointsX = []; pointsY = []
        for polygon in geom:
            geomInner = geom.GetGeometryRef(count)
            ring = geomInner.GetGeometryRef(0)
            numpoints = ring.GetPointCount()
            for p in range(numpoints):
                    lon, lat, z = ring.GetPoint(p)
                    pointsX.append(lon)
                    pointsY.append(lat)
            count += 1
    elif (geom.GetGeometryName() == 'POLYGON'):
        ring = geom.GetGeometryRef(0)
        numpoints = ring.GetPointCount()
        pointsX = []
        pointsY = []
        for p in range(numpoints):
                lon, lat, z = ring.GetPoint(p)
                pointsX.append(lon)
                pointsY.append(lat)

    else:
        sys.exit()

    xmin = min(pointsX)
    xmax = max(pointsX)
    ymin = min(pointsY)
    ymax = max(pointsY)

    # Specify offset and rows and columns to read
    xoff = int((xmin - xOrigin)/pixelWidth)
    yoff = int((yOrigin - ymax)/pixelWidth)
    xcount = int((xmax - xmin)/pixelWidth)+1
    ycount = int((ymax - ymin)/pixelWidth)+1

    # Create memory target raster
    target_ds = gdal.GetDriverByName('MEM').Create('', xcount, ycount, gdal.GDT_Byte)
    target_ds.SetGeoTransform((
        xmin, pixelWidth, 0,
        ymax, 0, pixelHeight,
    ))

    # Create for target raster the same projection as for the value raster
    raster_srs = osr.SpatialReference()
    raster_srs.ImportFromWkt(raster.GetProjectionRef())
    target_ds.SetProjection(raster_srs.ExportToWkt())

    # Rasterize zone polygon to raster
    gdal.RasterizeLayer(target_ds, [1], layer, burn_values=[1])

    # Read raster as arrays
    banddataraster = raster.GetRasterBand(1)
    dataraster = banddataraster.ReadAsArray(xoff, yoff, xcount, ycount).astype(numpy.float)

    bandmask = target_ds.GetRasterBand(1)
    datamask = bandmask.ReadAsArray(0, 0, xcount, ycount).astype(numpy.float)

    # Mask zone of raster, excluding places not in the polygon and places that the raster doesn't have data for
    zoneraster = numpy.ma.masked_array(dataraster,  numpy.logical_not(datamask) | [i == -9999 for i in dataraster])

    # Calculate statistics of zonal raster
    zoneMean = numpy.mean(zoneraster)

    # Set zone mean to mean of raster when zoneMean gets returned as '--' because zone is smaller than raster
    # resolution, in order to approximate a decent value for that zone
    if isinstance(zoneMean, numpy.ma.core.MaskedConstant):
        zoneMean = numpy.mean(dataraster)

    # Set zone mean to none when the above steps return an illogical value
    if zoneMean < minAllowableValue or zoneMean > maxAllowableValue:
        zoneMean = None

    return zoneMean


def loop_zonal_stats(shapefilePath, rasterPath, minAllowableValue, maxAllowableValue):
    """
    :param input_zone_polygon:
    :param input_value_raster:
    :return: Returns dictionary of Field ID: Mean raster value
    """

    # Open input files
    shapefile = ogr.Open(shapefilePath)
    raster = gdal.Open(rasterPath)

    # Set up all features layer
    allFeaturesLayer = shapefile.GetLayer()
    allFeaturesList = range(allFeaturesLayer.GetFeatureCount())
    allFeaturesLayerDefinition = allFeaturesLayer.GetLayerDefn()

    # Create in-memory data source for single-feature layer
    driver = ogr.GetDriverByName('Memory')
    dataSource = driver.CreateDataSource('inMemoryDataSource')
    singleFeatureLayer = dataSource.CreateLayer('singleFeatureLayer',
                                                srs = allFeaturesLayer.GetSpatialRef(),
                                                geom_type = allFeaturesLayerDefinition.GetGeomType())
    singleFeatureLayerDefinition = singleFeatureLayer.GetLayerDefn()
    vector.copyLayerFieldDefinitions(allFeaturesLayer, singleFeatureLayer)

    # Loop through each feature
    statDict = {}

    print 'Number of FIDs to process: ' + str(len(allFeaturesList))

    for FID in allFeaturesList:

        feature = allFeaturesLayer.GetFeature(FID)

        # Copy feature to layer
        singleFeature = ogr.Feature(singleFeatureLayerDefinition)
        singleFeature.SetGeometry(feature.GetGeometryRef())
        vector.copyFeatureFieldValues(feature, singleFeature)
        singleFeatureLayer.CreateFeature(singleFeature)

        # Get mean of raster within feature
        meanValue = zonal_stats(singleFeature, singleFeatureLayer, raster, minAllowableValue, maxAllowableValue)
        statDict[FID] = meanValue

        singleFeatureLayer.DeleteFeature(FID)
        singleFeature.Destroy()

        if FID % 25 == 0:
            print 'Processed FID: ' + str(FID)

    dataSource.Destroy()
    return statDict



