/*  file_Date = $Date: 2012-06-21 17:14:58 +0000 (Thu, 21 Jun 2012) $
 *  file_Revision = $Revision: 292 $
 *  file_Author = $Author: geoffc $
 *  file_HeadURL = $HeadURL: https://svn.ssec.wisc.edu/repos/geoffc/C/griddingAndGranulation/griddingAndGranulation.c $
 *  file_Id = $Id: griddingAndGranulation.c 292 2012-06-21 17:14:58Z geoffc $
 *
 *  ________________________
 *  griddingAndGranulation.c
 *  ________________________
 *
 *  A library for transforming datasets between regular and
 *  non-regular grids. 
 *
 *  Created by Geoff Cureton on 2011-06-27.
 *  Copyright (c) 2011 University of Wisconsin SSEC. All rights reserved.
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define TRUE 1
#define FALSE 0

int test_ctypes(double *dblArr, long longVar, int intVar)
{

    long idx;

    printf("longVar = %ld\n",longVar);
    printf("intVar = %d\n",intVar);

    printf("dblArr = \n");
    for (idx=0;idx<longVar;idx++){
        printf("%f\n",dblArr[idx]);
    }

    return(0);
}

/*
 *  _________
 *  gran2grid
 *  _________
 *
 *  This routine takes the data[nData] array, defined on the non-regular grid
 *  defined by lat[ndata] and lon[nData], and regrids to the regular grid 
 *  gridData[nGridRows*nGridCols] defined by gridLat[nGridRows*nGridCols] and 
 *  gridLon[nGridRows*nGridRows].
 *
 *  Inputs:
 *
 *  lat         : Non-regular latitude grid (nData)
 *  lon         : Non-regular longitude grid (nData)
 *  data        : Data array on non-regular grid defined by lat and lon (nData)
 *  nData       : Length of non-regular grid arrays
 *  gridLat     : Regular latitude grid (nGridRows*nGridCols)
 *  gridLon     : Regular longitude grid (nGridRows*nGridCols)
 *  nGridRows   : Number of grid rows of the regular grid
 *  nGridCols   : Number of grid columns of the regular grid
 *
 *  Outputs:
 *
 *  gridData    : Data array on grid defined by gridLat and gridLon 
 *                grid (nGridRows*nGridCols)
 *  gridDataIdx : Array of indices into data, defined by gridLat and gridLon 
 *                grid (nGridRows*nGridCols)
 *
 *
 */
int gran2grid(double *lat, 
              double *lon, 
              double *data, 
              long nData, 
              double *gridLat,
              double *gridLon,
              double *gridData,
              long *gridDataIdx,
              int nGridRows,
              int nGridCols
              )
{

    long idx;
    int latGridIdxLo, latGridIdxHi, lonGridIdxLo, lonGridIdxHi;

    double latVal,lonVal,dataVal;
    double gridLatInc,gridLonInc;

    int gridLatPt,gridLonPt;
    int gridLatPoints[4], gridLonPoints[4];
    double minDist,dist,latDist,lonDist;
    double gridLatVal,gridLonVal;
    int crnrPt,snapCrnrPt;
    int rowInBounds,colInBounds;

    printf("Shape of data is (%ld)\n",nData);
    printf("Shape of gridData is (%d, %d)\n", nGridRows,nGridCols);
    printf("TRUE = %d\n", TRUE);
    printf("FALSE = %d\n", FALSE);

    int numShow = 10;
    int i;
    printf("\nDisplaying lat,lon,data...\n");
    for (i=0;i<numShow;i++){
        printf("%6.1f %6.1f %6.1f\n",lat[i],lon[i],data[i]);
    }

    // Determine the lat and lon grid spacings
    printf("nGridRows = %d\n", nGridRows);
    printf("nGridCols = %d\n", nGridCols);
    printf("gridLat[nGridCols] = %f\n", gridLat[nGridCols]);
    printf("gridLat[0] = %f\n", gridLat[0]);
    printf("gridLon[1] = %f\n", gridLon[1]);
    printf("gridLon[0] = %f\n", gridLon[0]);

    gridLatInc = fabs(gridLat[nGridCols]-gridLat[0]);
    gridLonInc = fabs(gridLon[1]-gridLon[0]);

    printf("gridLatInc,gridLonInc = (%6.1f %6.1f)\n",gridLatInc,gridLonInc);

    // Loop through non-gridded data, find matching gridpoint and assign
    // data value to that gridpoint

    for (idx=0;idx<nData;idx++){
    /*for (idx=0;idx<numShow;idx++){*/
        latVal = lat[idx];
        lonVal = lon[idx];
        dataVal = data[idx];

        // Determine lat/lon grid indices which bound the non-gridded point
        latGridIdxLo = (int) floor((latVal-gridLat[0])/gridLatInc);
        latGridIdxHi = latGridIdxLo + 1;
        lonGridIdxLo = (int) floor((lonVal-gridLon[0])/gridLonInc);
        lonGridIdxHi = lonGridIdxLo + 1;

        rowInBounds = TRUE;
        colInBounds = TRUE;

        // If the grid indices bounding the non-gridded point are off the 
        // grid, mark this non-gridded point as out-of-bounds.
        if ((latGridIdxLo<0) || (latGridIdxHi>=nGridRows)){
            rowInBounds = FALSE;
            /*printf("Row idx out of bounds...\n");*/
        }
        if ((lonGridIdxLo<0) || (lonGridIdxHi>=nGridCols)){
            colInBounds = FALSE;
            /*printf("Column idx out of bounds...\n");*/
        }

        if (rowInBounds==FALSE){
            continue;
        }else if (colInBounds==FALSE){
            continue;
        }else{
            gridLatPoints[0] = latGridIdxLo;
            gridLatPoints[1] = latGridIdxLo;
            gridLatPoints[2] = latGridIdxHi;
            gridLatPoints[3] = latGridIdxHi;

            gridLonPoints[0] = lonGridIdxLo;
            gridLonPoints[1] = lonGridIdxHi;
            gridLonPoints[2] = lonGridIdxLo;
            gridLonPoints[3] = lonGridIdxHi;

            minDist = 1.e+30;
            snapCrnrPt = 0;

            // Loop through the corners bounding the non-gridded point
            for (crnrPt=0;crnrPt<4;crnrPt++){

                gridLatPt = (int) gridLatPoints[crnrPt];
                gridLonPt = (int) gridLonPoints[crnrPt];

                // Get the lat and lon values at this corner
                gridLatVal = gridLat[nGridCols * gridLatPt + gridLonPt];
                gridLonVal = gridLon[nGridCols * gridLatPt + gridLonPt];

                // The Pythagorean distance of this corner from the non-gridded
                // data point...
                latDist = latVal-gridLatVal;
                lonDist = lonVal-gridLonVal;
                dist = sqrt(latDist*latDist + lonDist*lonDist); 

                // If this distance is the smallest so far, save the corner
                // idx and distance...
                if (dist < minDist){
                    snapCrnrPt = crnrPt;
                    minDist = dist;
                }
            }

            // Get the grid indices of the corner closest to the non-gridded point
            gridLatPt = (int) gridLatPoints[snapCrnrPt];
            gridLonPt = (int) gridLonPoints[snapCrnrPt];

            // Assign the value of the non-gridded data point to the grid point 
            // closest to it.
            gridData[nGridCols * gridLatPt + gridLonPt] = dataVal;
            /*printf("gridDataIdx indices are (%d, %d)\n",gridLatPt,gridLonPt);*/

            // Save the index of the non-gridded data point to the nearest grid 
            // point.
            gridDataIdx[nGridCols * gridLatPt + gridLonPt] = idx;
            
            // TODO : Save subsequent data indices in the same gridcell
        }
    }

    return(0);
}

/*
 *  _________
 *  gran2grid_moments
 *  _________
 *
 *  This routine takes the data[nData] array, defined on the non-regular grid
 *  defined by lat[ndata] and lon[nData], and regrids to the regular grid 
 *  gridData[nGridRows*nGridCols] defined by gridLat[nGridRows*nGridCols] and 
 *  gridLon[nGridRows*nGridRows].
 *
 *  In order to facilitate computation of the various statistical descriptors, 
 *  we save on the defined grid the various sums (x, x^2, x^3, etc...), and the 
 *  number of observations added to each grid cell. From these datasets the various 
 *  gridded cumulants can be computed.
 *
 *  Inputs:
 *
 *  lat         : Non-regular latitude grid (nData)
 *  lon         : Non-regular longitude grid (nData)
 *  data        : Data array on non-regular grid defined by lat and lon (nData)
 *  nData       : Length of non-regular grid arrays
 *  gridLat     : Regular latitude grid (nGridRows*nGridCols)
 *  gridLon     : Regular longitude grid (nGridRows*nGridCols)
 *  nGridRows   : Number of grid rows of the regular grid
 *  nGridCols   : Number of grid columns of the regular grid
 *
 *  Outputs:
 *
 *  gridData    : Data array on grid defined by gridLat and gridLon 
 *                grid (nGridRows*nGridCols)
 *  gridDataIdx : Array of indices into data, defined by gridLat and gridLon 
 *                grid (nGridRows*nGridCols)
 *
 *
 */
int gran2grid_moments(double *lat, 
              double *lon, 
              double *data, 
              long nData, 
              double *gridLat,
              double *gridLon,
              double *gridData,
              long *gridDataIdx,
              int nGridRows,
              int nGridCols
              )
{

    long idx;
    int latGridIdxLo, latGridIdxHi, lonGridIdxLo, lonGridIdxHi;

    double latVal,lonVal,dataVal;
    double gridLatInc,gridLonInc;

    int gridLatPt,gridLonPt;
    int gridLatPoints[4], gridLonPoints[4];
    double minDist,dist,latDist,lonDist;
    double gridLatVal,gridLonVal;
    int crnrPt,snapCrnrPt;
    int rowInBounds,colInBounds;

    printf("Shape of data is (%ld)\n",nData);
    printf("Shape of gridData is (%d, %d)\n", nGridRows,nGridCols);
    printf("TRUE = %d\n", TRUE);
    printf("FALSE = %d\n", FALSE);

    int numShow = 10;
    int i;
    printf("\nDisplaying lat,lon,data...\n");
    for (i=0;i<numShow;i++){
        printf("%6.1f %6.1f %6.1f\n",lat[i],lon[i],data[i]);
    }

    // Determine the lat and lon grid spacings
    printf("nGridRows = %d\n", nGridRows);
    printf("nGridCols = %d\n", nGridCols);
    printf("gridLat[nGridCols] = %f\n", gridLat[nGridCols]);
    printf("gridLat[0] = %f\n", gridLat[0]);
    printf("gridLon[1] = %f\n", gridLon[1]);
    printf("gridLon[0] = %f\n", gridLon[0]);

    gridLatInc = fabs(gridLat[nGridCols]-gridLat[0]);
    gridLonInc = fabs(gridLon[1]-gridLon[0]);

    printf("gridLatInc,gridLonInc = (%6.1f %6.1f)\n",gridLatInc,gridLonInc);

    // Loop through non-gridded data, find matching gridpoint and assign
    // data value to that gridpoint

    for (idx=0;idx<nData;idx++){
        latVal = lat[idx];
        lonVal = lon[idx];
        dataVal = data[idx];

        // Determine lat/lon grid indices which bound the non-gridded point
        latGridIdxLo = (int) floor((latVal-gridLat[0])/gridLatInc);
        latGridIdxHi = latGridIdxLo + 1;
        lonGridIdxLo = (int) floor((lonVal-gridLon[0])/gridLonInc);
        lonGridIdxHi = lonGridIdxLo + 1;

        rowInBounds = TRUE;
        colInBounds = TRUE;

        // If the grid indices bounding the non-gridded point are off the 
        // grid, mark this non-gridded point as out-of-bounds.
        if ((latGridIdxLo<0) || (latGridIdxHi>=nGridRows)){
            rowInBounds = FALSE;
        }
        if ((lonGridIdxLo<0) || (lonGridIdxHi>=nGridCols)){
            colInBounds = FALSE;
        }

        if (rowInBounds==FALSE){
            continue;
        }else if (colInBounds==FALSE){
            continue;
        }else{
            gridLatPoints[0] = latGridIdxLo;
            gridLatPoints[1] = latGridIdxLo;
            gridLatPoints[2] = latGridIdxHi;
            gridLatPoints[3] = latGridIdxHi;

            gridLonPoints[0] = lonGridIdxLo;
            gridLonPoints[1] = lonGridIdxHi;
            gridLonPoints[2] = lonGridIdxLo;
            gridLonPoints[3] = lonGridIdxHi;

            minDist = 1.e+30;
            snapCrnrPt = 0;

            // Loop through the corners bounding the non-gridded point
            for (crnrPt=0;crnrPt<4;crnrPt++){

                gridLatPt = (int) gridLatPoints[crnrPt];
                gridLonPt = (int) gridLonPoints[crnrPt];

                // Get the lat and lon values at this corner
                gridLatVal = gridLat[nGridCols * gridLatPt + gridLonPt];
                gridLonVal = gridLon[nGridCols * gridLatPt + gridLonPt];

                // The Pythagorean distance of this corner from the non-gridded
                // data point...
                latDist = latVal-gridLatVal;
                lonDist = lonVal-gridLonVal;
                dist = sqrt(latDist*latDist + lonDist*lonDist); 

                // If this distance is the smallest so far, save the corner
                // idx and distance...
                if (dist < minDist){
                    snapCrnrPt = crnrPt;
                    minDist = dist;
                }
            }

            // Get the grid indices of the corner closest to the non-gridded point
            gridLatPt = (int) gridLatPoints[snapCrnrPt];
            gridLonPt = (int) gridLonPoints[snapCrnrPt];

            // Add the value of the non-gridded data point to the grid point 
            // closest to it.
            gridData[nGridCols * gridLatPt + gridLonPt] += dataVal;

            // Increment the number of obervations which have been attributed to this grid 
            // point.
            gridDataIdx[nGridCols * gridLatPt + gridLonPt] += 1;
            
        }
    }

    return(0);
}

/*
 *  _________
 *  grid2gran
 *  _________
 *
 *  This routine takes the regular grid gridData[nGridRows*nGridCols], 
 *  defined by gridLat[nGridRows*nGridCols] and gridLon[nGridRows*nGridRows],
 *  and regrids to the data[nData] array, defined on the non-regular grid
 *  defined by lat[ndata] and lon[nData].
 *
 *  Inputs:
 *
 *  gridLat     : Non-regular latitude grid (nGridRows*nGridCols)
 *  gridLon     : Non-regular longitude grid (nGridRows*nGridCols)
 *  gridData    : Data array on grid defined by gridLat and gridLon 
 *                grid (nGridRows*nGridCols)
 *  lat         : Non-regular latitude grid (nData)
 *  lon         : Non-regular longitude grid (nData)
 *  nData       : Length of non-regular grid arrays
 *  nGridRows   : Number of grid rows of the regular grid
 *  nGridCols   : Number of grid columns of the regular grid
 *
 *  Outputs:
 *
 *  data        : Data array on non-regular grid defined by lat and lon (nData)
 *  gridDataIdx : Array of indices into gridData, defined by lat and lon 
 *                non-regular grid (nData)
 *
 *
 */
int grid2gran(double *lat, 
              double *lon, 
              double *data, 
              long nData, 
              double *gridLat,
              double *gridLon,
              double *gridData,
              long *gridDataIdx,
              int nGridRows,
              int nGridCols
              )
{

    long idx;
    int latGridIdxLo, latGridIdxHi, lonGridIdxLo, lonGridIdxHi;

    double latVal,lonVal,dataVal;
    double gridLatInc,gridLonInc;

    int gridLatPt,gridLonPt;
    int gridLatPoints[4], gridLonPoints[4];
    double minDist,dist,latDist,lonDist;
    double gridLatVal,gridLonVal;
    int crnrPt,snapCrnrPt;
    int rowInBounds,colInBounds;

    /*printf("Shape of data is (%ld)\n",nData);*/
    /*printf("Shape of gridData is (%d, %d)\n", nGridRows,nGridCols);*/
    /*printf("TRUE = %d\n", TRUE);*/
    /*printf("FALSE = %d\n", FALSE);*/

    int numShow = 10;
    int i;
    /*printf("\nDisplaying lat,lon,data...\n");*/
    /*for (i=0;i<numShow;i++){*/
        /*printf("%6.1f %6.1f %6.1f\n",lat[i],lon[i],data[i]);*/
    /*}*/

    // Determine the lat and lon grid spacings
    /*printf("nGridRows = %d\n", nGridRows);*/
    /*printf("nGridCols = %d\n", nGridCols);*/
    /*printf("gridLat[nGridCols] = %f\n", gridLat[nGridCols]);*/
    /*printf("gridLat[0] = %f\n", gridLat[0]);*/
    /*printf("gridLon[1] = %f\n", gridLon[1]);*/
    /*printf("gridLon[0] = %f\n", gridLon[0]);*/

    gridLatInc = fabs(gridLat[nGridCols]-gridLat[0]);
    gridLonInc = fabs(gridLon[1]-gridLon[0]);

    /*printf("gridLatInc,gridLonInc = (%6.1f %6.1f)\n",gridLatInc,gridLonInc);*/

    // Loop through non-gridded data points, find matching gridpoint, and assign
    // gridpoint data value to that non-gridded data point.

    for (idx=0;idx<nData;idx++){
    /*for (idx=0;idx<numShow;idx++){*/
        latVal = lat[idx];
        lonVal = lon[idx];

        // Determine lat/lon grid indices which bound the non-gridded point
        latGridIdxLo = (int) floor((latVal-gridLat[0])/gridLatInc);
        latGridIdxHi = latGridIdxLo + 1;
        lonGridIdxLo = (int) floor((lonVal-gridLon[0])/gridLonInc);
        lonGridIdxHi = lonGridIdxLo + 1;

        rowInBounds = TRUE;
        colInBounds = TRUE;

        // If the grid indices bounding the non-gridded point are off the 
        // grid, mark this non-gridded point as out-of-bounds.
        if ((latGridIdxLo<0) || (latGridIdxHi>=nGridRows)){
            rowInBounds = FALSE;
            /*printf("Row idx out of bounds...\n");*/
        }
        if ((lonGridIdxLo<0) || (lonGridIdxHi>=nGridCols)){
            colInBounds = FALSE;
            /*printf("Column idx out of bounds...\n");*/
        }

        if (rowInBounds==FALSE){
            continue;
        }else if (colInBounds==FALSE){
            continue;
        }else{
            gridLatPoints[0] = latGridIdxLo;
            gridLatPoints[1] = latGridIdxLo;
            gridLatPoints[2] = latGridIdxHi;
            gridLatPoints[3] = latGridIdxHi;

            gridLonPoints[0] = lonGridIdxLo;
            gridLonPoints[1] = lonGridIdxHi;
            gridLonPoints[2] = lonGridIdxLo;
            gridLonPoints[3] = lonGridIdxHi;

            minDist = 1.e+30;
            snapCrnrPt = 0;

            // Loop through the corners bounding the non-gridded point
            for (crnrPt=0;crnrPt<4;crnrPt++){

                gridLatPt = (int) gridLatPoints[crnrPt];
                gridLonPt = (int) gridLonPoints[crnrPt];

                // Get the lat and lon values at this corner
                gridLatVal = gridLat[nGridCols * gridLatPt + gridLonPt];
                gridLonVal = gridLon[nGridCols * gridLatPt + gridLonPt];

                // The Pythagorean distance of this corner from the non-gridded
                // data point...
                latDist = latVal-gridLatVal;
                lonDist = lonVal-gridLonVal;
                dist = sqrt(latDist*latDist + lonDist*lonDist); 

                // If this distance is the smallest so far, save the corner
                // idx and distance...
                if (dist < minDist){
                    snapCrnrPt = crnrPt;
                    minDist = dist;
                }
            }

            // Get the grid indices of the corner closest to the non-gridded point
            gridLatPt = (int) gridLatPoints[snapCrnrPt];
            gridLonPt = (int) gridLonPoints[snapCrnrPt];

            // Assign the value of the grid point closest to the non-gridded data 
            // point to that non-gridded data point .
            data[idx] = gridData[nGridCols * gridLatPt + gridLonPt];
            /*printf("gridDataIdx indices are (%d, %d)\n",gridLatPt,gridLonPt);*/

            // Save the index of the closest grid data point to the non-gridded 
            // data point.
            gridDataIdx[idx] = nGridCols * gridLatPt + gridLonPt;

            // TODO : Save subsequent data indices in the same gridcell
        }
    }

    return(0);
}
