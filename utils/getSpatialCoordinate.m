function x = getSpatialCoordinate(gridObj)

    x = gridObj.parentGrid.tPFVgeometry.cells.centroids(gridObj.mappings.cellmap);
    x = x(:)';

end
