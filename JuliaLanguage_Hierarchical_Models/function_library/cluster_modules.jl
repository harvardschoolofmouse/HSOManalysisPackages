#
#  The following functions allow us to get distances between points and such
#
#   getEuclideanDistance(pos1::Tuple, pos2::Tuple)
#   isclose(pos1, pos2, threshold)
#

function getEuclideanDistance(pos1::Tuple, pos2::Tuple)
    x1 = pos1[1];
    x2 = pos2[1];
    y1 = pos1[2];
    y2 = pos2[2];
    distance = sqrt((x2-x1)^2 + (y2-y1)^2)
    return distance 
end;
function isclose(pos1, pos2, threshold)
    d = getEuclideanDistance(pos1, pos2)
    if d < threshold
        return true
    else
        return false
    end
end;