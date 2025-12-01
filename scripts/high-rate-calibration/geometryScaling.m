function y = geometryScaling(geom, x)

    charLength = sqrt((geom.L/2)^2 + (geom.h/2)^2);
    factor = geom.t / charLength;
    y = factor * x;

end
