"""`semimajoraxis(x, y, z, u, v, w, m1, m2)`

Calculates semimajor axis for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function semimajoraxis(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    return 1/( (2/r)-(vsq/M) )
end

"""`eccentricity(x, y, z, u, v, w, m1, m2)`

Calculates eccentricity for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies
with masses `m1` and `m2`."""
function eccentricity(x, y, z, u, v, w, m1, m2)
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    a = semimajoraxis(x, y, z, u, v, w, m1, m2)
    M = m1+m2
    quotient = hsq/( M*a )
    return sqrt(1 - quotient)
end

"""`inclination(x, y, z, u, v, w)`

Calculates inclination for the two
body problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two bodies."""
function inclination(x, y, z, u, v, w)
    hz = x*v-y*u
    h = sqrt( (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2 )
    return acos(hz/h)
end

"""`aei(x, y, z, u, v, w, m1, m2)`

Returns semimajor axis `a`, eccentricity `e` and inclination `I` for the two-body
problem defined by the relative position `(x,y,z)` and velocity `(u,v,w)`vectors
between two bodies with masses `m1` and `m2`."""
function aei(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    a = 1/( (2/r)-(vsq/M) )
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    h = sqrt(hsq)
    quotient = hsq/( M*a )
    e = sqrt(1 - quotient)
    hz = x*v-y*u
    inc = acos(hz/h)
    return a, e, inc
end

"""`ae(x, y, z, u, v, w)`

Returns semimajor axis `a` and eccentricity `e` for the two-body problem defined
by the relative position `(x,y,z)` and velocity `(u,v,w)` vectors between two
bodies with masses `m1` and `m2`."""
function ae(x, y, z, u, v, w, m1, m2)
    r = sqrt(x^2+y^2+z^2)
    vsq = u^2+v^2+w^2
    M = m1+m2
    a = 1/( (2/r)-(vsq/M) )
    hsq = (y*w-z*v)^2 + (z*u-x*w)^2 + (x*v-y*u)^2
    quotient = hsq/( M*a )
    e = sqrt(1 - quotient)
    return a, e
end
