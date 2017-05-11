"""
Cython kernels for the fatiando.gravmag.tesseroid module.

Used to optimize some slow tasks and compute the actual gravitational fields.
"""
from __future__ import division
import numpy

from .. import constants

from libc.math cimport sin, cos, sqrt, atan2, acos
# Import Cython definitions for numpy
cimport numpy
cimport cython

# To calculate sin and cos simultaneously
cdef extern from "math.h":
    void sincos(double x, double* sinx, double* cosx)

cdef:
    double d2r = numpy.pi/180.
    double[::1] nodes
    double MEAN_EARTH_RADIUS = constants.MEAN_EARTH_RADIUS
    double G = constants.G
    ctypedef double (*kernel_func)(double, double, double, double, double,
                                   double[::1], double[::1], double[::1],
                                   double[::1])
nodes = numpy.array([-0.577350269189625731058868041146,
                     0.577350269189625731058868041146])


@cython.boundscheck(False)
@cython.wraparound(False)
cdef inline double scale_nodes(
    double w, double e, double s, double n, double top, double bottom,
    double[::1] lonc,
    double[::1] sinlatc,
    double[::1] coslatc,
    double[::1] rc):
    "Put GLQ nodes in the integration limits for a tesseroid"
    cdef:
        double dlon, dlat, dr, mlon, mlat, mr, scale, latc
        unsigned int i
    dlon = e - w
    dlat = n - s
    dr = top - bottom
    mlon = 0.5*(e + w)
    mlat = 0.5*(n + s)
    mr = 0.5*(top + bottom + 2.*MEAN_EARTH_RADIUS)
    # Scale the GLQ nodes to the integration limits
    for i in range(2):
        lonc[i] = d2r*(0.5*dlon*nodes[i] + mlon)
        latc = d2r*(0.5*dlat*nodes[i] + mlat)
        sinlatc[i] = sin(latc)
        coslatc[i] = cos(latc)
        rc[i] = (0.5*dr*nodes[i] + mr)
    scale = d2r*dlon*d2r*dlat*dr*0.125
    return scale


@cython.boundscheck(False)
@cython.wraparound(False)
def potential(
    numpy.ndarray[double, ndim=1] bounds,
    double density,
    double ratio,
    int stack_max,
    numpy.ndarray[double, ndim=1] lons,
    numpy.ndarray[double, ndim=1] sinlats,
    numpy.ndarray[double, ndim=1] coslats,
    numpy.ndarray[double, ndim=1] radii,
    numpy.ndarray[double, ndim=1] result):
    """
    Calculate this gravity field of a tesseroid at given locations.
    """
    with_rediscretization(bounds, density, ratio, stack_max, lons, sinlats,
                          coslats, radii, result, kernelpot)


@cython.boundscheck(False)
@cython.wraparound(False)
@cython.cdivision(True)
cdef inline double kernelpot(
    double lon, double sinlat, double coslat, double radius, double scale,
    double[::1] lonc, double[::1] sinlatc, double[::1] coslatc,
    double[::1] rc):
    cdef:
        unsigned int i, j, k
        double kappa, r_sqr, coslon, l_sqr, cospsi
        double result
    r_sqr = radius**2
    result = 0
    for i in range(2):
        coslon = cos(lon - lonc[i])
        for j in range(2):
            cospsi = sinlat*sinlatc[j] + coslat*coslatc[j]*coslon
            for k in range(2):
                l_sqr = r_sqr + rc[k]**2 - 2*radius*rc[k]*cospsi
                kappa = (rc[k]**2)*coslatc[j]
                result += kappa/sqrt(l_sqr)
    return result*scale


@cython.boundscheck(False)
@cython.wraparound(False)
cdef with_rediscretization(
    numpy.ndarray[double, ndim=1] bounds,
    double density,
    double ratio,
    int stack_max,
    numpy.ndarray[double, ndim=1] lons,
    numpy.ndarray[double, ndim=1] sinlats,
    numpy.ndarray[double, ndim=1] coslats,
    numpy.ndarray[double, ndim=1] radii,
    numpy.ndarray[double, ndim=1] result,
    kernel_func kernel):
    """
    Calculate the given kernel function on the computation points by
    rediscretizing the tesseroid when needed.
    """
    cdef:
        unsigned int l, size
        double[::1] lonc =  numpy.empty(2, numpy.float)
        double[::1] sinlatc =  numpy.empty(2, numpy.float)
        double[::1] coslatc =  numpy.empty(2, numpy.float)
        double[::1] rc =  numpy.empty(2, numpy.float)
        double[::1] bounds
        double scale
        unsigned int i, j, k
        int nlon, nlat, nr
        int stack_size
        double[:, ::1] stack =  numpy.empty((stack_max, 6), numpy.float)
        double w, e, s, n, top, bottom
        double lon, sinlat, coslat, radius
        double dlon, dlat, dr, distance
        double res
    size = len(result)
    for l in range(size):
        lon = lons[l]
        sinlat = sinlats[l]
        coslat = coslats[l]
        radius = radii[l]
        res = 0
        for i in range(6):
            stack[0, i] = bounds[i]
        stack_size = 0
        while stack_size > 0:
            stack_size = stack_size - 1
            w = stack[stack_size, 0]
            e = stack[stack_size, 1]
            s = stack[stack_size, 2]
            n = stack[stack_size, 3]
            top = stack[stack_size, 4]
            bottom = stack[stack_size, 5]
            distance = distance_n_size(w, e, s, n, top, bottom, lon, sinlat,
                                       coslat, radius, &dlon, &dlat, &dr)
            # Check which dimensions I have to divide
            nlon = 1
            nlat = 1
            nr = 1
            if distance < ratio*dlon:
                nlon = 2
            if distance < ratio*dlat:
                nlat = 2
            if distance < ratio*dr:
                nr = 2
            if nlon == 1 and nlat == 1 and nr == 1:
                # Put the nodes in the current range
                scale = scale_nodes(w, e, s, n, top, bottom, lonc, sinlatc,
                                    coslatc, rc)
                res += kernel(lon, sinlat, coslat, radius, scale, lonc,
                              sinlatc, coslatc, rc)
            else:
                if stack_size + nlon*nlat*nr > stack_max:
                    raise ValueError('Tesseroid stack overflow')
                dlon = (e - w)/nlon
                dlat = (n - s)/nlat
                dr = (top - bottom)/nr
                for i in xrange(nlon):
                    for j in xrange(nlat):
                        for k in xrange(nr):
                            stack[stack_size, 0] = w + i*dlon
                            stack[stack_size, 1] = w + (i + 1)*dlon
                            stack[stack_size, 2] = s + j*dlat
                            stack[stack_size, 3] = s + (j + 1)*dlat
                            stack[stack_size, 4] = bottom + (k + 1)*dr
                            stack[stack_size, 5] = bottom + k*dr
                            stack_size = stack_size + 1
        result[l] += density*res


@cython.cdivision(True)
cdef inline double distance_n_size(
    double w, double e, double s, double n, double top, double bottom,
    double lon, double sinlat, double coslat, double radius,
    double* dlon, double* dlat, double* dr):
    cdef:
        double rt, rtop, lont, latt, sinlatt, coslatt, cospsi, distance
    # Calculate the distance to the observation point
    rt = 0.5*(top + bottom) + MEAN_EARTH_RADIUS
    lont = d2r*0.5*(w + e)
    latt = d2r*0.5*(s + n)
    sinlatt = sin(latt)
    coslatt = cos(latt)
    cospsi = sinlat*sinlatt + coslat*coslatt*cos(lon - lont)
    distance = sqrt(radius**2 + rt**2 - 2*radius*rt*cospsi)
    # Calculate the dimensions of the tesseroid in meters
    rtop = top + MEAN_EARTH_RADIUS
    dlon[0] = rtop*acos(sinlatt**2 + (coslatt**2)*cos(d2r*(e - w)))
    dlat[0] = rtop*acos(sin(d2r*n)*sin(d2r*s) + cos(d2r*n)*cos(d2r*s))
    dr[0] = top - bottom
    return distance
