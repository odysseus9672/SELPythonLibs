#from math import *
from numpy import *
import scipy.interpolate as spintp

#  Change of bases matrices - left multiply on a column vector to get a column
# result. matrix[row][column] x = 0, y = 1, z = 2
deg_per_rad = 180 / pi
rad_per_deg = pi / 180.0
__c2 = 1.0 / 60.0
__c3 = pi / 12.0
TwoPi = 2.0 * pi
HalfPi = 0.5 * pi
ThreeHalfPi = 3.0 * HalfPi

def DMS(deg, min=0.0, sec=0.0):
    """Converts degrees, minutes and seconds of arc into radians."""
    result = rad_per_deg * ( abs(deg) + __c2 * ( abs(min) + __c2 * abs(sec) ) )
    return copysign(result, deg)

def ReadDMS(x, sep=":"):
    """Converts a string containing degrees, minutes and seconds of arc
    into radians."""
    parts = x.split(sep)
    degs = int(parts[0])
    mins = int(parts[1])
    secs = float(parts[2])

    return DMS(degs, min=mins, sec=secs)

def HMS(hr, min=0.0, sec=0.0):
    """Converts hours, minutes and seconds of arc into radians."""
    result = __c3 * ( abs(hr) + __c2 * ( abs(min) + __c2 * abs(sec) ) )
    return copysign(result, hr)

def ReadHMS(x, sep=":"):
    """Converts a string containing hours, minutes and seconds of arc
    into radians."""
    parts = x.split(sep)
    hrs = int(parts[0])
    mins = int(parts[1])
    secs = float(parts[2])

    return HMS(hrs, min=mins, sec=secs)

def RadtoHMS(x, sep=":", prec=3):
    """Converts an arc length in radians to a string containing hours,
    minutes, and seconds separated by the optional arguments sep (defaults to
    ':') and with prec specifying the number of decimal places of precision
    to use in printing the seconds quantity (defaults to 3)."""
    #convert radians to hours
    resid = x / __c3
    hrs = int(resid)
    resid -= hrs
    #dig out minutes
    resid *= 60.0
    mins = int( resid )
    resid -= mins
    #What's left is seconds
    resid *= 60.0

    fmt = sep.join(["{h:0>2d}", "{m:0>2d}",
                    "{s:0>" + str(prec+3) + "." + str(prec) + "f}"])
    #sys.stderr.write(fmt + "\n")
    return( fmt.format( h=hrs, m=mins, s=resid ) )

def RadtoDMS(x, sep=":", prec=2):
    """Converts an arc length in radians to a string containing degrees,
    arcminutes, and arcseconds separated by the optional arguments sep
    (defaults to ':') and with prec specifying the number of decimal places
    of precision to use in printing the seconds quantity (defaults to 2)."""
    #convert radians to degrees > 0
    resid = abs(x) * deg_per_rad
    deg = int(resid)
    resid -= deg
    #convert to arcmin
    resid *= 60
    mins = int(resid)
    resid -= mins
    #What's left is seconds
    resid *= 60

    if ( x >=0 ):
        sgn = "+"
    else:
        sgn = "-"

    fmt = sep.join(["{d:0>2d}", "{m:0>2d}",
                    "{s:0>" + str(prec+3) + "." + str(prec) + "f}"])

    return( sgn + fmt.format( d=deg, m=mins, s=resid ) ) 

def __Transpose(matrix):
    """Produces a new array with the outer two indices switched.
    This implementation is inefficient and slow, and should not be used
    for anything performance critical."""
    result = []
    for i in range(len(matrix[0])):
        new_row = []
        for row in matrix:
            new_row.append(row[i])
        result.append(new_row)
    
    return result

def LonLat_to_UVec(coord):
    """Converts a pair of spherical coordinates, in (longitude, latitude)
    format and radians units, to a unit vector."""
    cd = cos(coord[1]);
    
    result = [0.0, 0.0, 0.0]
    result[0] = cd * cos(coord[0]);
    result[1] = cd * sin(coord[0]);
    result[2] = sin(coord[1]);

    return result

#Rotation matrix for converting to Galactic coordinates from Equatorial
#based on J2000 coordinates from:
#https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm
gal_from_eq = [LonLat_to_UVec( [ HMS(17.76033), DMS(-28.93618) ] ), \
               LonLat_to_UVec( [ HMS(21.20029), DMS(48.32960) ] ), \
               LonLat_to_UVec( [ HMS(12.85730), DMS(27.12830) ] )]

eq_from_gal = __Transpose(gal_from_eq)

#for converting to Galactic from Ecliptic
gal_from_ec = [LonLat_to_UVec( [ DMS(266.83949), DMS(266.83949) ] ), \
               LonLat_to_UVec( [ DMS(347.33989), DMS(59.57417) ] ), \
               LonLat_to_UVec( [ DMS(180.02319), DMS(29.81149) ] )]

ec_from_gal = __Transpose(gal_from_ec)

#for converting to Ecliptic from Equatorial
ec_from_eq = [[ 1.0, 0.0, 0.0 ], \
              [ 0.0,  cos(DMS(23.43929)), sin(DMS(23.43929)) ], \
              [ 0.0,  -sin(DMS(23.43929)), cos(DMS(23.43929)) ]]

eq_from_ec = __Transpose(ec_from_eq)


def Dot_C3V(vec1, vec2):
    """Simple function for computing the dot product between two
    three dimensional vectors (index-able objects of length 3)."""
    return vec1[0] * vec2[0] + vec1[1] * vec2[1] + vec1[2] * vec2[2]

def Length_C3V(vec):
    """Simple function for computing the length of a three dimensional
    vector. Note that any other implementation of a "hypotenuse" function
    will likely perform better."""
    return sqrt(Dot_C3V(vec, vec));

def Hypot_C3V(vec):
    """Computes the hypotenuse of the 3d vector, vec, using a numerically
    optimal algorithm. The argument, vec, must be an indexable object of
    length 3."""
    mags = map(abs, vec)
    mags.sort()
    intermed = 1.0 + (mags[0]/mags[1])**2
    return mags[2] * sqrt(1.0 + intermed * (mags[1]/mags[2])**2)

def UVec_to_LonLat(invec):
    """Convert a unit vector to a (longitude, latitude) pair, in radians."""
    length = Hypot_C3V(invec);
    
    z_norm = invec[2]/length;
    
    alpha = arctan2(invec[1], invec[0]);
    delta = arcsin(z_norm);
    
    if alpha < 0 :
        alpha += TwoPi;
    
    return [alpha, delta]

def Change_Coord(lonlat, ro_matrix):
    """Function for applying a rotation matrix to a (longitude, latitude)
    coordinate pair. The coordinate pair must be in radians."""
    vec = LonLat_to_UVec(lonlat)
    ro_vec = []
    for row in ro_matrix:
        dot = row[0] * vec[0]
        dot += row[1] * vec[1]
        dot += row[2] * vec[2]
        ro_vec.append(dot)

    return UVec_to_LonLat(ro_vec)

def Ec_to_Gal( ec_coord ):
    """Function for converting an Ecliptic coordinate pair to Galactic
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(ec_coord, gal_from_ec)

def Gal_to_Ec( gal_coord ):
    """Function for converting a Galactic coordinate pair to Ecliptic
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(gal_coord, ec_from_gal)

def Eq_to_Gal( eq_coord ):
    """Function for converting an Equatorial coordinate pair to Galactic
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(eq_coord, gal_from_eq)

def Gal_to_Eq( gal_coord ):
    """Function for converting a Galactic coordinate pair to Equatorial
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(gal_coord, eq_from_gal)

def Eq_to_Ec( eq_coord ):
    """Function for converting an Equatorial coordinate pair to Ecliptic
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(eq_coord, ec_from_eq)

def Ec_to_Eq( ec_coord ):
    """Function for converting an Equatorial coordinate pair to Ecliptic
    coordinates. The coordinate pair must be in radians.
    Based on data from:
    https://lambda.gsfc.nasa.gov/toolbox/tb_coordconv.cfm"""
    return Change_Coord(ec_coord, eq_from_ec)

def Sqr_Chord_Len_C3V(vec1, vec2):
    """Computes the square length of the difference between the two arguments.
    The arguments are assumed to be 3d vectors (indexable items of length 3)."""
    diff0 = vec1[0] - vec2[0];
    diff1 = vec1[1] - vec2[1];
    diff2 = vec1[2] - vec2[2];

    result = diff0 * diff0 + diff1 * diff1 + diff2 * diff2;
    return result

def Haversine_LonLat( ra_dec1, ra_dec2 ):
    """Computes the Haversine of the arc connecting the point on a sphere at
    (longitude, latitude) type point ra_dec1 to ra_dec2. The points must
    be in radians."""
    d_lat = ( ra_dec1[1] - ra_dec2[1]) * .5
    d_lon = ( ra_dec1[0] - ra_dec2[0]) * .5

    result = sin(d_lat)**2 + cos(ra_dec1[1]) *  cos(ra_dec2[1]) * \
             sin(d_lon)**2
    return result

def Haversine_Ang( theta ):
    """Computes the Haversine of its argument (argument must be in radians)."""
    return sin(theta * .5)**2

def aHaversine( x ):
    """Computes the inverse Haversine (arcHaversine) of its argument.
    Result is in radians."""
    return 2.0 * arcsin(sqrt(x))

def Arc_Len_Equatorial( ra_dec1, ra_dec2 ):
    """Computes the arc length along the great circle connecting the
    (longitude, latitude) type point ra_dec1 to ra_dec2 using the
    Haversine algorithm (optimal for small arc lengths, sub-optimal for
    large ones). Input coordinates must be in radians."""
    diff_delta = ra_dec1[1] - ra_dec2[1]
    diff_alpha = ra_dec1[0] - ra_dec2[0]
    result = 2.0 * arcsin(sqrt( sin(diff_delta / 2.0)**2 + cos(ra_dec1[1]) *
                                cos(ra_dec2[1]) * sin(diff_alpha / 2.0)**2))
    return result

def ArcLen_LonLat( coord1, coord2 ):
    """Computes the arc length along the great circle connecting the
    (longitude, latitude) type point coord1 to coord2 using the
    Haversine algorithm (optimal for small arc lengths, sub-optimal for
    large ones). Input coordinates must be in radians."""
    diff_lat = coord1[1] - coord2[1]
    diff_lon = coord1[0] - coord2[0]
    result = 2.0 * arcsin(sqrt( sin(diff_lat / 2.0)**2 + cos(coord1[1]) *
                                cos(coord2[1]) * sin(diff_lon / 2.0)**2))
    return result

def ArcLen_LonLat_acos( coord1, coord2 ):
    """Computes the arc length along the great circle connecting the
    (longitude, latitude) type point coord1 to coord2 using the naive
    dot product approach. Most appropriate for arc lengths around pi/2.
    Input coordinates must be in radians."""
    uvec1 = LonLat_to_UVec( coord1 )
    uvec2 = LonLat_to_UVec( coord2 )

    return arccos( Dot_C3V( uvec1, uvec2 ) )

def ArcLen_LonLat_Kahan( coord1, coord2 ):
    """Computes the arc length along the great circle connecting the
    (longitude, latitude) type point coord1 to coord2 using the
    Kahan's formulae (numerical stable for all angle). Input coordinates must
    be in radians.
    See:
    https://scicomp.stackexchange.com/questions/27689/numerically-stable-way-of-computing-angles-between-vectors
    Pg 15 of: https://people.eecs.berkeley.edu/~wkahan/MathH110/Cross.pdf
    """
    uvec1 = LonLat_to_UVec( coord1 )
    uvec2 = LonLat_to_UVec( coord2 )
    sumvec, difvec = zip( *[ (x + y, x - y) for x, y in zip(uvec1, uvec2) ] )

    return 2.0 * arctan2( Hypot_C3V(difvec), Hypot_C3V(sumvec) )

def PA_LonLat( coord1, coord2 ):
    """Calculate the heading angle of the path to coord2 from the position of
    coord1. All quantities in radians."""
    a0, d0 = coord1[0:2]
    a1, d1 = coord2[0:2]

    dalpha = a0 - a1
    cd1 = cos( d1 )
    y = cd1 * sin( dalpha )
    x = - sin( d1 ) * cos( d0 ) + cd1 * sin( d0 ) * cos( dalpha )
    
    angle = arctan2( y, x )
    if angle < 0:
        angle += 2*pi
        
    return angle

def DeltaLon( lonlat1, lonlat2 ):
    """Calculates the difference in longitude between to points. All
    quantities in radians."""
    diff = lonlat2[0] - lonlat1[0]
    sgn = copysign( 1.0, diff )
    diff = abs( diff )
    return sgn * min( diff, 2*pi - diff )

def ArcLenSqrApprox1_LonLat( coord1, coord2 ):
    """Calculates the square arc length between (longitude, latitude) type
    points coord1 and coord2, subject to the approximation that the
    sine of the arc length is small. All quantites in radians."""
    diff_lat = coord1[1] - coord2[1]
    diff_lon = coord1[0] - coord2[0]
    result = diff_lat**2 + (cos(coord1[1]) * diff_lon)**2
    return result

def ArcLenApprox1_LonLat( coord1, coord2 ):
    """Calculates the arc length between (longitude, latitude) type
    points coord1 and coord2, subject to the approximation that the
    sine of the arc length is small. All quantites in radians."""
    diff_lat = coord1[1] - coord2[1]
    diff_lon = coord1[0] - coord2[0]
    result = diff_lat**2 + (cos(coord1[1]) * diff_lon)**2
    return sqrt(result)

def ArcLenApprox2_LonLat( coord1, coord2 ):
    """Calculates an upper bound to the arc length between
    (longitude, latitude) type points coord1 and coord2 using a Manhatten
    style approximation. The longitude is traversed as near the equator as
    possible. All quantities in radians."""
    diff_lon = abs( coord1[0] - coord2[0] )
    diff_lon = min( diff_lon, 2*pi - diff_lon )
    diff_lat = abs( coord1[1] - coord2[1] )

    c = cos( min( abs(coord[0]), abs(coord[1]) ) )

    return diff_lat + c * diff_lon

def MaxLonSep( maxarc, baselat ):
    """Calculates the maximum separation in longitude that a point can have
    from a reference point at latitude baselat and still be within a given
    great circle arc length, maxarc, of the reference point. All quantities
    in radians."""
    if abs(baselat) + maxarc <= 0.5 * pi:
        #result = asin( abs( sin(maxarc) ) / cos( baselat ) )
        #result = acos(sqrt(cos(baselat)**2 - sin(maxarc)**2)/cos(baselat))
        c = cos( baselat )
        s = abs( sin( maxarc ) )
        y = s 
        x = sqrt( ( c + s ) * ( c - s ) )
        result = atan2( y, x )
    else:
        result = pi
    return result

def MaxLonSepQuad( maxarc, baselat ):
    """Calculates the same quantity as MaxLonSep, but provides an upper
    bound based on a quadratic approximation instead of calling the exact
    trig functions. All quantites in radians."""
    if abs(baselat) + maxarc <= 0.5 * pi:
        result =  baselat**2 / ( 0.5 * pi - maxarc )
        #Ensure upper limit when pi/2 - maxarc = pi/2
        result *= 1. + 2.0**-52 
        result += maxarc
    else:
        result = pi
    return result

def MaxLonSepLin( maxarc, baselat ):
    """Calculates the same quantity as MaxLonSep, but provides an upper
    bound based on a linear approximation instead of calling the exact trig
    functions. All quantities in radians."""
    if abs(baselat) + maxarc <= 0.5 * pi:
        result = maxarc + abs(baselat)
    else:
        result = pi
    return result

def Cross_C3V(vec1, vec2, result):
    """Computes the cross product between the 3d vectors vec1 and vec2
    (indexable objects of length 3) and stores the result in the third
    argument, named result."""
    x1, y1, z1 = vec1;
    x2, y2, z2 = vec2;

    result[0] = y1 * z2 - z1 * y2;
    result[1] = z1 * x2 - x1 * z2;
    result[2] = x1 * y2 - y1 * x2;

    return None

def Arc_Len_C3V(vec1, vec2):
    """Computes the arc length of the great circle connecting unit
    vectors vec1 and vec2."""
    c_theta = Dot_C3V(vec1, vec2)
    s_theta_vec = [0.0, 0.0, 0.0]
    Cross_C3V(vec1, vec2, s_theta_vec)
    s_theta = Hypot_C3V(s_theta_vec)

    theta = arctan2(s_theta, c_theta)
    #s_theta >= 0 so theta > 0 always
##    if theta < 0.0:
##        theta = abs(theta)

    return theta

def Normalize_C3V(vec):
    """Scales the argument by its Euclidean length, assuming its a 3d vector
    (indexable object of length 3)."""
    length = Hypot_C3V(vec);
    vec[0] /= length
    vec[1] /= length
    vec[2] /= length

    return None;

def Rotate_C3V(vec, axis, theta, result):
    """Rotates the 3d vector, vec (arg 1), around the 3d vector in
    axis (arg 2) by the angle theta (arg 3, in radians), and stores the
    result in the 3d vector result (arg 4). vec is unchanged, and all
    3d vectors are length 3 indexable ojbects."""
    vec3 = [0.0, 0.0, 0.0];

    coef1 = cos(theta);

    coef2 = sin(theta * 0.5);
    coef2 *= 2.0 * coef2 * Dot_C3V(axis, vec);

    coef3 = sin(theta);

    Cross_C3V(axis, vec, vec3);

    result[0] = coef1 * vec[0] + coef2 * axis[0] + coef3 * vec3[0];
    result[1] = coef1 * vec[1] + coef2 * axis[1] + coef3 * vec3[1];
    result[2] = coef1 * vec[2] + coef2 * axis[2] + coef3 * vec3[2];

    return None

def Orthog_Project_C3V(projectee, ortho, result):
    """Computes the orthogonal projection of the 3d vector, projectee
    (arg 1), along the unit 3d vector, ortho (arg 2), and stores the answer
    in result (arg 3). All 3d vectors must be length 3 indexable ojbects,
    and ortho must have unit length."""
    intermediate = [0.0, 0.0, 0.0];
    Cross_C3V(projectee, ortho, intermediate);
    Cross_C3V(ortho, intermediate, result);
    
    return None

def Pick_Orthog_C3V(original, result):
    """Picks a 3d vector orthogonal to 3d unit vector original (arg 1) and
    stores it in the 3d vector result (arg 2). The algorithm is to take a
    basis vector along original's smallest component, calculate its
    orthogonal projection along original, and then normalize the projected
    vector. The vector in original must have length 1. All 3d vectors must
    be length 3 indexable ojbects."""
    mag = map(abs, original)
    smallest = min(mag)
    result[0] = result[1] = result[2] = 0.0
    result[mag.index(smallest)] = 1.0

    Orthog_Project_C3V(result, original, result)
    Normalize_C3V(result)

    return None

def Solid_Angle(arc_radius):
    """Computes the solid angle subtended by a circular section of a unit
    sphere with the given arc radius (in radians). Result is in steradians."""
    return 4.0 * pi * sin(arc_radius / 2.0)**2

def __LonShift( lon, lon0 ):
    #Map longitudes to -pi to pi
    if lon > pi:
        lon -= TwoPi

    #Shift the longitudes and remap to -pi to pi
    lon -= lon0
    if lon >= TwoPi or lon <= TwoPi:
        lon -= floor( lon / TwoPi ) * TwoPi

    if lon > pi:
        lon -= TwoPi
    elif lon <= -pi:
        lon += TwoPi

    return lon

#Implement the Mollweide projection
#http://en.wikipedia.org/wiki/Mollweide_projection
__t0 = linspace( 0.0, .5*pi, 2048 )
__MollweideThetaDat = ( arcsin( (2.0*__t0 + sin(2.0*__t0)) / pi ), __t0 )
__MollweideTheta = spintp.UnivariateSpline( __MollweideThetaDat[0],
                                            __MollweideThetaDat[1],
                                            k=3, s=0 )
__MwCx = 2.0 * sqrt(2.0) / pi
__MwCy = sqrt(2.0)

def MollweideProj( lonlat, lon0=0.0 ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the Mollweide
    projection of that point. lon0 is the central longitude for the
    projection process (in radians). The x coordinate is in the range
    [-2*\sqrt(2), 2*sqrt(2)] and the y coordinate is in the range
    [-sqrt(2), sqrt(2)]."""
    lon, lat = lonlat
    lon = __LonShift( lon, lon0 )
    
    t = __MollweideTheta( abs(lat) )
    t = copysign( t, lat )
        
    return ( __MwCx * lon * cos(t), __MwCy * sin(t) )

#Implement the Lambert azimuthal equal-area projection
#http://en.wikipedia.org/wiki/Lambert_azimuthal_equal-area_projection

__LamConst = 2.0

def LambertAzEqAreaProj( lonlat ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the Lambert azimuthal
    equal-area projection of that point, with the center of projection being
    the North pole of the lonlat system. The x and y
    coordinates are contained in a circle of radius 2.0."""
    lon, lat = lonlat
    
    chalfpolang = cos(0.5*(HalfPi + lat))
    return ( __LamConst * cos(lon) * chalfpolang,
             __LamConst * sin(lon) * chalfpolang )
    

#Implement the Hammer projection (often incorrectly called Aitoff)
#http://en.wikipedia.org/wiki/Hammer_projection

__HmCx = 2.0 * sqrt(2.0)
__HmCy = sqrt(2.0)

def HammerProj( lonlat, lon0=0.0 ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the Hammer
    projection of that point. lon0 is the central longitude for the
    projection process (in radians). The x coordinate is in the range
    [-2*\sqrt(2), 2*sqrt(2)] and the y coordinate is in the range
    [-sqrt(2), sqrt(2)]."""
    lon, lat = lonlat
    lon = __LonShift( lon, lon0 )
    
    clat = cos(lat)
    shalflon = sin(0.5 * lon)
    chalflon = cos(0.5 * lon)
    norm = 1.0 / sqrt( 1.0 + clat * chalflon )
    
    return ( __HmCx * clat * shalflon * norm, __HmCy * sin(lat) * norm )

#Implement the Azimuthal equidistant projection
#http://en.wikipedia.org/wiki/Azimuthal_equidistant_projection
def AzEqDistProj( lonlat ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the azimuthal
    equidistant projection of that point, with the center of projection being
    the North pole of the lonlat system. The x and y
    coordinates are contained in a circle of radius pi."""
    lon, lat = lonlat
    polang = HalfPi - lat
    
    return( polang * cos(lon), polang * sin(lon)  )

#Implement the Aitoff projection
#http://en.wikipedia.org/wiki/Aitoff_projection
def __invsincu( x ): #Unnormalized multiplicative inverse of the sinc function
    if x != 0.0:
        return x / sin(x)
    else:
        return 1.0
    
def AitoffProj( lonlat, lon0=0.0 ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the Aitoff
    projection of that point. lon0 is the central longitude for the
    projection process (in radians). The x coordinate is in the range
    [-pi, pi] and the y coordinate is in the range [-pi/2, pi/2]."""
    lon, lat = lonlat
    lon = __LonShift( lon, lon0 )
    
    clat = cos(lat)
    invSincAlpha = __invsincu(arccos( clat * cos( 0.5*lon ) ))

    return( 2.0 * clat * sin( 0.5*lon ) * invSincAlpha,
            sin( lat ) * invSincAlpha )

#Implement the Winkel tripel projection
#http://en.wikipedia.org/wiki/Winkel_Tripel_projection

def WinkelTripelProj( lonlat, StandardLat=arccos(2.0/pi), lon0=0.0 ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the Winkel tripel
    projection of that point. lon0 is the central longitude for the
    projection process (in radians). With the default settings, the x
    coordinate is in the range (-pi/2 - 1, pi/2 + 1) and the y coordinate is
    in the range [-pi/2, pi/2]. Modifying StandardLat can change this."""
    lon, lat = lonlat
    lon = __LonShift( lon, lon0 )
    
    clat = cos(lat)
    isinc = __invsincu( arccos( clat * cos(0.5*lon) ) )
    
    return( 0.5 * ( lon * cos(StandardLat)
                    + 2.0 * clat * sin(0.5 * lon) * isinc ),
            0.5 * ( lat + sin(lat) * isinc ) )

#Implement Cylindrical equal area projection
#http://en.wikipedia.org/wiki/Cylindrical_equal-area_projection
def CylindEqualAreaProj( lonlat, StandardLat=0.25*pi, lon0=0.0 ):
    """Takes a point at the (longitude, latitude) point in lonlat (arg 1,
    in radians) and returns the (x,y) coordinates of the cylindrical equal area
    projection of that point. lon0 is the central longitude for the
    projection process (in radians). With the default settings, the x
    coordinate is in the range
    [-pi*cos(StandardLat), pi*cos(StandardLat)) and the y coordinate is
    in the range [-1/cos(StandardLat), 1/cos(StandardLat)]"""
    lon, lat = lonlat
    lon = __LonShift( lon, lon0 )

    clat0 = cos(StandardLat)
    return( lon * clat0, sin(lat) / clat0 )
