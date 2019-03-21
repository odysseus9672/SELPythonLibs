# Libary that supplies numerical integration routines for higher dimensional
# spaces. The "Simple" functions are basically the trapezoid rule, and the 
# "Romberg" functions use the Romberg algorithm to extrapolate to higher 
# orders than the trapezoid rule.

import math
#import sys
from NevInterp import NevilleEvenRound



class BugError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)

class NotImplementedError(Exception):
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


def __MkAccumulator( function, weight, divs ):
	def Acc( tot, idx ):
		f = function
		w = weight
		d = float(divs)
		return tot + w * f( float(idx) / d )

	return Acc


def Romberg3dO1( fun ):
	"""Performs 3 dimensional Romberg integration based on an order 1
	estimator. The order one estimator is to divide the region of integration
	into cubic cells and each cell contributes fun(center) * h**3 to the
	integral, where h is the side length of the cube. The region of
	integration assumed is the octant 1 unit
	cube: x in [0,1], y in [0,1], and z in [0,1]. The user should
	perform the needed changes of variables to make this so. The integrand
	must also be finite in the entire region. """

	#Note: The sequences increase by a factor of 3 because then the samples from
	#      the smaller number of divs are a subset of the higher ones.
	#      For this function the individual sequences must be in increasing order.
	divseqs = [ ( 9, ), ( 6, 18 ), ( 7, 21 ), ( 12, ), ( 5, 15 ) ]

	estimates = []
	scales = []
	funcalls = 0
	for s in divseqs:
		cur_ests = [ 0.0 for x in range( len(s) ) ]
		bigdivs = s[-1]
		cellwidth = 1.0 / bigdivs
		startval = 0.5 * cellwidth

		cyclic_maxes = [ bigdivs / x for x in s[:-1] ]

		cyclic_idxes = [ [ y/2 + 1, y/2 + 1, y/2 + 1 ] for y in cyclic_maxes ]
		weights = [ 1.0 / (x * x * x) for x in s ]

		coords = [ startval, startval, startval ]
		for zi in range( bigdivs ):

			coords[1] = startval
			for i in range( len(cyclic_idxes) ):
				cyclic_idxes[i][1] = cyclic_maxes[i] / 2 + 1

			for yi in range( bigdivs ):

				coords[0] = startval
				for i in range( len(cyclic_idxes) ):
					cyclic_idxes[i][0] = cyclic_maxes[i] / 2 + 1

				for xi in range( bigdivs ):
					val = fun( coords[0], coords[1], coords[2] )
					funcalls += 1
					cur_ests[-1] += val * weights[-1]

					coords[0] += cellwidth
					for i in range( len(cyclic_idxes) ):
						if all([ x == 0 for x in cyclic_idxes[i] ]):
							cur_ests[i] += val * weights[i]

						cyclic_idxes[i][0] += 1
						if cyclic_idxes[i][0] >= cyclic_maxes[i]:
							cyclic_idxes[i][0] = 0

				coords[1] += cellwidth
				for i in range( len(cyclic_idxes) ):
					cyclic_idxes[i][1] += 1

					if cyclic_idxes[i][1] >= cyclic_maxes[i]:
						cyclic_idxes[i][1] = 0

			coords[2] += cellwidth
			for i in range( len(cyclic_idxes) ):
				cyclic_idxes[i][2] += 1

				if cyclic_idxes[i][2] >= cyclic_maxes[i]:
					cyclic_idxes[i][2] = 0

		estimates.extend( cur_ests )
		scales.extend( [ 1.0 / (x*x*x) for x in s ] )

	sortinglist = zip( scales, estimates )
	sortinglist.sort( key= lambda x: x[0] )
	scales, estimates = zip( *sortinglist )
	#print(funcalls)
	return NevilleEvenRound( scales, estimates, 0.0 )


def Romberg3dO3( fun ):
	"""Performs 3 dimensional Romberg integration based on an order 3
	estimator. The order three estimator is to divide the region of integration
	into cubic cells and each cell contributes
	(2<f(cell edge centers)> - <f(cell corners)>) * h**3 to the
	integral, where h is the side length of the cube. The region of
	integration assumed is the octant 1 unit
	cube: x in [0,1], y in [0,1], and z in [0,1]. The user should
	perform the needed changes of variables to make this so. The integrand
	must also be finite in the entire region. """
	raise NotImplementedError("Romberg3dO3 not implemented yet. Probably not ever.")
	#Note: The sequences increase by a factor of 2 because then samples from
	#      the smaller number of divs are a subset of the higher ones.
	#      For this function, the individual sequences must be in decreasing order.
	divseqs = [ ( 2, ), ( 12, ), ( 6, ), ( 3, ), ( 5, ), (10,), ( 9,), (18,) ]

	estimates = []
	scales = []
	funcalls = 0
	for s in divseqs:
		cur_ests = [ 0.0 for x in range( len(s) ) ]
		maingridsize = 2 * s[0] + 1
		edgeidx = maingridsize - 1
		stepsize = 0.5 / s[0]

		base_weights = [ 1.0 / (x * x * x) for x in s ]
		indicators = [ [ 0, 0, 0 ] for x in s ]

		coords = [ 0.0, 0.0, 0.0 ]
		for zi in range( maingridsize ):
			for i in range(len(indicators)):
				indicators[i][2] = zi & (2**i)

			coords[1] = 0.0
			for yi in range( maingridsize ):
				for i in range(len(indicators)):
					indicators[i][1] = yi & (2**i)

				coords[0] = 0.0
				for xi in range( maingridsize ):
					val = None
					lasttp = 0

					for i in range(len(indicators)):
						indicators[i][0] = xi & (2**i)

						tp = sum( [ x / (2**i) for x in indicators[i] ] )
						if lasttp != 0:
							break

						lasttp = tp
						#tp == 0 => corner
						#tp == 1 => edge center
						#tp == 2 => face center
						#tp == 3 => center
						locs = [ x > 0 and x < edgeidx for x in (xi, yi, zi) ]
						locs.sort()

						#decide the weight, or if the loop can exit early
						if tp == 1: #Cell edge
							#interior
							if all(locs):
								w = 2.0/3

							#face
							elif all(locs[1:]):
								w = 1.0/3

							#edge
							elif locs[2]:
								w = 1.0/6

							else: #This would be a cell edge on the region corner
								raise BugError("You shouldn't see this unless there was a bug or data corruption.")

						elif tp == 0: #Cell corner
							#interior
							if all(locs):
								w = -1.0

							#face
							elif all(locs[1:]):
								w = -0.5

							#edge
							elif locs[2]:
								w = -0.25

							#corner
							else:
								w = -0.125

						else: #Point is on a cell face or center, not relevant anymore
							break

						if val == None:
							val = fun( coords[0], coords[1], coords[2] )
							funcalls += 1

						cur_ests[i] += val * base_weights[i] * w

					coords[0] += stepsize

				coords[1] += stepsize

			coords[2] += stepsize

		estimates.extend( cur_ests )
		scales.extend( [ 1.0 / (x*x*x) for x in s ] )

	#print( scales )
	#print( estimates )

	sortinglist = zip( scales, estimates )
	sortinglist.sort( key= lambda x: x[0] )
	scales, estimates = zip( *sortinglist )
	print(funcalls)
	return NevilleEvenRound( scales, estimates, 0.0 )

def Simple2dO1( fun, divs ):
	"""Performs simple 2 dimensional trapezoidal integration on fun on the
	closed unit rectangle: x in [0,1], y in [0,1]. Divs must be a tuple/list
	of positive integers with length 2 and it sets the number of divisions in each
	dimension."""
	if len(divs) != 2:
		raise IndexError( "Simple2dO1: divs must have length 2." )

	divs = tuple([ int(d) for d in divs ])

	if divs[0] <= 0 or divs[1] <= 0:
		raise ValueError( "Simple2dO1: elements in divs must be > 0." )

	xdivs, ydivs = divs
	basearea = 1.0 / (xdivs * ydivs)

	#Add corners
	result = 0.25 * basearea * ( fun(0.0, 0.0) + fun(0.0, 1.0) +
								 fun(1.0, 0.0) + fun(1.0, 1.0) )

	#Add x-like sides, using a function closure to reduce lookup times for symbols in the reduce
	def MkAcc( function, weight, divs ):
		def Acc( tot, xi ):
			f = function
			w = weight
			d = float(divs)
			x = float(xi)/divs
			return tot + w * (f(x, 0.0) + f(x, 1.0))
		return Acc

	result += reduce( MkAcc( fun, 0.5 * basearea, xdivs ), xrange( 1, xdivs ), 0.0 )

	#Add y-like sides, using a function closure to reduce lookup times for symbols in the reduce
	def MkAcc( function, weight, divs ):
		def Acc( tot, yi ):
			f = function
			w = weight
			d = float(divs)
			y = float(yi)/divs
			return tot + w * (f(0.0, y) + f(1.0, y))
		return Acc

	result += reduce( MkAcc( fun, 0.5 * basearea, ydivs ), xrange( 1, ydivs ), 0.0 )

	#Add center region:
	def MkAcc( function, weight, xval, ydivs ):
		def Acc( tot, yi ):
			f = function
			w = weight
			d = float(ydivs)
			x = xval
			return tot + w * f(x, float(yi)/d)
		return Acc

	for xi in xrange( 1, xdivs):
		x = float(xi) / xdivs

		result += reduce( MkAcc( fun, basearea, x, ydivs ), xrange( 1, ydivs ), 0.0 )

	return result

def Simple2dO2( fun, divs ):
	"""Performs simple 2 dimensional second order integration on fun on the
	closed unit rectangle: x in [0,1], y in [0,1]. The user should
	perform the needed changes of variables to make this so. The integrand
	must also be finite in the entire region. Divs must be a tuple/list
	of integers (>= 5) with length 2 and it sets the number of divisions in each
	dimension. Based on equation 4.1.14 of 'Numerical Recipes' Third Edition."""
	raise NotImplementedError("Simple2dO2 not implemented yet. Possibly not ever.")
	if len(divs) != 2:
		raise IndexError( "Simple2dO2: divs must have length 2." )

	divs = tuple([ int(d) for d in divs ])

	if divs[0] < 5 or divs[1] < 5:
		raise ValueError( "Simple2dO2: elements in divs must be >= 5." )

	wx = 1.0 / divs[0]
	wy = 1.0 / divs[1]
	basearea = wx * wy

	#Border values
	xb = ( 0.0, wx, 2.0*wx, 1.0 - 2.0*wy, 1.0 - wx, 1.0 )
	yb = ( 0.0, wy, 2.0*wy, 1.0 - 2.0*wy, 1.0 - wy, 1.0 )

	#Border weghts
	bw = ( 3.0 / 8, 7.0 / 6, 23.0 / 24 )

	#Add corner regions
	result = bw[0]*bw[0] * basearea * ( fun(xb[ 0], yb[0]) + fun(xb[ 0], yb[-1]) +
										fun(xb[-1], yb[0]) + fun(xb[-1], yb[-1]) )

	result += bw[0]*bw[1] * basearea * ( fun(xb[ 0], yb[ 1]) + fun(xb[ 1], yb[ 0]) +
										 fun(xb[ 0], yb[-2]) + fun(xb[ 1], yb[-1]) +
										 fun(xb[-2], yb[ 0]) + fun(xb[-1], yb[ 1]) +
										 fun(xb[-1], yb[-2]) + fun(xb[-2], yb[-1]) )

	result += bw[0]*bw[2] * basearea * ( fun(xb[ 0], yb[ 2]) + fun(xb[ 2], yb[-1]) +
										 fun(xb[ 0], yb[-3]) + fun(xb[ 2], yb[-1]) +
										 fun(xb[-3], yb[ 0]) + fun(xb[-1], yb[ 2]) +
										 fun(xb[-1], yb[-3]) + fun(xb[-3], yb[-1]) )

	result += bw[1]*bw[1] * basearea * ( fun(xb[ 1], yb[ 1]) + fun(xb[ 1], yb[-2]) +
										 fun(xb[-1], yb[-2]) + fun(xb[-2], yb[-2]) )

	result += bw[1]*bw[2] * basearea * ( fun(xb[ 1], yb[ 2]) + fun(xb[ 2], yb[ 1]) +
										 fun(xb[ 1], yb[-3]) + fun(xb[ 2], yb[-2]) +
										 fun(xb[-3], yb[ 1]) + fun(xb[-2], yb[ 2]) +
										 fun(xb[-2], yb[-3]) + fun(xb[-3], yb[-2]) )

	result += bw[2]*bw[2] * basearea * ( fun(xb[ 2], yb[ 2]) + fun(xb[ 2], yb[-3]) +
										 fun(xb[-3], yb[ 2]) + fun(xb[-3], yb[-3]) )

	#Add side stripes
	for xi in xrange(3, divs[0]-2):
		x = xi * wx
		result += basearea * ( bw[0] * (fun(x, yb[0]) + fun(x, yb[-1])) +
							   bw[1] * (fun(x, yb[1]) + fun(x, yb[-2])) +
							   bw[2] * (fun(x, yb[2]) + fun(x, yb[-3])) )

	for yi in xrange(3, divs[1]-2):
		y = yi * wy
		result += basearea * ( bw[0] * (fun(xb[0], y) + fun(xb[-1], y)) +
							   bw[1] * (fun(xb[1], y) + fun(xb[-2], y)) +
							   bw[2] * (fun(xb[2], y) + fun(xb[-3], y)) )

	#Add center region
	for xi in xrange( 4, divs[0]-3 ):
		x = xi * wx
		for yi in xrange( 4, divs[1]-3 ):
			result += basearea * fun( x, yi * wy )

	return result


def Romberg2dO1( fun, coarsedivs=( (4,4), (5,5) ), addlevels=( 3, 3 ) ):
	"""Performs 2 dimensional Romberg integration based on an order 1
	estimator. The order one estimator is to divide the region of integration
	into square cells and each cell contributes 0.25 * sum(fun(corners)) * hx * hy
	to the integral, where hx and hy are the width and height of the square.
	The region of integration assumed is the quadrant 1 unit square:
	x in [0,1], and y in [0,1]. The user should
	perform the needed changes of variables to make this so. The integrand
	must also be finite in the entire region.
	coarsedivs is a list of integer pairs specifying the number of divisions
	in the x direction and the number in y (all must therefore be >= 1).
	addlevels is a list of integers (all must be >=0), one for each pair in
	coarsedivs, and each specifies how many additional factors of 2 to apply
	to the divisions in coarsedivs."""

	if not len(coarsedivs) >= 1:
		raise IndexError( "Romberg2dO1: coarsedivs must be a sequence of 1 or more elements." )

	coarsedivs = list( coarsedivs )
	for i, d in enumerate(coarsedivs):
		if len(d) != 2:
			raise IndexError( "Romberg2dO1: elements of coarsedivs must have length 2." )

		coarsedivs[i] = tuple([ int(j) for j in d ])
		d = coarsedivs[i]

		if d[0] <= 0 or d[1] <= 0:
			raise ValueError( "Romberg2dO1: all integers in coarsedivs must be > 0." )
	coarsedivs = tuple( coarsedivs )

	if len(addlevels) != len(coarsedivs):
		raise IndexError( "Romberg2dO1: addlevels must have the same length as coarsedivs." )

	addlevels = tuple([ int(i) for i in addlevels ])

	for l in addlevels:
		if l < 0:
			raise ValueError( "Romberg2dO1: all elements of addlevels must be >= 0." )

	#Types and bounds out of the way, begin calculations
	estimates = []
	unitareas = []

	aveCorners = 0.25 * ( fun( 0.0, 0.0 ) + fun( 1.0, 0.0 ) +
						  fun( 0.0, 1.0 ) + fun( 1.0, 1.0 ) )

	for bdivs, adlvls in zip( coarsedivs, addlevels ):
		#Handle the base division pattern
		xdivs, ydivs = bdivs
		uarea = 1.0 / ( xdivs * ydivs )
		unitareas.append( uarea )

		curest = aveCorners * uarea

		#Add x-like sides, using a function closure to reduce lookup times for symbols in the reduce
		def MkAcc( function, weight, divs ):
			def Acc( tot, xi ):
				f = function
				w = weight
				d = float(divs)
				x = float(xi)/divs
				return tot + w * (f(x, 0.0) + f(x, 1.0))
			return Acc

		curest += reduce( MkAcc( fun, 0.5 * uarea, xdivs ) , xrange( 1, xdivs ), 0.0 )

		#Add y-like sides, using a function closure to reduce lookup times for symbols in the reduce
		def MkAcc( function, weight, divs ):
			def Acc( tot, yi ):
				f = function
				w = weight
				d = float(divs)
				y = float(yi)/divs
				return tot + w * (f(0.0, y) + f(1.0, y))
			return Acc

		curest += reduce( MkAcc( fun, 0.5 * uarea, ydivs ), xrange( 1, ydivs ), 0.0 )

		#Add center points
		def MkAcc( function, weight, xval, ydivs ):
			def Acc( tot, yi ):
				f = function
				w = weight
				d = float(ydivs)
				x = xval
				return tot + w * f(x, float(yi)/d)
			return Acc

		for xi in xrange( 1, xdivs ):
			x = float( xi ) / xdivs

			curest += reduce( MkAcc( fun, uarea, x, ydivs ),
							  xrange( 1, ydivs ), 0.0 )

		estimates.append( curest )

		#Handle additional levels
		for levnum in xrange( 1, adlvls + 1 ):
			xdivs <<= 1
			ydivs <<= 1
			uarea *= 0.25
			unitareas.append( uarea )

			curest = 0.0

			#Add x-like sides
			def MkAcc( function, weight, divs ):
				def Acc( tot, xi ):
					f = function
					w = weight
					d = float(divs)
					x = float(xi)/divs
					return tot + w * (f(x, 0.0) + f(x, 1.0))
				return Acc

			curest += reduce( MkAcc( fun, 0.5 * uarea, xdivs ), xrange( 1, xdivs, 2 ), 0.0 )

			#Add y-like sides
			def MkAcc( function, weight, divs ):
				def Acc( tot, yi ):
					f = function
					w = weight
					d = float(divs)
					y = float(yi)/divs
					return tot + w * (f(0.0, y) + f(1.0, y))
				return Acc

			curest += reduce( MkAcc( fun, 0.5 * uarea, ydivs ), xrange( 1, ydivs, 2 ), 0.0 )

			#Add center points on xi odd lines
			def MkAcc( function, weight, xval, ydivs ):
				def Acc( tot, yi ):
					f = function
					w = weight
					d = float(ydivs)
					x = xval
					return tot + w * f(x, float(yi)/d)
				return Acc

			for xi in xrange( 1, xdivs, 2 ):
				x = float( xi ) / xdivs

				curest += reduce( MkAcc( fun, uarea, x, ydivs ),
								  xrange( 1, ydivs ), 0.0 )

			#Add center points on xi even lines
			for xi in xrange( 2, xdivs, 2 ):
				x = float( xi ) / xdivs

				curest += reduce( MkAcc( fun, uarea, x, ydivs ),
								  xrange( 1, ydivs, 2 ), 0.0 )

			estimates.append( 0.25 * estimates[-1] + curest )

	sortinglist = zip( unitareas, estimates )
	sortinglist.sort( key= lambda x: x[0] )
	unitareas, estimates = zip( *sortinglist )

	return NevilleEvenRound( unitareas, estimates, 0.0 )


def Romberg1dO1( fun, coarsedivs=( 4, 5 ), addlevels=( 3, 3 ) ):
	"""Performs 1 dimensional Romberg integration based on an order 1
	estimator. The order one estimator is to divide the region of integration
	into cells and each cell contributes 0.5 * sum(fun(edges)) * hx
	to the integral, where hx is the width of the cell.
	The region of integration assumed is the unit interval:
	x in [0,1]. The user should
	perform the needed changes of variables to make this so. The integrand
	must also be finite in the entire region.
	coarsedivs is a list of integers specifying numbers of divisions
	in the x direction (all must therefore be >= 1).
	addlevels is a list of integers (all must be >=0), one for each element in
	coarsedivs, and each specifies how many additional factors of 2 to apply
	to the divisions in coarsedivs."""
	if not len(coarsedivs) >= 1:
		raise IndexError( "Romberg1dO1: coarsedivs must be a sequence of 1 or more elements." )

	coarsedivs = list( coarsedivs )
	for i, d in enumerate(coarsedivs):
		if d <= 0:
			raise ValueError( "Romberg1dO1: all integers in coarsedivs must be > 0." )
	coarsedivs = tuple( coarsedivs )

	if len(addlevels) != len(coarsedivs):
		raise IndexError( "Romberg1dO1: addlevels must have the same length as coarsedivs." )

	addlevels = tuple([ int(i) for i in addlevels ])

	for l in addlevels:
		if l < 0:
			raise ValueError( "Romberg1dO1: all elements of addlevels must be >= 0." )

	#Types and bounds out of the way, begin calculations
	estimates = []
	unitlens = []

	aveEdges = 0.5 * ( fun(0.0) + fun(1.0) )

	for bdiv, adlvls in zip( coarsedivs, addlevels ):
		#Handle the base division pattern
		xdivs = bdiv
		ulen = 1.0 / (xdivs)
		unitlens.append( ulen )

		curest = aveEdges * ulen

		#Add center points, using a function closure to reduce lookup times for symbols in the reduce
		def MkAcc( function, weight, divs ):
			def Acc( tot, xi ):
				f = function
				w = weight
				d = float(divs)
				return tot + w * f( float(xi)/divs )
			return Acc

		curest += reduce( MkAcc( fun, ulen, xdivs ), xrange( 1, xdivs ), 0.0 )
		estimates.append( curest )

		#Handle additional levels
		for levnum in xrange( 1, adlvls + 1 ):
			xdivs <<= 1
			ulen *= 0.5
			unitlens.append( ulen )

			#Add center points, using a function closure to reduce lookup times for symbols in the reduce
			def MkAcc( function, weight, divs ):
				def Acc( tot, xi ):
					f = function
					w = weight
					d = float(divs)
					return tot + w * f( float(xi)/divs )
				return Acc

			curest = reduce( MkAcc( fun, ulen, xdivs ), xrange( 1, xdivs, 2 ), 0.0 )

			estimates.append( 0.5 * estimates[-1] + curest )

	sortinglist = zip( unitlens, estimates )
	sortinglist.sort( key= lambda x: x[0] )
	unitlens, estimates = zip( *sortinglist )

	return NevilleEvenRound( unitlens, estimates, 0.0 )
