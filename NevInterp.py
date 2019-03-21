# NevInterp.py
#
# by Sean E. Lake
# Originally written at UCLA some time between 2010 and 2017 as part of my
# thesis project.
#
# This small library contains implementations of Neville's algorithm that are 
# different from the usual ones in textbooks. First, no assumption is made about
# the spacing of the x-values, so execution speed will be suffer in those special
# cases, but the applicability of the functions is broader. Second, when the 
# function being interpolated is known to be even, you can apply this library's
# 'Even' functions to get improved numerical performance.
#
# The primary application this library was designed for is performing Romberg
# integration with optimal numerical accuracy and as few decimations of the 
# interval as possible. It is the intrinsic evenness of numerical integration
# schemes, when properly implemented, that inspired the Even functions. 
# Of primary use is the 'Round' functions, for their robust performance in the
# face of divergence of the algorithm.

import math


def NevilleEven( inxvals, inyvals, xeval ):
	"""Implementation of Neville's algorithm for evaluating an nth order polynomial
	that goes through n+1 points (n+1 is the length of inpoints) at the point x.
	This implementation is modified by the assumption that the polynomial is purely
	even, so the differences in Neville's algorithm are implemented as differences of
	squares, ie (a+b)*(a-b) to preserve numerical accuracy.
	Also uses the last iterations to produce an error estimate. Returns
	the difference between the final answer and the closer of the two penultimate
	estimates as an error estimate. This is optimistic in the best cases, but there
	is no single factor that corrects the estimate.
	See:
	Weisstein, Eric W. "Neville's Algorithm." From MathWorld--A Wolfram Web Resource.
	http://mathworld.wolfram.com/NevillesAlgorithm.html"""
	av = sum(inyvals) / float(len(inyvals))
	var = reduce( lambda x, y : x + ( y - av )**2, inyvals, 0.0 )
	var /= float(len(inyvals))

	#print(av, var**.5)
	if len(inyvals) > 2:
		n = len( inxvals )
		l = len( inyvals )
		m = n - l

		newyvals = [ ( (xeval - x2)*(xeval + x2)*y1 - (xeval - x1)*(xeval + x1)*y2 )
					 / ( (x1 + x2)*(x1 - x2) )
					 for x1, x2, y1, y2 in zip( inxvals[:l-1], inxvals[m+1:m+l],
					 							inyvals[:-1], inyvals[1:] ) ]

		return NevilleEven( inxvals, newyvals, xeval=xeval )

	elif len(inyvals) == 2:
		x1 = inxvals[0]
		x2 = inxvals[-1]

		yret = (( (xeval - x2)*(xeval + x2)*inyvals[0] -
				  (xeval - x1)*(xeval + x1)*inyvals[1] ) /
				( (x1 + x2)*(x1 - x2) ) )
		yerr = max( abs(inyvals[0] - yret), abs(inyvals[1] - yret) )

		return( yret, yerr )

	else:
		return ( inyvals[0], None )

def Neville( inxvals, inyvals, xeval ):
	"""Implementation of Neville's algorithm for evaluating an nth order polynomial
	that goes through n+1 points (n+1 is the length of inpoints) at the point x.
	Also uses the last iterations to produce an error estimate. Returns
	the difference between the final answer and the closer of the two penultimate
	estimates as an error estimate. This is optimistic in the best cases, but there
	is no single factor that corrects the estimate.
	See:
	Weisstein, Eric W. "Neville's Algorithm." From MathWorld--A Wolfram Web Resource.
	http://mathworld.wolfram.com/NevillesAlgorithm.html"""
	if len(inyvals) > 2:
		n = len( inxvals )
		l = len( inyvals )
		m = n - l

		newyvals = [ ( (xeval - x2)*y1 - (xeval - x1)*y2 ) / (x1 - x2)
					 for x1, x2, y1, y2 in zip( inxvals[:l-1], inxvals[m+1:m+l],
					 							inyvals[:-1], inyvals[1:] ) ]

		return Neville( inxvals, newyvals, xeval=xeval )

	elif len(inyvals) == 2:
		x1 = inxvals[0]
		x2 = inxvals[-1]

		yret = (( (xeval - x2)*inyvals[0] - (xeval - x1)*inyvals[1] ) /
				(x1 - x2) )
		yerr = max( abs(inyvals[0] - yret), abs(inyvals[1] - yret) )

		return( yret, yerr )

	else:
		return ( inyvals[0], None )

def NevilleEvenRound( inxvals, inyvals, xeval ):
	"""Implementation of Neville's algorithm for evaluating an nth order polynomial
	that goes through n+1 points (n+1 is the length of inpoints) at the point x.
	This implementation is modified first by the assumption that the polynomial is purely
	even, so the differences in Neville's algorithm are implemented as differences of
	squares, ie (a+b)*(a-b) to preserve numerical accuracy. Also monitors the estimates
	for the effect of roundoff or overfitting (ie the variance stops decreasing) and
	reports out the average of the set with the minimum variances and the variance
	as an error estimate or uses the last iterations to produce an error estimate.
	Returns the difference between the final answer and the closer of the two penultimate
	estimates as an error estimate. This is optimistic in the best cases, but there
	is no single factor that corrects the estimate.
	See:
	Weisstein, Eric W. "Neville's Algorithm." From MathWorld--A Wolfram Web Resource.
	http://mathworld.wolfram.com/NevillesAlgorithm.html"""
	xvals = [ x for x in inxvals ]
	lastyvals = [ y for y in inyvals ]
	lastave = sum(lastyvals) / float( len( lastyvals ) )
	lastvar = reduce( lambda x,y : x + ( y - lastave )**2, lastyvals, 0.0 )
	lastvar /= float( len(lastyvals) )

	n = len( xvals )

	#Sort on distance from xeval to put the optimal extrapolation at 0
	sortArr = zip( [ abs(x - xeval) for x in xvals ], xvals, lastyvals )
	sortArr.sort( key=lambda x: x[0] )
	b, xvals, lastyvals = zip( *sortArr )

	del( b, sortArr )

	while len(lastyvals) > 2:
		l = len( lastyvals )
		m = n - l

		newyvals = [ ( (xeval - x2)*(xeval + x2)*y1 - (xeval - x1)*(xeval + x1)*y2 )
					 / ( (x1 + x2)*(x1 - x2) )
					 for x1, x2, y1, y2 in zip( xvals[:l-1], xvals[m+1:m+l],
					 							lastyvals[:-1], lastyvals[1:] ) ]
		newave = sum( newyvals ) / float( len(newyvals) )
		newvar = reduce( lambda x,y : x + ( y - newave )**2, newyvals, 0.0 )
		newvar /= float( len(newyvals) )

		if newvar <= lastvar:
			lastvar = newvar
			lastave = newave
			lastyvals = newyvals

		else:
			break

	if len(lastyvals) == 2:
		x1 = xvals[0]
		x2 = xvals[-1]

		yret = (( (xeval - x2)*(xeval + x2)*lastyvals[0] -
				  (xeval - x1)*(xeval + x1)*lastyvals[1] ) /
				  ( (x1 + x2)*(x1 - x2) ) )
		yerr = max( abs(lastyvals[0] - yret), abs(lastyvals[1] - yret) )

	elif len(lastyvals) > 2:
		yret = lastyvals[0] #lastave #the closest extrapolant is better than lastave
		yerr = math.sqrt( lastvar )

	elif len(lastyvals) == 1:
		yret = lastyvals[0]
		yerr = float("nan")

	else:
		raise ValueError("Neville algorithm not defined on empty lists.")

	return( yret, yerr )


def NevilleRound( inxvals, inyvals, xeval ):
	"""Implementation of Neville's algorithm for evaluating an nth order polynomial
	that goes through n+1 points (n+1 is the length of inpoints) at the point x.
	This implementation is modified to monitor the estimates
	for the effect of roundoff or overfitting (ie the variance stops decreasing) and
	reports out the average of the set with the minimum variances and the variance
	as an error estimate or uses the last iterations to produce an error estimate.
	Returns the difference between the final answer and the closer of the two penultimate
	estimates as an error estimate. This is optimistic in the best cases, but there
	is no single factor that corrects the estimate.
	See:
	Weisstein, Eric W. "Neville's Algorithm." From MathWorld--A Wolfram Web Resource.
	http://mathworld.wolfram.com/NevillesAlgorithm.html"""
	xvals = [ x for x in inxvals ]
	lastyvals = [ y for y in inyvals ]
	lastave = sum(lastyvals) / float( len( lastyvals ) )
	lastvar = reduce( lambda x,y : x + ( y - lastave )**2, lastyvals, 0.0 )
	lastvar /= float( len(lastyvals) )

	n = len( xvals )

	#Sort on distance from xeval to put the optimal extrapolation at 0
	sortArr = zip( [ abs(x - xeval) for x in xvals ], xvals, lastyvals )
	sortArr.sort( key=lambda x: x[0] )
	b, xvals, lastyvals = zip( *sortArr )

	del( b, sortArr )

	while len(lastyvals) > 2:
		l = len( lastyvals )
		m = n - l

		newyvals = [ ( (xeval - x2)*y1 - (xeval - x1)*y2 ) / (x1 - x2)
					 for x1, x2, y1, y2 in zip( xvals[:l-1], xvals[m+1:m+l],
					 							lastyvals[:-1], lastyvals[1:] ) ]
		newave = sum( newyvals ) / float( len(newyvals) )
		newvar = reduce( lambda x,y : x + ( y - newave )**2, newyvals, 0.0 )
		newvar /= float( len(newyvals) )

		if newvar <= lastvar:
			lastvar = newvar
			lastave = newave
			lastyvals = newyvals

		else:
			break

	if len(lastyvals) == 2:
		x1 = xvals[0]
		x2 = xvals[-1]

		yret = ( (xeval - x2)*lastyvals[0] - (xeval - x1)*lastyvals[1] ) / (x1 - x2)
		yerr = max( abs(lastyvals[0] - yret), abs(lastyvals[1] - yret) )

	elif len(lastyvals) > 2:
		yret = lastyvals[0] #lastave #the closest extrapolant is better than lastave
		yerr = math.sqrt( lastvar )

	elif len(lastyvals) == 1:
		yret = lastyvals[0]
		yerr = float("nan")

	else:
		raise ValueError("Neville algorithm not defined on empty lists.")

	return( yret, yerr )