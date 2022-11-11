#!/usr/bin/env python3.8 -i

"""PanStarrsSkyCellTools: a few tools to comput the coordinates of PanStarrs PS1 skycells, and find the skycell that contains any given coordinate.

First Written: 2022-11-11
Author: Sean E. Lake
Institution: NAOC"""


import os, re
from bisect import bisect_right

import numpy as np
import astropy.io.fits as pf
import astropy.wcs as wcslib

Projcellsize = 4.0
minra, maxra = (0.0, 360.0)
mindec, maxdec = (-90.0, 90.0)
radperdeg = np.pi / 180.0
pixsize = 0.25 / 3600.0 # degrees per pixel
skycellspercell = 10 # width/height of a projection cell in sky cells
bufferpix = 240 # additional pixels added to each image edge to create overlap
skycellformat = re.compile("^\d{4}\.0\d{2}$")
filenameregex = re.compile("skycell\.\d{4}\.0\d{2}\.")


# The file ps1grid.fits comes from the PanStarrs team. It was downloaded from:
# https://outerspace.stsci.edu/download/attachments/10257181/ps1grid.fits?version=3&modificationDate=1532367528459&api=v2
__PS1SkycellCache = os.path.join(os.path.split(__file__)[0], "ps1grid.fits")

with pf.open(__PS1SkycellCache) as ff:
    __PS1SkycellData = ff[1].data

del(ff)

minzone = np.amin(__PS1SkycellData["ZONE"])
maxzone = np.amax(__PS1SkycellData["ZONE"])
maxcellnum = np.amax(__PS1SkycellData["PROJCELL"]) # the pole has one cell


def CoordinateToProjcell(ra, dec):
    """"""
    if not (minra <= ra <= maxra):
        raise ValueError(
            f"CoordinateToProjcell: ra {ra} out of bounds [{minra}, {maxra}].")
    
    if not (mindec <= dec <= maxdec):
        raise ValueError(
            f"CoordinateToProjcell: dec {dec} out of bounds [{mindec}, {maxdec}].")
    
    # first guess at the best sone id
    zone = int(np.round((dec - mindec) / Projcellsize))
    cacherow = zone - minzone
    if cacherow < 0:
        return None

    zonelo = __PS1SkycellData["DEC_MIN"][cacherow]
    zonehi = __PS1SkycellData["DEC_MAX"][cacherow]
    if dec < zonelo:
        zone -= 1
        cacherow -= 1
        if cacherow < 0:
            return None
    elif dec > zonehi:
        zone += 1
        cacherow += 1

    zonecells = __PS1SkycellData["NBAND"][cacherow]
    cellnum0 = __PS1SkycellData["PROJCELL"][cacherow]

    cellnum = int(np.round((ra - minra) / (maxra - minra) * zonecells))
    cellnum %= zonecells # catch the coordinates with ra = maxra - epsilon
    cellnum += cellnum0
    return cellnum


def ProjcellToCacheRow(cellnum):
    return bisect_right(__PS1SkycellData["PROJCELL"], cellnum) - 1


def ProjcellToCenter(cellnum):
    try:
        cellnum = int(cellnum)
    except ValueError:
        raise TypeError(
            "ProjcellToCenter: cellnum must be (convertable to) an int")

    if cellnum < 0 or cellnum > maxcellnum:
        raise ValueError(
            f"ProjcellToCenter: cellnum {cellnum} out of bounds [0,{maxcellnum}]")
    
    cacherow = ProjcellToCacheRow(cellnum)
    zonecells = __PS1SkycellData["NBAND"][cacherow]
    cellnum0 = __PS1SkycellData["PROJCELL"][cacherow]
    zonedec = __PS1SkycellData["DEC"][cacherow]
    zoneid = cellnum - cellnum0
    cellra = float(zoneid) / zonecells * (maxra - minra) + minra
    return (cellra, zonedec)


def SkycellDimToProjcellDim(size):
    return (size - 2*bufferpix) / skycellspercell + 2*bufferpix


def ProjcellDimToSkycellDim(size):
    return (size - 2*bufferpix) * skycellspercell + 2*bufferpix


def ProjcellWCS(Pcell):
    Pra, Pdec = ProjcellToCenter(Pcell)
    cacherow = ProjcellToCacheRow(Pcell)
    skycellwid, skycellhei = [
        __PS1SkycellData[s][cacherow] for s in ("XCELL", "YCELL") ]
    Pcellwid, Pcellhei = [
        ProjcellDimToSkycellDim(i) for i in (skycellwid, skycellhei) ]

    projhdr = pf.Header()
    projhdr["BITPIX"] = -32
    projhdr["NAXIS"] = 2
    projhdr["NAXIS1"] = Pcellwid
    projhdr["NAXIS2"] = Pcellhei
    projhdr["CDELT1"] = pixsize
    projhdr["CDELT2"] = pixsize
    projhdr["CTYPE1"] = "RA---TAN"
    projhdr["CTYPE2"] = "DEC--TAN"
    projhdr["CRVAL1"] = Pra
    projhdr["CRVAL2"] = Pdec
    projhdr["CRPIX1"] = 0.5 * Pcellwid + 1.0  # fits convention
    projhdr["CRPIX2"] = 0.5 * Pcellhei + 1.0  # fits convention
    
    projhdr["PC1_1"] = -1.0
    projhdr["PC1_2"] = 0.0
    projhdr["PC2_1"] = 0.0
    projhdr["PC2_2"] = 1.0

    Pcellwcs = wcslib.WCS(header=projhdr)

    return Pcellwcs


def CoordinateToSkycell(ra, dec):
    if not (minra <= ra <= maxra):
        raise ValueError(
            f"CoordinateToSkycell: ra {ra} out of bounds [{minra}, {maxra}].")
    
    if not (mindec <= dec <= maxdec):
        raise ValueError(
            f"CoordinateToSkycell: dec {dec} out of bounds [{mindec}, {maxdec}].")

    Pcell = CoordinateToProjcell(ra, dec)
    wcs = ProjcellWCS(Pcell)
    
    cacherow = ProjcellToCacheRow(Pcell)
    skycellwid, skycellhei = [
        __PS1SkycellData[s][cacherow] for s in ("XCELL", "YCELL") ]
    
    Px, Py = wcs.all_world2pix(np.array([[ra, dec]]), 1)[0]
    
    icell = int(np.floor((Px - bufferpix) / (skycellwid - 2*bufferpix)))
    jcell = int(np.floor((Py - bufferpix) / (skycellhei - 2*bufferpix)))

    return (f"{Pcell:04d}.0{jcell}{icell}")


def SkycellToCoordinate(skycell):
    if not isinstance(skycell, str):
        raise TypeError(
            "SkycellToCoordinate: skycell must be a string")

    if not skycellformat.match(skycell):
        raise ValueError(
            "SkycellToCoordinate: skycell must have the format 'dddd.0dd'")

    Pcell, skycellcode = skycell.split(".")
    Pcell = int(Pcell)
    
    jcell, icell = [ int(s) for s in skycellcode[1:] ]
    cacherow = ProjcellToCacheRow(Pcell)
    skycellwid, skycellhei = [
        __PS1SkycellData[s][cacherow] for s in ("XCELL", "YCELL") ]

    Xcell = (icell + 0.5) * (skycellwid - 2*bufferpix) + bufferpix
    Ycell = (jcell + 0.5) * (skycellhei - 2*bufferpix) + bufferpix

    wcs = ProjcellWCS(Pcell)
    ra, dec = wcs.all_pix2world(np.array([[Xcell, Ycell]]), 1)[0]
    return (ra, dec)


def ExtractSkycell(filename):
    if not isinstance(filename, str):
        raise TypeError("ExtractSkycell: filename must be a string.")
    
    sch = filenameregex.search(filename)
    if sch is None:
        raise ValueError(
            f"ExtractSkycell: filename '{filename}' was formatted incorrectly.")

    skcstart, skcend = sch.span()
    skcend -= 1 # drop trailing period
    skcstart += len("skycell.") # drop leading "skycell."
    return filename[skcstart:skcend]
    

if __name__ == "__main__":
    pass
