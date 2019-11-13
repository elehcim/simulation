import os


def maps(orbit_sideon):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    return appendix + '_maps.fits'

def maps_hi(orbit_sideon):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    return appendix + '_HI_maps.fits'


def maps_allband(orbit_sideon):
    appendix = "" if not orbit_sideon else "_orbit_sideon"
    return appendix + '_maps_allbands.fits'
