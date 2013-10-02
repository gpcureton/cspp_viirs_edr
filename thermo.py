#!/usr/bin/env python
# encoding: utf-8
"""
thermo.py

Various thermodynamic relationship for ice and water.

Created by Geoff Cureton <geoff.cureton@ssec.wisc.edu> on 2009-04-04.
Copyright (c) 2009-2013 University of Wisconsin Regents. All rights reserved.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

file_Date = '$Date$'
file_Revision = '$Revision$'
file_Author = '$Author$'
file_HeadURL = '$HeadURL$'
file_Id = '$Id$'

__author__ = 'G.P. Cureton <geoff.cureton@ssec.wisc.edu>'
__version__ = '$Id$'
__docformat__ = 'Epytext'

from scipy import log10

def rh_to_mr( rh, p, t) :
	'''
	Returns mixing ratio, in g/kg, given relative humidity in %, 
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmix(p, t)

def rh_to_mr_wat( rh, p, t) :
	'''
	Returns mixing ratio over water, in g/kg, given relative humidity in %, 
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmixwat(p, t)

def rh_to_mr_ice( rh, p, t) :
	'''
	Returns mixing ratio over ice, in g/kg, given relative humidity in %, 
	pressure in hPa and temperature in K.
	'''
	return rh * 0.01 * satmixice(p, t)

def mr_to_rh( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio in g/kg,  
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmix(p, t)

def mr_to_rh_wat( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio over water in g/kg,  
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmixwat(p, t)

def mr_to_rh_ice( mr,  p,  t) :
	'''
	Returns relative humidity in %, given the mixing ratio over ice in g/kg,  
	pressure in hPa and temperature in K.
	'''
	return mr * 100. / satmixice(p, t)

def satmix( p, t) :
	'''
	Returns saturation mixing ratio in g/kg, given pressure in hPa and
	temperature in K.
	'''
	if (t > 253.) :
		return satmixwat(p, t)
	else :
		return satmixice(p, t)

def satmixwat( p,  t) :
	'''
	Returns saturation mixing ratio over water, in g/kg, given pressure in hPa and
	temperature in K.
	'''
	es = svpwat(t)
	return (622. * es)/p

def satmixice( p, t) :
	'''
	Returns saturation mixing ratio over ice, in g/kg, given pressure in hPa and
	temperature in K.
	'''
	es = svpice(t);
	return (622. * es) / p;


def svpwat(t) :
	'''
	Returns saturation vapor pressure over water, in hPa, given temperature in K.
	
	'''

	a0 =  0.999996876e0
	a1 = -0.9082695004e-2
	a2 =  0.7873616869e-4
	a3 = -0.6111795727e-6
	a4 =  0.4388418740e-8
	a5 = -0.2988388486e-10
	a6 =  0.2187442495e-12
	a7 = -0.1789232111e-14
	a8 =  0.1111201803e-16
	a9 = -0.3099457145e-19
	b = 0.61078e+1
	t -= 273.16
	return (b / ((a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+t*(a6+t*(a7+t*(a8+t*a9)))))))))**8.))

def svpice( t) :
	'''
	Returns saturation vapor pressure over ice, in hPa, given temperature in K.
	The Goff-Gratch equation (Smithsonian Met. Tables,  5th ed., pp. 350, 1984)
	'''
	a = 273.16 / t
	exponent = -9.09718 * (a - 1.) - 3.56654 * log10(a) + 0.876793 * (1. - 1./a) + log10(6.1071)

	return 10.0**exponent
