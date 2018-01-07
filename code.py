from openpyxl import load_workbook
from astropy.coordinates import AltAz,BaseRADecFrame
from astropy import units as u
from math import sin,cos,asin,acos,atan,degrees,radians,sqrt
import ephem
import numpy as np
import sys

latitude = 22.572645
longitude = 88.3639
vMin = 0.428988806046
theta = 23.5
minDist = {}
bodies = {}
allTags = []
SMALL_NUMBER = sys.float_info.epsilon

mars = ephem.Mars()

def norm(v):
	return sqrt(np.dot(v,v))

def param(r, v):

	h = np.cross(r, v)
	n = np.cross([0, 0, 1], h)

	mu = 66.74

	#eccentricity
	e = (1.0/mu)*((norm(v) ** 2 - mu / norm(r))*np.array(r) - np.dot(r, v) * np.array(v))
	E = np.linalg.norm(v) ** 2 / 2 - mu / np.linalg.norm(r)

	#semi major axis
	a = mu / (2 * E)

	p = np.dot(h, h) / mu

	e1 = sqrt(np.dot(e, e))

	#inclination
	i = np.arccos(h[2]/ np.linalg.norm(h))

	#right ascension of ascending node angle
	raan = np.arccos(n[0]/np.linalg.norm(n))

	#angle of periapsis
	arg_pe = np.arccos(np.dot(n , e)/(np.linalg.norm(n) * np.linalg.norm(e)))

	if abs(e1 - 0) < SMALL_NUMBER:
		if abs(i - 0) < SMALL_NUMBER:
			f = acos(r.x / np.linalg.norm(r))
			if v.x > 0: f = 2 * np.pi - f
		else:
			f = acos(dot(n, r) / (np.linalg.norm(n) * np.linalg.norm(r)))
			if np.dot(n, v) > 0: f = 2 * np.pi - f
	else:
		if e[2] < 0: arg_pe = 2 * np.pi - arg_pe
		f = acos(np.dot(e, r) / (np.linalg.norm(e) * np.linalg.norm(r)))
		if np.dot(r, v) < 0: f = 2 * np.pi - f

 	return a, e1, i, raan, arg_pe, f 

def getVal(oldRa,oldDec,oldRadius):
	X1 = oldRadius*sin(oldRa)*cos(oldDec)
	Y1 = oldRadius*sin(oldRa)*sin(oldDec)
	Z1 = oldRadius*cos(oldRa)
	Z0 = 149598000
	X2,Y2,Z2 = X1,Y1,Z1-Z0
	return X2,Y2,Z2

def modifier(x,date):
	m = ephem.Mars(date)
	x.compute(date)
	r1,theta1,phi1 = float(m.earth_distance),float(m.ra),float(m.dec)
	# print r1,theta1,phi1
	# print m
	r2,theta2,phi2 = float(x.radius),float(x.ra),float(x.dec)
	# print r2,"o"
	# print r1*r1+r2*r2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2))
	value = sqrt(r1*r1+r2*r2-2*r1*r2*(sin(theta1)*sin(theta2)*cos(phi1-phi2)+cos(theta1)*cos(theta2)))
	return value

wb = load_workbook('all_asteroid_data.xlsx')
sheet = wb.get_sheet_by_name('Sheet1')
sheet['B1'],sheet['C1'] = "ra(degrees)","dec(degrees)",
for i in range(2,402):
	allTags.append(sheet["A"+str(i)].value)
	minDist[sheet["A"+str(i)].value] = (100000000000.0,ephem.Date('2018/01/05'))	
	aznum,altnum = "B"+str(i),"C"+str(i)
	az,alt = sheet[aznum].value,sheet[altnum].value
	sinD = sin(radians(alt))*sin(radians(latitude))+cos(radians(alt))*cos(radians(latitude))*cos(radians(az))
	dec = degrees(asin(sinD))
	sinH = -(sin(radians(az))*cos(radians(alt)))/cos(radians(dec))
	H = degrees(asin(sinH))
	ra = longitude-H
	dec -= 23.5
	sheet[aznum],sheet[altnum] = ra,dec
	sheet['F'+str(i)] = sheet['F'+str(i)].value + vMin
	vx,vz = sheet['E'+str(i)].value,sheet['G'+str(i)].value
	sheet['E'+str(i)] = (vx*cos(radians(theta)))-(vz*sin(radians(theta)))
	sheet['G'+str(i)] = (vx*sin(radians(theta)))+(vz*cos(radians(theta)))
	someX,someY,someZ = getVal(ra,dec,sheet['F'+str(i)].value)
	vx,vy,vz = sheet['E'+str(i)].value,sheet['F'+str(i)].value,sheet['G'+str(i)].value
	x = ephem.EllipticalBody(ephem.Date('2018/01/05'))
	x._a, x._e, x._inc, x._Om, x._om, x._M = param([someX,someY,someZ], [vx,vy,vz])
	bodies[sheet["A"+str(i)].value] = x
wb.save(filename='newdata3.xlsx')

startDate = ephem.Date('2018/01/05')
endDate = ephem.Date('2023/01/05')
currentDate = startDate

for tag in allTags:
	while currentDate < endDate:
		value = modifier(bodies[tag],currentDate)
		# print value
		if value < minDist[tag][0]:
			minDist[tag] = (value,ephem.Date(currentDate))
		currentDate += 1
	print tag+","+minDist[tag][0]+","+minDist[tag][1]