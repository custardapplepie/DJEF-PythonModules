# RADEC module

# dec2deg, ra2deg Written by Enno Middelberg 2001
# Made to module by David Floyd

import sys
import string

def dec2deg(s):
	print s
	inp = s.split(":")
	print inp
	if len(inp)<3:
		raise ValueError("Too few fields in H:M:S string.")
	elif len(inp) > 3:
		raise ValueError("Too many fields in H:M:S string.")
	hh=abs(float(inp[0]))
	mm=float(inp[1])/60.
	ss=float(inp[2])/3600.
	if float(inp[0]) < 0 or inp[0]=='-00' or inp[0]=='-0':
		print" That's -"+str(hh+mm+ss)+" degrees"
		dec = -(hh+mm+ss)
	else:
		print" That's "+str(hh+mm+ss)+" degrees"
		dec = hh+mm+ss
	return (dec)


def ra2deg(s):
	print s
	inp = s.split(":")
	if len(inp)<3:
		raise ValueError("Too few fields in H:M:S string.")
	elif len(inp) > 3:
		raise ValueError("Too many fields in H:M:S string.")
	hh=(float(inp[0]))*15.
       	mm=(float(inp[1])/60.)*15.
	ss=(float(inp[2])/3600.)*15.
	print" That's "+str(hh+mm+ss)+" degrees"
	ra = hh+mm+ss
	return (ra)
