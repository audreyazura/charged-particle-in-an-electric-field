#Particle Movement: calculate the movement of a charged particle in an electric field given as input
#    Copyright (C) 2020  Alban Lafuente-Sampietro
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <https://www.gnu.org/licenses/>.

from FuncUtil import *
from math import sqrt
from numpy.random import normal
from os import path, mkdir, remove
from re import compile, search
from multiprocessing import Pool, current_process
#~ from progressbar import ProgressBar

#number of processus we want to spawn for that calculation
nproc = 6

#electron
particle = "Electron"
EfieldName = "Etotal"
maxSteps = 100000
nvelocity = 10000
m = 9.10938188*10**(-31)*0.089	#kg, electron effective mass, from this article https://www.sciencedirect.com/science/article/pii/S0921510710001960
q = -1.60217733*10**(-19)		#C, elementary charge

#~# #hole
#~ particle = "Hole"
#~ EfieldName = "Einternal"
#~ maxSteps = 50000
#~ nvelocity = 1000
#~ m = 9.10938188*10**(-31)*0.693		#kg, hole effective mass
#~ q = 1.60217733*10**(-19)			#C, elementary charge

#for notch position variations
inputFolder = "".join(["Extract/Notch/", EfieldName, "/"])
fullOutputGenFolder = "".join(["Extract/Notch/", particle, "Movement/MovingCalc/"])
accelerationGenFolder = "".join(["Extract/Notch/", particle, "Movement/Acceleration/"])
notchPositions = ["0", "200", "400", "600", "800", "1000", "1200", "1300", "1400", "1500", "1600", "1700", "1750", "1800", "1850", "1900", "1950", "2000"]
xinits = [1000.0, 1250.0, 1500.0, 1550.0, 1600.0, 1650.0, 1700.0, 1750.0, 1800.0, 1850.0]
#~ Efields = ["0.0", "0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8"]
Efields = ["0.8"]


#physical constant and conversions
kB = 1380649*10**(-23)				#J/K
T = 300								#K
dt = 10**(-12)						#step every fs
vth = sqrt(kB*T/m)					#m/s

#definition of initial velocities: if a file with them exists, we take it (if it has the right size), otherwise, we generate and save them
vinits = []
if path.isfile("ElectronSpeeds.dat"):
	with open("ElectronSpeeds.dat", "r") as speedFile:
		for speed in speedFile:
			vinits.append(float(speed))
	if not(len(vinits) == nvelocity):
		remove("ElectronSpeeds.dat")
	
if not(path.isfile("ElectronSpeeds.dat")):
	vinits = [1.0, 1.0]
	while(not(len(vinits)==len(set(vinits)))):
		vinits = normal(0, vth, nvelocity)
	with open("ElectronSpeeds.dat", "w") as speedFile:
		for speed in vinits:
			print(speed, file=speedFile)

#function definition
def CalculateMovement (x0, v0, maxLoop, absisse, Efield):
	xt = x0
	vx = v0
	ax = 0.0
	
	xcalc = []
	acalc = []
	vcalc = []
	tcalc = []
	
	for nsteps in range(maxLoop):
		try:
			Ex = FindValue(xt/1000, absisse, Efield)
			ax = q*Ex/m
			acalc.append(ax)
			
			vchange = (vx+ax*dt)/2
			vx += ax*dt
			vcalc.append(vx)

			xcalc.append(xt)
			xt *= 10**(-6)	#converting x to meter
			xt += vchange*dt	#moving the position
			xt *= 10**6		#converting x back to micrometer
			
			tcalc.append(nsteps*dt*10**9)
			nsteps += 1
		except ValueError as err:
			if not(err.args[0] == "max" or err.args[0] == "min"):
				raise
			
			nsteps = maxLoop+1
			return xcalc, vcalc, acalc, tcalc, err.args[0]
	return xcalc, vcalc, acalc, tcalc, 0

def WorkerCalc (E, pos):
	accelerationFolder = "".join([accelerationGenFolder, "E", E, "V/"])
	if not(path.isdir(accelerationFolder)):
		mkdir(accelerationFolder)
	
	fullOutputEfieldFolder = "".join([fullOutputGenFolder, "E", E, "V/"])
	if not(path.isdir(fullOutputEfieldFolder)):
		mkdir(fullOutputEfieldFolder)

	inputFile = "".join([inputFolder, "Light/E", E, "V/E", E, "V_Notch", pos, "nm_Light.sim"])
	accelerationFile = "".join([accelerationFolder, "Notch", pos, "nm_Acceleration.sim"])
	
	notchFolder = "".join([fullOutputEfieldFolder, "Notch", pos, "nm/"])
	if not(path.isdir(notchFolder)):
		mkdir(notchFolder)
	
	(depth, Efield) = ExtractTables(inputFile, [1,3], float)
	indexFront = FindIndex(2.05, depth)
	EfieldNC = [E*100 for E in Efield]		#converting the electric from V/cm to N/C
	
	with open(accelerationFile, "w") as accFile:
		print("Depth (microm)\tAcceleration (m/s)", file=accFile)
		for i in range(len(depth)):
			acc = q*EfieldNC[i]/m
			print("\t".join([str(depth[i]), str(acc), str(q), str(EfieldNC[i])]), file=accFile)
	print(" ".join([current_process().name, ": File", accelerationFile, "created."]))
	
	EfieldtruncNC = EfieldNC[:indexFront]		#truncating Efield so the electron stays in the absorber and don't see the border of the sample
	depthTrunc = depth[:indexFront]				#truncating depth so it has the same size as Efield and wont cause problem when trying to find Efield values
	
	posFloat = float(pos)
	xcopy = xinits.copy()
	if not(posFloat in xinits):
		xcopy.append(posFloat)
	
	for xi in xcopy:
		fullOutputFolder = "".join([notchFolder, "xi", str(int(xi)), "nm/"])
		if not(path.isdir(fullOutputFolder)):
			mkdir(fullOutputFolder)
		
		exitFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_Exit.sim"])
		meanFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_Mean.sim"])
		frontFastFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_FrontFast.sim"])
		frontSlowFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_FrontSlow.sim"])
		frontMeanFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_FrontMean.sim"])
		backFastFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_BackFast.sim"])
		backSlowFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_BackSlow.sim"])
		backMeanFile = "".join([fullOutputFolder, "Notch", pos, "nm_xi=", str(int(xi)), "nm_BackMean.sim"])
		
		if not(path.isfile(backSlowFile)):
			xMean = []
			vMean = []
			aMean = []
			tMean = []
			
			xFrontFast = []
			vFrontFast = []
			aFrontFast = []
			tFrontFast = []
			
			xFrontSlow = []
			vFrontSlow = []
			aFrontSlow = []
			tFrontSlow = []
			
			xFrontMean = []
			vFrontMean = []
			aFrontMean = []
			tFrontMean = []
			
			xBackFast = []
			vBackFast = []
			aBackFast = []
			tBackFast = []
			
			xBackSlow = []
			vBackSlow = []
			aBackSlow = []
			tBackSlow = []
			
			xBackMean = []
			vBackMean = []
			aBackMean = []
			tBackMean = []
			
			frontCounter = 0
			backCounter = 0
			
			frontFirst = True
			backFirst = True
			
			print("".join([current_process().name, " : Beginning of velocity calculation for x_notch = ", pos, " nm and x_initial = ", str(xi), "nm."]))
			#~ pbar = ProgressBar()
			#~ for vi in pbar(vinits):
			for vi in vinits:
				
				(xtab, vtab, atab, ttab, lifterr) = CalculateMovement(xi, vi, maxSteps, depthTrunc, EfieldtruncNC)
				tSize = len(ttab)
				
				#Calculating the mean path for all particles
				tMeanSize = len(tMean)
				for i in range(max(tSize, tMeanSize)):
					if i < tMeanSize:
						#len(ttab) can be superior or inferior to len(tsMean), so we have to test both cases
						if i < tSize:
							xMean[i] += xtab[i]
							vMean[i] += vtab[i]
							aMean[i] += atab[i]
						if i >= tSize:
							xMean[i] += xtab[tSize-1]
							vMean[i] += vtab[tSize-1]
							aMean[i] += atab[tSize-1]
					if i >= tMeanSize:
						#here, we have to be in the case len(ttab) > len(tMean), therefore i will go to len(ttab)-1 and cannot be above it, so we don't have to test
						tMean.append(ttab[i])
						if(tMeanSize == 0):
							xMean.append(xtab[i])
							vMean.append(vtab[i])
							aMean.append(atab[i])
						else:
							xMean.append(xMean[tMeanSize-1] + xtab[i])
							vMean.append(vMean[tMeanSize-1] + vtab[i])
							aMean.append(aMean[tMeanSize-1] + atab[i])
				
				if lifterr == "max":
					#Increasing the number of electron having went out at the front
					frontCounter += 1
					frontFirst = (len(tFrontFast) == 0)
					
					#if nothing has been put in tFrontFast, we put the first that arrive. If not, we test if this electron have been faster, and if so, we change teh FrontFast table
					if frontFirst:
						xFrontFast = xtab
						vFrontFast = vtab
						aFrontFast = atab
						tFrontFast = ttab
					if tSize < len(tFrontFast):
						xFrontFast = xtab
						vFrontFast = vtab
						aFrontFast = atab
						tFrontFast = ttab
						
					#No need to initialize xFrontSlow: if the list is empty, the first table to come will always be longer than it
					if tSize > len(tFrontSlow):
						xFrontSlow = xtab
						vFrontSlow = vtab
						aFrontSlow = atab
						tFrontSlow = ttab
					
					#Calculating the mean. If the new electron stayed longer that the longest electron in the mean, we have to make the table bigger.
					tMeanSize = len(tFrontMean)
					for i in range(max(tSize, tMeanSize)):
						if i < tMeanSize:
							#len(ttab) can be superior or inferior to len(tFrontMean), so we have to test both cases
							if i < tSize:
								xFrontMean[i] += xtab[i]
								vFrontMean[i] += vtab[i]
								aFrontMean[i] += atab[i]
							if i >= tSize:
								xFrontMean[i] += depth[indexFront]*1000					#depth[indexFront] is in micrometer, not nm
								vFrontMean[i] += vtab[tSize-1]
								aFrontMean[i] += atab[tSize-1]
						if i >= tMeanSize:
							#here, we have to be in the case len(ttab) > len(tFrontMean), therefore i will go to len(ttab)-1 and cannot be above it, so we don't have to test
							tFrontMean.append(ttab[i])
							if frontFirst:
								xFrontMean.append(xtab[i])
								vFrontMean.append(vtab[i])
								aFrontMean.append(atab[i])
							elif not(frontFirst):
								#print(str(len(tFrontMean))+"\t"+str(len(ttab))+"\t"+str(len(xtab)))
								xFrontMean.append(depth[indexFront]*1000*(frontCounter) + xtab[i])	#depth[indexFront] is in micrometer, not nm
								vFrontMean.append(vFrontMean[tMeanSize-1] + vtab[i])
								aFrontMean.append(aFrontMean[tMeanSize-1] + atab[i])		
					
				#We do the same thing with the back, see comment above
				if lifterr == "min":
					backCounter += 1
					backFirst = (len(tBackFast) == 0)
					if backFirst:
						xBackFast = xtab
						vBackFast = vtab
						aBackFast = atab
						tBackFast = ttab
					if tSize < len(tBackFast):
						xBackFast = xtab
						vBackFast = vtab
						aBackFast = atab
						tBackFast = ttab
					if tSize > len(tBackSlow):
						xBackSlow = xtab
						vBackSlow = vtab
						aBackSlow = atab
						tBackSlow = ttab
					tMeanSize = len(tBackMean)
					for i in range(max(tSize, tMeanSize)):
						if i < tMeanSize:
							if i < tSize:
								xBackMean[i] += xtab[i]
								vBackMean[i] += vtab[i]
								aBackMean[i] += atab[i]
							if i >= tSize:
								xBackMean[i] += 0
								vBackMean[i] += vtab[tSize-1]
								aBackMean[i] += atab[tSize-1]
						if i >= tMeanSize:
							tBackMean.append(ttab[i])
							if backFirst:
								xBackMean.append(xtab[i])
								vBackMean.append(vtab[i])
								aBackMean.append(atab[i])
							elif not(backFirst):
								xBackMean.append(xtab[i])
								vBackMean.append(vBackMean[i-1] + vtab[i])
								aBackMean.append(aBackMean[i-1] + atab[i])
			
			print("".join([current_process().name, " : Calculation ended for x_notch = ", pos, " nm and x_initial = ", str(xi), "nm."]))
			
			for i in range(len(tMean)):
				xMean[i] /= nvelocity
				vMean[i] /= nvelocity
				aMean[i] /= nvelocity
			
			for i in range(len(tFrontMean)):
				xFrontMean[i] /= frontCounter
				vFrontMean[i] /= frontCounter
				aFrontMean[i] /= frontCounter
			
			for i in range(len(tBackMean)):
				xBackMean[i] /= backCounter
				vBackMean[i] /= backCounter
				aBackMean[i] /= backCounter
			
			with open(exitFile, "w") as outFile:
				print("Absisse\tPosition\tNumber of carrier on this side", file=outFile)
				print("\t".join(["0", "Not accounted", str(nvelocity - frontCounter - backCounter)]), file=outFile)
				print("\t".join(["1", "Front", str(frontCounter)]), file=outFile)
				print("\t".join(["2", "Back", str(backCounter)]), file=outFile)
			print(" ".join([current_process().name, ": File", exitFile, "created."]))
			
			with open(meanFile, "w") as outMean:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outMean)
				for j in range(len(tMean)):
					print("\t".join([str(tMean[j]), str(xMean[j]), str(vMean[j]), str(aMean[j])]), file=outMean)
			print(" ".join([current_process().name, ": File", meanFile, "created."]))
			
			with open(frontFastFile, "w") as outFrontFast:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outFrontFast)
				for j in range(len(tFrontFast)):
					print("\t".join([str(tFrontFast[j]), str(xFrontFast[j]), str(vFrontFast[j]), str(aFrontFast[j])]), file=outFrontFast)
			print(" ".join([current_process().name, ": File", frontFastFile, "created."]))
			
			with open(frontSlowFile, "w") as outFrontSlow:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outFrontSlow)
				for j in range(len(tFrontSlow)):
					print("\t".join([str(tFrontSlow[j]), str(xFrontSlow[j]), str(vFrontSlow[j]), str(aFrontSlow[j])]), file=outFrontSlow)
			print(" ".join([current_process().name, ": File", frontSlowFile, "created."]))
			
			with open(frontMeanFile, "w") as outFrontMean:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outFrontMean)
				for j in range(len(tFrontMean)):
					print("\t".join([str(tFrontMean[j]), str(xFrontMean[j]), str(vFrontMean[j]), str(aFrontMean[j])]), file=outFrontMean)
			print(" ".join([current_process().name, ": File", frontMeanFile, "created."]))
				
			with open(backFastFile, "w") as outBackFast:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outBackFast)
				for j in range(len(tBackFast)):
					print("\t".join([str(tBackFast[j]), str(xBackFast[j]), str(vBackFast[j]), str(aBackFast[j])]), file=outBackFast)
			print(" ".join([current_process().name, ": File", backFastFile, "created."]))
			
			with open(backSlowFile, "w") as outBackSlow:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outBackSlow)
				for j in range(len(tBackSlow)):
					print("\t".join([str(tBackSlow[j]), str(xBackSlow[j]), str(vBackSlow[j]), str(aBackSlow[j])]), file=outBackSlow)
			print(" ".join([current_process().name, ": File", backSlowFile, "created."]))
			
			with open(backMeanFile, "w") as outBackMean:
				print("Time (ns)\tPosition (microm)\tSpeed (m/s)\tAcceleration (m^2/s)", file=outBackMean)
				for j in range(len(tBackMean)):
					print("\t".join([str(tBackMean[j]), str(xBackMean[j]), str(vBackMean[j]), str(aBackMean[j])]), file=outBackMean)
			print(" ".join([current_process().name, ": File", backMeanFile, "created."]))

if not(path.isdir("".join(["Extract/Notch/", particle, "Movement"]))):
	mkdir("".join(["Extract/Notch/", particle, "Movement"]))
if not(path.isdir(fullOutputGenFolder)):
	mkdir(fullOutputGenFolder)
if not(path.isdir(accelerationGenFolder)):
	mkdir(accelerationGenFolder)

print("")

#starting parallelization
workerpool = Pool(nproc)
valuelist = [(E, notch) for E in Efields for notch in notchPositions]
workerpool.starmap(WorkerCalc, valuelist)
