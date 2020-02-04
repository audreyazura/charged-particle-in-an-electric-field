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

#Extract columns of data from a file
#Warning: only work for numerical data
#Warning: ignores no-utf8 characters. If you want them to raise an error, remove the errors="ignore" parameter in the file opening
def ExtractTables(fileName, column, dataType):
	from sys import exit
	
	with open(fileName, 'r', errors="ignore") as Input:
	
		extractedLists = []
		nColumn = len(column)
		printing = False
		cont = True
	
		for i in range(nColumn):
			extractedLists.append([])
	
		for line in (y for y in Input if cont):			
			cleanLine = line.strip()
			if printing:
				cont = not(cleanLine == "")
				printing = False
			
			splitLine = cleanLine.split("\t")
			selectLine = TestFunc(dataType, splitLine[0])
			
			if selectLine and cont:
				printing = True
				for i in range(nColumn):
					extractedLists[i].append(dataType(splitLine[column[i]]))
					
	return extractedLists

#return the index of a number in a list, or the first index containing a number above (if inf = False) or below (if inf = True) the number in value.
#Warning: will only find the first occurence of the number
def FindIndex(value, table, inf=False):
	index = 0
	
	if value < min(table):
		index = None
		raise ValueError("min")
	
	elif value > max(table):
		index = None
		raise ValueError("max")
	
	elif value > min(table) and value < max (table):
		find = False
		sortedList = table.copy()
		sortedList.sort()
		tableSize = len(table)
		i = 0
		if inf:
			while not(find) and i < tableSize:
				if sortedList[i] >= value:
					tableValue = sortedList[i]
					find = True
				elif sortedList[i] < value:
					i += 1
		#to find the highest value in table that is still lower than value, we find the first value in table that is above value, and take it
		elif not(inf):
			while not(find) and index < (tableSize-1):
				if sortedList[i+1] >= value:
					tableValue = sortedList[i]
					find = True
				elif sortedList[i+1] < value:
					i += 1
		
		if find:
			index = table.index(tableValue)
		if not(find):
			index = None
			raise ValueError("not found")
			
	return index

#return the value of a math function given as a table at an absisse find in an abisse table, either by finding the index of the absisse, or by doing a linear interpolation from the closest absisses in the table
def FindValue(absisse, absTable, funcTable):
	try:
		indAbsInf = FindIndex(absisse, absTable, False)
		indAbsSup = FindIndex(absisse, absTable, True)
	except ValueError:
		raise
	
	if indAbsInf == indAbsSup:
		value = funcTable[indAbsInf]
	
	elif not(indAbsInf == indAbsSup):
		a = (funcTable[indAbsSup] - funcTable[indAbsInf])/(absTable[indAbsSup] - absTable[indAbsInf])
		b = funcTable[indAbsSup] - a*absTable[indAbsSup]
		
		value = a*absisse + b
	
	return value
	
#Test if the user input correspond to something expected
def TestFunc(passFunc, data):
	try:
		passFunc(data)
		if(isinstance(passFunc(data), bool)):
			return passFunc(data)
	except ValueError:
		return False
	return True
