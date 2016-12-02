# this defines an amino acid class

# This file is for reading a pdb file and extract the atom coordinate information.
# The neccesary information is then writen to a text file witht the same name as
# the input file with different location.


# can improve later
import pickle
import os
import sys
import csv
import time 
from numpy.linalg import norm
import numpy as np


#_______________________________________________________________________________________________________________________________________________________

def get_angle(v1, v2, v3):
	#print vec_2, vec_3,vec_4
	normal_1 = np.cross(v1, v2)
	normal_2 = np.cross(v2, v3)
	angle = np.arctan2(np.dot(np.cross(normal_1,normal_2),v2)/norm(v2),np.dot(normal_1,normal_2))
	angle =  180*angle/3.14
	return round(angle,2)
	
#_______________________________________________________________________________________________________________________________________________________

def get_dihedral_angles(protein_sequence):
	#this method takes atomic coordinates of the protein
	#sequence and returns the dihedral angles of each amino acid in the sequence.
	#If the dihedral angles can\'t be computed then default value is set to 1000.
	#The computed angles are also written to a file named phi_psi_angles.csv. "
	
	C1 = []
	C2 = []
	CA = []
	N1 = []
	N2 = []
	angles = []
	angle_sum = []
	#print protein_sequence
	#print "length" , len(protein_sequence)
	atom1 = protein_sequence[0]
	atom2 = protein_sequence[1]
	atom3 = protein_sequence[2]
	atom4 = protein_sequence[3]
	N1 = []
	N1.append(atom1[-3])
	N1.append(atom1[-2])
	N1.append(atom1[-1])
	
	
	CA = []
	CA.append(atom2[-3])
	CA.append(atom2[-2])
	CA.append(atom2[-1])
	amino_name = atom2[2]
	#if amino_name != "CA":
		#return
	
	C2 = []
	C2.append(atom3[-3])
	C2.append(atom3[-2])
	C2.append(atom3[-1])
	
	N2 = []
	N2.append(atom4[-3])
	N2.append(atom4[-2])
	N2.append(atom4[-1])
	#print protein_sequence[6]
	
	#Get the vector for the base plane
	vec_2 =list((float(CA[0]-N1[0]),float(CA[1]-N1[1]),float(CA[2]-N1[2])))
	vec_3 =list((float(C2[0]-CA[0]),float(C2[1]-CA[1]),float(C2[2]-CA[2])))
	vec_4 =list((float(N2[0]-C2[0]),float(N2[1]-C2[1]),float(N2[2]-C2[2])))
	#normalize vectors
	#print vec_2
	
	
	angle_psi = get_angle(vec_2, vec_3,vec_4)
	
	angles.append((amino_name, 1000, angle_psi))
	#for x in protein_sequence:
		#print x
	
	#print CA1, CA2, CA3, CA4
	for i in range(4, len(protein_sequence)-2,3):	
		
		C1 = C2
		N1 = N2
		vec_1 = vec_4
		atom2 = protein_sequence[i]
		atom3 = protein_sequence[i+1]
		atom4 = protein_sequence[i+2]
		#print "in loop", protein_sequence[i]
		#print "in loop", protein_sequence[i+1]
		#print "in loop", protein_sequence[i+2]
		#print "jsfjhsgui next"
	
		
		CA = []
		CA.append(atom2[-3])
		CA.append(atom2[-2])
		CA.append(atom2[-1])
		amino_name = atom2[2]
		
		C2 = []
		C2.append(atom3[-3])
		C2.append(atom3[-2])
		C2.append(atom3[-1])
	
		N2 = []
		N2.append(atom4[-3])
		N2.append(atom4[-2])
		N2.append(atom4[-1])
		
		#Get the vector for the base plane
		
		vec_2 =list((float(CA[0]-N1[0]),float(CA[1]-N1[1]),float(CA[2]-N1[2])))/norm(vec_2)
		vec_3 =list((float(C2[0]-CA[0]),float(C2[1]-CA[1]),float(C2[2]-CA[2])))/norm(vec_3)
		vec_4 =list((float(N2[0]-C2[0]),float(N2[1]-C2[1]),float(N2[2]-C2[2])))/norm(vec_4)
		#print norm(vec_1),norm(vec_2),norm(vec_3),norm(vec_4)
		# NORMALS
		normal_vector1 = np.cross(vec_1,vec_2)
		normal_base = np.cross(vec_2, vec_3)
		normal_vector2 = np.cross(vec_3, vec_4)
		
		angle_phi = get_angle(vec_1, vec_2, vec_3)
		angle_psi = get_angle(vec_2,vec_3,vec_4)
		#print angle_phi,angle_psi, angle_phi + angle_psi
		angles.append((amino_name,angle_phi,angle_psi))
	
	# read last two atoms to get phi angle
	atom2 = protein_sequence[-2]
	atom3 = protein_sequence[-1]
	#print "sairam", protein_sequence[-2]
	#print "sairam", protein_sequence[-1]
	C1 = C2
	N1 = N2
	vec_1 = vec_4
	
	CA = []
	CA.append(atom2[-3])
	CA.append(atom2[-2])
	CA.append(atom2[-1])
	amino_name = atom2[2]
	#print atom2
	
	C2 = []
	C2.append(atom3[-3])
	C2.append(atom3[-2])
	C2.append(atom3[-1])
	
	vec_2 =list((float(CA[0]-N1[0]),float(CA[1]-N1[1]),float(CA[2]-N1[2])))
	vec_3 =list((float(C2[0]-CA[0]),float(C2[1]-CA[1]),float(C2[2]-CA[2])))
	#print "sairam",CA, N1,vec_2,vec_3
	if( 0.0 != norm(np.cross(vec_1,vec_2)) and 0.0 !=norm(np.cross(vec_2, vec_3))):
		angle_phi = get_angle(vec_1,vec_2,vec_3)
		#print angle_phi
	else:
		angle_phi = 1000
	angles.append((amino_name,angle_phi, 1000))
	angle_sum.append((0,0,1000))
	angle_sum.append((angles[1][1]+ angles[0][2],angles[1][2],angles[0][1]+ angles[0][2]))
	for i in range(2,len(angles)-1):
		angle_sum.append((round(angles[i][1]+angles[i-1][2],2),round(angles[i][2]+ angles[i-1][1],2),round(angles[i][1]+angles[i][2],2)))
		
	angle_sum.append((round(angles[-1][1]+angles[-2][2],2),angles[-2][1],1000))
	#for y in angles:
		#print y
	#for x in angle_sum:
		#print x
	ptr = open("phi_psi_angles.csv",'wb')
	try:
		writer = csv.writer(ptr)
		writer.writerow(('amino name','psi','phi'))
		for x in angles:
			writer.writerow(x)
	finally:
		ptr.close()
	
	
	return angles

#------------------------------------------------------------------------------------------#


		

def get_sequence(filename):
	# This function takes a pdb filename and 
	#returns the atomic coordinates for the protein
	#sequence and its the length of the sequence."
	sequence = [] 
	names = ['CA','N','C']
	x_coordinate, y_coordinate,z_coordinate = 0,0,0
	starting_atom  = True
	
	ptr = open(filename,"r")
	lines = ptr.readlines()
	for line in lines:
		temp_string = line.split()
		if temp_string[0] == "ATOM":
			atom_name = str(line[12:16]).strip()
			
			if  True == starting_atom and atom_name != "N":
				continue
			else:
				starting_atom = False
			if atom_name in names and (line[16] == " " or line[16] == "A"):
				acid_name = line[17:20].strip()
				chain_name = line[21].strip()
				x_coordinate = float(line[30:38].strip())
				y_coordinate = float(line[38:46].strip())
				z_coordinate = float(line[46:54].strip())
				t = (chain_name,atom_name,acid_name,x_coordinate,y_coordinate,z_coordinate)
				sequence.append(t) 
		if temp_string[0] == "ENDMDL" or temp_string[0]== "TER": #if more than one model are present then consider only one of them
			break
	
	return  len(sequence)/3, sequence



def display_dihedral_angles(angles):
	for an_angle in angles:
		print an_angle


# main execution starts here. 
# the program works on all the files that are present in the current directory.
def main():
	protein_sequence = []
	min_length = 10
	#start_path = '../new_data_1' #current directory
	start_path = '.' #current directory
	count  = 0
	for path,dirs,files in os.walk(start_path):
		for filename in files: #for each file check if it is a pdb file.
			if filename[-4:] == '.pdb' or filename[-4:]=='.ent': #if it is a protein file
				os.path.join(path,filename)
				absolute_path = path + "/"+ filename
				print filename
				seq_length, protein_sequence = get_sequence(absolute_path)
				#for y in protein_sequence:
					#print y, seq_length
				
				# if the sequence leagnth is not long enough 
				# then either it is not the protein structure 
				# or there is some wrong input.
				# so we assume the lenght of the protein to be atleast 10 
				if seq_length >= min_length:
					count += 1
					angle_list = get_dihedral_angles(protein_sequence)
					display_dihedral_angles(angle_list)
				else:
					print "The sequence length is too short! \n"
	if count == 0 :
		print "No protein file available! "
	
if __name__ == '__main__': main()	