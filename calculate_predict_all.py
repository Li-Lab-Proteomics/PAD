import string
from math import *
import math
from numpy import *
import numpy as np
import gc
import sys

file1 = "f_disease_list.txt"
file2 = "f_protein_list.txt"
file3 = "f_drug_list.txt"

file4 = "f_gene_disease.txt"
file5 = "f_protein_protein.txt"
file6 = "f_drug_target.txt"
file7 = "f_drug_drug_similarity.txt"
file8 = "f_disease_disease_similarity.txt"
file11 = "f_drug_disease.txt"

file9 = "final_matrix_" + sys.argv[1] + "_" + sys.argv[2] + "_" + sys.argv[3] + "_" + sys.argv[4] + ".txt"
file10 = "process_" + sys.argv[1] + "_" + sys.argv[2] + "_" + sys.argv[3] + "_" + sys.argv[4] + ".txt"


infile1 = open(file9,"w")
infile2 = open(file10,"w")

#--------------------------------------------------

#Read in parameters
a = float(sys.argv[1])
b = float(sys.argv[2])
c = float(sys.argv[3])
r = float(sys.argv[4])

#Change between p1 and p2 measured by L1 norm
def L1_norm(p1,p2):
	p = p1 - p2
	p = p.tolist()
	sum = 0
	for i in range(len(p)):
		sum += abs(p[i][0])
	return sum

#Return all related targets of certain drug
def get_target(drug):
	drugname = drugs[i]
	target = []
	fromfile8 = open(file6,"r")
	for line in fromfile8:
		if line.split():
			temp = line[:-1].split("\t")
			if temp[0] == drugname:
				target.append(temp[1])
	return target

#disease list
diseases = []
fromfile1 = open(file1,"r")
for line in fromfile1:
	if line.split():
		diseases.append(line[:-1])
fromfile1.close()
disease_num = len(diseases)
infile2.write("Complete reading disease_list\n")
infile2.flush()

#protein list
proteins = []
fromfile2 = open(file2,"r")
for line in fromfile2:
	if line.split():
		proteins.append(line[:-1])
fromfile2.close()
protein_num = len(proteins)
infile2.write("Complete reading protein_list\n")
infile2.flush()

#drug list
drugs = []
fromfile3 = open(file3,"r")
for line in fromfile3:
	if line.split():
		drugs.append(line[:-1])
fromfile3.close()
drug_num = len(drugs)
infile2.write("Complete reading drug_list\n")
infile2.flush()

#s5,m7(disease -> disease)
fromfile9 = open(file8,"r")
s5 = [([0]*disease_num) for i in range(disease_num)]
m7 = [([0]*disease_num) for i in range(disease_num)]
for line in fromfile9:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = diseases.index(temp[0])
		row_num = diseases.index(temp[1])
		s5[col_num][row_num] = float(temp[2])
		s5[row_num][col_num] = float(temp[2])

for i in range(disease_num):
	colsum = sum(s5[i])
	for j in range(disease_num):
		m7[i][j] = a * s5[i][j]/colsum

m7 = np.matrix(m7)
infile2.write("Complete constructing m7\n")
infile2.flush()
del s5

#s1,m1(disease->protein)
fromfile4 = open(file4,"r")
s1 = [([0]*protein_num) for i in range(disease_num)]
m1 = [([0]*protein_num) for i in range(disease_num)]
for line in fromfile4:
	if line.split():
		temp = line[:-1].split("\t")
		row_num = diseases.index(temp[1])
		col_num = proteins.index(temp[0])
		s1[row_num][col_num] = 1

for i in range(disease_num):
	colsum = sum(s1[i])
	for j in range(protein_num):
		m1[i][j] = b * s1[i][j]/colsum

m1 = np.matrix(m1)
fromfile4.close()
infile2.write("Complete constructing m1\n")
infile2.flush()

#protein to everything,colsums first
#s1_t colsum
s1_t = [([0]*disease_num) for i in range(protein_num)]
s1_t_colsum = []
for i in range(protein_num):
	for j in range(disease_num):
		s1_t[i][j] = s1[j][i]

for i in range(protein_num):
	colsum = sum(s1_t[i])
	s1_t_colsum.append(colsum)

fromfile4.close()
del s1

#s3 colsum
fromfile5 = open(file5,"r")
s3 = [([0]*protein_num) for i in range(protein_num)]
s3_colsum = []
for line in fromfile5:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = proteins.index(temp[0])
		row_num = proteins.index(temp[1])
		s3[col_num][row_num] = float(temp[2])
		s3[row_num][col_num] = float(temp[2])

for i in range(protein_num):
	colsum = sum(s3[i])
	s3_colsum.append(colsum)

#s2_t colsum
fromfile6 = open(file6,"r")
s2_t = [([0]*drug_num) for i in range(protein_num)]
s2_t_colsum = []
for line in fromfile6:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = drugs.index(temp[0])
		row_num = proteins.index(temp[1])
		s2_t[row_num][col_num] = 1

for i in range(protein_num):
	colsum = sum(s2_t[i])
	s2_t_colsum.append(colsum)

#--------------------------------s6_drug_disease----------------------------
#s6
fromfile11 = open(file11, "r")
s6 = [([0]*disease_num) for i in range(drug_num)]
for line in fromfile11:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = diseases.index(temp[1])
		row_num = drugs.index(temp[0])
		s6[row_num][col_num] = 1.0
fromfile11.close()

#s6_t
fromfile11 = open(file11, "r")
s6_t = [([0]*drug_num) for i in range(disease_num)]
for line in fromfile11:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = drugs.index(temp[0])
		row_num = diseases.index(temp[1])
		s6_t[row_num][col_num] = 1.0

fromfile11.close()
#---------------------------------------------------------------------------
#m2
m2 = [([0]*disease_num) for i in range(protein_num)]
for i in range(protein_num):
	colsum1 = s1_t_colsum[i]
	colsum2 = s3_colsum[i]
	colsum3 = s2_t_colsum[i]

	if colsum1 == 0:
		for j in range(disease_num):
			m2[i][j] = 0
	elif (colsum2 + colsum3 == 0):
		for j in range(disease_num):
			m2[i][j] = s1_t[i][j]/colsum1
	else:
		for j in range(disease_num):
			m2[i][j] = a * s1_t[i][j]/colsum1

m2 = np.matrix(m2)
del s1_t
infile2.write("Complete constructing m2\n")
infile2.flush()

#m3
m3 = [([0]*protein_num) for i in range(protein_num)]
for i in range(protein_num):
	colsum1 = s1_t_colsum[i]
	colsum2 = s3_colsum[i]
	colsum3 = s2_t_colsum[i]

	if colsum2 == 0:
		for j in range(protein_num):
			m3[i][j] = 0
	elif (colsum1 + colsum3 == 0):
		for j in range(protein_num):
			m3[i][j] = s3[i][j]/colsum2
	else:
		for j in range(protein_num):
			m3[i][j] = b * s3[i][j]/colsum2

m3 = np.matrix(m3)
del s3
infile2.write("Complete constructing m3\n")
infile2.flush()

#m4
m4 = [([0]*drug_num) for i in range(protein_num)]
for i in range(protein_num):
	colsum1 = s1_t_colsum[i]
	colsum2 = s3_colsum[i]
	colsum3 = s2_t_colsum[i]

	if colsum3 == 0:
		for j in range(drug_num):
			m4[i][j] = 0
	elif (colsum1 + colsum2 == 0):
		for j in range(drug_num):
			m4[i][j] = s2_t[i][j]/colsum3
	else:
		for j in range(drug_num):
			m4[i][j] = 2 * s2_t[i][j]/colsum3

m4 = np.matrix(m4)
infile2.write("Complete constructing m4\n")
infile2.flush()
del s1_t_colsum
del s3_colsum
del s2_t_colsum

#m5
s2 = [([0]*protein_num) for i in range(drug_num)]
m5 = [([0]*protein_num) for i in range(drug_num)]
for i in range(drug_num):
	for j in range(protein_num):
		s2[i][j] = s2_t[j][i]

for i in range(drug_num):
	colsum = sum(s2[i])
	for j in range(protein_num):
		m5[i][j] = b * s2[i][j]/colsum
m5 = np.matrix(m5)
infile2.write("Complete constructing m5\n")
infile2.flush()
del s2_t

#m6
s4 = [([0]*drug_num) for i in range(drug_num)]
m6 = [([0]*drug_num) for i in range(drug_num)]
fromfile7 = open(file7,"r")
for line in fromfile7:
	if line.split():
		temp = line[:-1].split("\t")
		col_num = drugs.index(temp[0])
		row_num = drugs.index(temp[1])
		s4[col_num][row_num] = float(temp[2])
		s4[row_num][col_num] = float(temp[2])

for i in range(drug_num):
	colsum = sum(s4[i])
	for j in range(drug_num):
		m6[i][j] = c * s4[i][j]/colsum

m6 = np.matrix(m6)
infile2.write("Complete constructing m6\n")
infile2.flush()
del s4

#-------------------------m8-------and-------m9----------------------
m8 = [([0]*drug_num) for i in range(disease_num)]
for i in range(disease_num):
	rowsum = sum(s6_t[i])
	for j in range(drug_num):
		if rowsum == 0:
			m8[i][j] = 0
		else:
		 	m8[i][j] = c * s6_t[i][j]/rowsum

m9 = [([0]*disease_num) for i in range(drug_num)]
for i in range(drug_num):
	rowsum = sum(s6[i])
	for j in range(disease_num):
		if rowsum == 0:
			m9[i][j] = 0
		else:
			m9[i][j] = a * s6[i][j]/rowsum
#---------------------------------------------------------------------

#make M
hang1 = hstack((m7,m1,m8))
del m7,m1,m8
gc.collect()

hang2 = hstack((m2,m3,m4))
del m2,m3,m4
gc.collect()

kong3 = m9
hang3 = hstack((m9,m5,m6))
del m9,m5,m6

M = vstack((hang1,hang2,hang3))
del hang1,hang2,hang3
gc.collect()

#normalization
Mc = M
M = (Mc - np.min(Mc))/(np.max(Mc)-np.min(Mc))
del Mc

infile2.write("M has normalized\n")
infile2.flush()

#start from drug
for i in range(drug_num):
	infile2.write(str(i)+"\n")
	infile2.flush()
	choose_drug = i
	choose_target = get_target(choose_drug)

	u = [([0]) for i in range(disease_num)]
	u = np.matrix(u)

	v = [([0]) for i in range(protein_num)]
	for i in choose_target:
		v[proteins.index(i)][0] = 1
	v = np.matrix(v)

	h = [([0]) for i in range(drug_num)]
	h[choose_drug][0] = 1
	h = np.matrix(h)

	#calculate
	p0 = vstack((a*u,b*v,c*h))

	p1 = p0
	p2 = (1-r)*M.T*p1 + r*p0
	
	j = 0

	#if not converge in 10000, break
	while(L1_norm(p1,p2) > 10**(-10)):
		p3 = (1-r)*M.T*p2 + r*p0
		p1 = p2
		p2 = p3
		j = j+1
		if (j == 10000):
			break

	infile2.write(str(j)+"\n")
	p2 = p2.tolist()
	infile2.write("Start writing result\n")
	infile2.flush()
	infile1.write(str(p2[0])[1:-1])
	for j in range(drug_num + protein_num + disease_num - 1):
		infile1.write("\t"+str(p2[j+1])[1:-1])

	infile1.write("\n")

	infile2.write("Writing done\n")
	infile2.flush()

infile1.close()
infile2.close()




