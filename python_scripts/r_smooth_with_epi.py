#!/bin/python
import numpy as np
class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    numbd=int()
    valence=int()
    typep=int()
    x=float()
    y=float()
    z=float()
    
class bond:
    '''class bond contains all info about bond: which beads connected by this bond'''
    first=int()
    last=int()
    
class chain:
    '''class chain has two lists of beads and bonds and general info about system such as total number of particles, density and box size along each axis'''
    def __init__(self):
        self.bd=[]
        self.bnd=[]
        self.number_of_beads=float()
        self.number_of_bonds=float()
        self.density=float()
        self.xbox=float()
        self.ybox=float()
        self.zbox=float()
    
def read_rst (f, polymer):
    '''function to read the restart file to the class chain, get path to file and name of chain'''
    one=bead()
    sb=bond()
    bnd=False
    for i,line in enumerate(f):
        if i==0:
            head,tail=line.split()
            polymer.number_of_beads = float(head)
            polymer.density = float(tail)
        elif i==1:
            xbox, ybox, zbox=line.split()
            polymer.xbox=float(xbox)
            polymer.ybox=float(ybox)
            polymer.zbox=float(zbox)
        elif 'bonds:' in line:
            bnd=True
        elif i>1 & bnd==False:
            one=bead()
            numbd,valence,typep,x,y,z=line.split()
            one.numbd=int(numbd)
            one.valence=int(valence)
            one.typep=int(typep)
            one.x=float(x)
            one.y=float(y)
            one.z=float(z)
            if one.typep==2:
                polymer.bd.append(one)
        elif 'angles' in line:
            print ('Reading finished')
        elif bnd:
            sb=bond()
            head,tail=line.split()
            sb.first=int(head)
            sb.last=int(tail)
            polymer.bnd.append(sb)

def sortcl (polymer):
    '''function sorts the chain elements, only beads'''
    emp=bead();
    for i in range(len(polymer.bd)-1):
        for j in range (len(polymer.bd)-i-1):
            if polymer.bd[j].numbd > polymer.bd[j+1].numbd:
                #polymer.bd[i], polymer.bd[j] = polymer.bd[j], polymer.bd[i]
                emp=polymer.bd[j]
                polymer.bd[j]=polymer.bd[j+1]
                polymer.bd[j+1]=emp

def removePBC (polymer):
    '''revome periodic boundary conditions in case of single chain and numbers of beads corresspond to the global numbers and start from 1'''
    itx=0
    ity=0
    itz=0
    for i in range(len(polymer.bd)-1):
        if polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x > polymer.xbox / 2:
            itx = itx + 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        elif polymer.bd[i].x - itx * polymer.xbox - polymer.bd[i+1].x < -polymer.xbox / 2:
            itx = itx - 1
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
        else:
            polymer.bd[i+1].x = polymer.bd[i+1].x + itx*polymer.xbox
            
        if polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y > polymer.ybox / 2:
            ity = ity + 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        elif polymer.bd[i].y - ity * polymer.ybox - polymer.bd[i+1].y < -polymer.ybox / 2:
            ity = ity - 1
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
        else:
            polymer.bd[i+1].y = polymer.bd[i+1].y + ity*polymer.ybox
            
        if polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z > polymer.zbox / 2:
            itz = itz + 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        elif polymer.bd[i].z - itz * polymer.zbox - polymer.bd[i+1].z < -polymer.zbox / 2:
            itz = itz - 1
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz*polymer.zbox
        else:
            polymer.bd[i+1].z = polymer.bd[i+1].z + itz * polymer.zbox
            
def coarsing(pol1, pol2, n=11):
    
    for i in range(len(pol1.bd)):
        one=bead()
        one.x=0
        one.y=0
        one.z=0
        pol2.bd.append(one)
        for j in range(max(i-int((n-1)/2),0),min(i+int((n-1)/2),len(pol1.bd))):
            pol2.bd[i].x += pol1.bd[j].x
            pol2.bd[i].y += pol1.bd[j].y
            pol2.bd[i].z += pol1.bd[j].z
        pol2.bd[i].x /= min(i+int((n-1)/2),len(pol1.bd))-max(i-int((n-1)/2),0)
        pol2.bd[i].y /= min(i+int((n-1)/2),len(pol1.bd))-max(i-int((n-1)/2),0)
        pol2.bd[i].z /= min(i+int((n-1)/2),len(pol1.bd))-max(i-int((n-1)/2),0)
    for i in range(len(pol1.bnd)):
        pol2.bnd.append(pol1.bnd[i])

def writeMol2 (polymer, path):
    bstr='1  ala'
    #the_file=open(path, 'w')
    with open(path, 'w') as the_file:
        the_file.write('@<TRIPOS>MOLECULE\n')
    with open(path, 'a') as the_file:
        the_file.write('mol_name\n')
    with open(path, 'a') as the_file:
        the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(len(polymer.bd), len(polymer.bd)-1, '0', '0', '0'))
    with open(path, 'a') as the_file:
        the_file.write('SMALL\n')
    with open(path, 'a') as the_file:
        the_file.write('USER_CHARGES\n')
    with open(path, 'a') as the_file:
        the_file.write('@<TRIPOS>ATOM\n')
    ty='C'
    for i in range(len(polymer.bd)):
        if polymer.bd[i].typep==1:
            ty='C'
        elif polymer.bd[i].typep==2:
            ty='O'
        elif polymer.bd[i].typep==3:
            ty='S'
        elif polymer.bd[i].typep==4:
            ty='N'
        with open(path, 'a') as the_file:
            the_file.write('%d \t %s \t %f \t %f \t %f \t %s \t %s \t %f \n' %(i+1, ty, polymer.bd[i].x, polymer.bd[i].y, polymer.bd[i].z, ty, bstr, float(i)))
    with open(path, 'a') as the_file:
        the_file.write('@<TRIPOS>BOND\n')
    newit=1
    for i in range(len(polymer.bnd)):
        if abs(polymer.bnd[i].first-polymer.bnd[i].last)==1:
            with open(path, 'a') as the_file:
                the_file.write('%d \t %d \t %d \t %s \n' %(newit, polymer.bnd[i].first, polymer.bnd[i].last, '1'))
            newit+=1
def f_writeMol2 (polymer, path):
    bstr='1  ala'
    the_file=open(path, 'w')
    the_file.write('@<TRIPOS>MOLECULE\n')
    the_file.write('mol_name\n')
    the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(len(polymer.bd), len(polymer.bd)-1, '0', '0', '0'))
    the_file.write('SMALL\n')
    the_file.write('USER_CHARGES\n')
    the_file.write('@<TRIPOS>ATOM\n')
    ty='C'
    for i in range(len(polymer.bd)):
        if polymer.bd[i].typep==1:
            ty='C'
        elif polymer.bd[i].typep==2:
            ty='O'
        elif polymer.bd[i].typep==3:
            ty='S'
        elif polymer.bd[i].typep==4:
            ty='N'

        the_file.write('%d \t %s \t %f \t %f \t %f \t %s \t %s \t %f \n' %(i+1, ty, polymer.bd[i].x, polymer.bd[i].y, polymer.bd[i].z, ty, bstr, float(i)))
    the_file.write('@<TRIPOS>BOND\n')
    newit=1
    for i in range(len(polymer.bnd)):
        if abs(polymer.bnd[i].first-polymer.bnd[i].last)==1:
            the_file.write('%d \t %d \t %d \t %s \n' %(newit, polymer.bnd[i].first, polymer.bnd[i].last, '1'))
            newit+=1
    the_file.close()

def apply_epi(polymer, path="path/to/annotation/colors_bybin_chrX.tsv"):
    f=open(path)
    for i,line in enumerate(f):
        if i==0:
            print ("epi start")
        if i!=0:
            chr,num,epitype=line.split()
            num=int(num)
            if epitype=="active":
                polymer.bd[num].typep=1
            elif epitype=="inactive":
                polymer.bd[num].typep=2
            elif epitype=="neutral":
                polymer.bd[num].typep=3
            elif epitype=="polycomb":
                polymer.bd[num].typep=4

def apply_epi_expression(polymer, path):
    f=open(path)
    for i,line in enumerate(f):
        if i==0:
            print ("epi start")
        if i!=0:
            a,b,c,epitype,e=line.split()
            if epitype=="0":
                polymer.bd[i-1].typep=1
            elif epitype=="4":
                polymer.bd[i-1].typep=2
            elif epitype=="1" or epitype=="2" or epitype=="3":
                polymer.bd[i-1].typep=3
                
def apply_epi_expression_biased(polymer, path, bias):
    f=open(path)
    for i,line in enumerate(f):
        if i==0:
            print ("epi start")
        elif i>bias:
            a,b,c,epitype,e=line.split()
            if epitype=="0":
                polymer.bd[i-1-bias].typep=1
            elif epitype=="4":
                polymer.bd[i-1-bias].typep=2
            elif epitype=="1" or epitype=="2" or epitype=="3":
                polymer.bd[i-1-bias].typep=3
    for i in range(len(polymer.bd)-bias,len(polymer.bd)):
        polymer.bd[i].typep=3
        
                
def unit_vector(vector):
    """ Returns the unit vector of the vector.  """
    return vector / np.linalg.norm(vector)

def angle_between(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'::

            >>> angle_between((1, 0, 0), (0, 1, 0))
            1.5707963267948966
            >>> angle_between((1, 0, 0), (1, 0, 0))
            0.0
            >>> angle_between((1, 0, 0), (-1, 0, 0))
            3.141592653589793
    """
    v1_u = unit_vector(v1)
    v2_u = unit_vector(v2)
    return np.arccos(np.clip(np.dot(v1_u, v2_u), -1.0, 1.0))

def main():
    pa = 'path/to/restarts/'
    #cells=['a5-b19', 'a5-b26', 'a6-a8', 'a6-b31', 'a8-b26', 'a8-b31', 'b31-sc23']
    #cells=['a5', 'a6', 'a8', 'b19', 'b26', 'b31', 'sc23']
    biases=[2,3,2,1,4,2,1,2,6,1,1,2,1,4,3,2,23,0,3,29]
    cells=['a2', 'a3', 'a5', 'a6', 'a8', 'a9', 'b3', 'b6', 'b15', 'b16', 'b19', 'b26', 'b31', 'sc1', 'sc16', 'sc19', 'sc21', 'sc23', 'sc24', 'sc29']
    for i in range(len(cells)):
        path =pa+cells[i]+'/restart3.dat'
        f0=open(path)
        polymer0=chain()
        pol2=chain()
        read_rst(f0,polymer0)
        sortcl (polymer0)
        removePBC (polymer0)
        coarsing(polymer0, pol2, n=15)
        #apply_epi(pol2, "C:/Kos/Analysis/Bio/single_cell/best_res/epigenetics/colors_bybin_chrX.tsv")
        apply_epi_expression_biased(pol2, "path/to/epigenetics/table_expression.tsv",biases[i])
        f_writeMol2(pol2, pa+'bee_'+cells[i]+'.mol2')
main()