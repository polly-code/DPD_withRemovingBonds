import numpy as np
from scipy.spatial import ConvexHull
from PyGEL3D import gel
import matplotlib.pyplot as plt

class bead:
    ''' class bead contains all info about bead in chain: global number, 
    the remaining valence, type of the bead and coordinates '''
    numbd=int()
    valence=int()
    typep=int()
    x=float()
    y=float()
    z=float()
    neighbors=[]
    def __lt__(self, other):
        return self.numbd < other.numbd
    
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
            if one.typep==2:# or one.typep==2:
                polymer.bd.append(one)
        elif 'angles' in line:
            print ('Reading finished')
        elif bnd:
            sb=bond()
            head,tail=line.split()
            sb.first=int(head)
            sb.last=int(tail)
            polymer.bnd.append(sb)
            
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
            
def apply_epi(polymer, path='path/to/epigenetics/colors_bybin_chrX.tsv'):
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
                
def dist(hull, points):
    # Construct PyGEL Manifold from the convex hull
    m = gel.Manifold()
    for s in hull.simplices:
        m.add_face(hull.points[s])

    dist = gel.MeshDistance(m)
    res = []
    for p in points:
        # Get the distance to the point
        # But don't trust its sign, because of possible
        # wrong orientation of mesh face
        d = dist.signed_distance(p)

        # Correct the sign with ray inside test
        if dist.ray_inside_test(p):
            if d > 0:
                d *= -1
        else:
            if d < 0:
                d *= -1
        res.append(d)
    return np.array(res)

def writeMol2 (polymer, path):
    bstr='1  ala'
    the_file=open(path, 'w')
    the_file.write('@<TRIPOS>MOLECULE\n')
    the_file.write('mol_name\n')
    the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(len(polymer.bd), len(polymer.bnd), '0', '0', '0'))
    #the_file.write('\t %d \t %d \t %s \t %s \t %s \n' %(polymer.number_of_beads, polymer.number_of_beads-1, '0', '0', '0'))
    the_file.write('SMALL\n')
    the_file.write('USER_CHARGES\n')
    the_file.write('@<TRIPOS>ATOM\n')
    for i in range(len(polymer.bd)):

        if polymer.bd[i].typep==1:
            ty='C'
        elif polymer.bd[i].typep==2:
            ty='O'
        elif polymer.bd[i].typep==3:
            ty='S'
        elif polymer.bd[i].typep==4:
            ty='N'
        else:
            print('Some strange type.. ', polymer.bd[i].typep)
        the_file.write('%d \t %s \t %f \t %f \t %f \t %s \t %s \t %f \n' %(i+1, ty, polymer.bd[i].x, polymer.bd[i].y, polymer.bd[i].z, ty, bstr, float(i)))
    the_file.write('@<TRIPOS>BOND\n')
    k=0
    for i in range(len(polymer.bnd)):
        
        if abs(polymer.bnd[i].first-polymer.bnd[i].last)>0:        
            the_file.write('%d \t %d \t %d \t %s \n' %(k+1, polymer.bnd[i].first, polymer.bnd[i].last, '1'))
            k+=1
#        elif (polymer.bnd[i].first>799 and polymer.bnd[i].first<901) and (polymer.bnd[i].last>799 and polymer.bnd[i].last<901):
#            the_file.write('%d \t %d \t %d \t %s \n' %(k+1, polymer.bnd[i].first, polymer.bnd[i].last, '2'))
#            k+=1
    the_file.close()
#mat = np.random.rand(100, 3)
#hull = ConvexHull(mat)
#points = np.random.rand(10, 3)
#print(dist(hull, points))
def calc_cm(coords, part, path):
    cm=(np.average(coords[:,0]),np.average(coords[:,1]),np.average(coords[:,2]))
    dst2cm=[]
    print(cm)
    for i in range(len(part)):
        dst2cm.append(np.sqrt((part[i][0]-cm[0])**2+(part[i][1]-cm[1])**2+(part[i][2]-cm[2])**2))
    np.savetxt(path,dst2cm)
    plt.boxplot(dst2cm)
    plt.show()
    
def calc_cm_dif(coords, part):
    cm=(np.average(coords[:,0]),np.average(coords[:,1]),np.average(coords[:,2]))
    cmp=(np.average(part[:,0]),np.average(part[:,1]),np.average(part[:,2]))
    print(cm, cmp)
    return (np.sqrt((cmp[0]-cm[0])**2+(cmp[1]-cm[1])**2+(cmp[2]-cm[2])**2))

def distance (a, b):
    '''calculate distance between two beads'''
    return np.sqrt((a.x-b.x)**2+(a.y-b.y)**2+(a.z-b.z)**2)

def cm_epi (polymer):
    cm_inact=[0,0,0,0]
    cm_act=[0,0,0,0]
    cm=[0,0,0,0]
    f=open('C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/sizes_chr_terr.dat','a')
    for i in range(len(polymer.bd)):
        if polymer.bd[i].typep==1:
            cm_act[0]+=polymer.bd[i].x
            cm_act[1]+=polymer.bd[i].y
            cm_act[2]+=polymer.bd[i].z
            cm_act[3]+=1
        if polymer.bd[i].typep==2:
            cm_inact[0]+=polymer.bd[i].x
            cm_inact[1]+=polymer.bd[i].y
            cm_inact[2]+=polymer.bd[i].z
            cm_inact[3]+=1
        cm[0]+=polymer.bd[i].x
        cm[1]+=polymer.bd[i].y
        cm[2]+=polymer.bd[i].z
        cm[3]+=1
    a=bead()
    b=bead()
    c=bead()
    a.x,a.y,a.z=cm_act[0]/cm_act[3], cm_act[1]/cm_act[3], cm_act[2]/cm_act[3]
    b.x,b.y,b.z=cm_inact[0]/cm_inact[3], cm_inact[1]/cm_inact[3], cm_inact[2]/cm_inact[3]
    c.x,c.y,c.z=cm[0]/cm[3], cm[1]/cm[3], cm[2]/cm[3]
    
    print ("center of mass active (", cm_act[3]," beads): ", cm_act[0]/cm_act[3], cm_act[1]/cm_act[3], cm_act[2]/cm_act[3])
    print ("center of mass inactive (", cm_inact[3]," beads): ", cm_inact[0]/cm_inact[3], cm_inact[1]/cm_inact[3], cm_inact[2]/cm_inact[3])
    print ("center of mass (", cm[3]," beads): ", cm[0]/cm[3], cm[1]/cm[3], cm[2]/cm[3])
    print("distance between CoMs active - inactive is: ",distance(a,b))
    print("distance between CoMs active - all is: ",distance(a,c))
    print("distance between CoMs inactive - all is: ",distance(c,b))
    
    #----
    #relate distance
    #----
    ra=bead()
    rb=bead()
    ra.x,ra.y,ra.z=a.x-c.x,a.y-c.y,a.z-c.z
    rb.x,rb.y,rb.z=b.x-c.x,b.y-c.y,b.z-c.z
    print ("relative center of mass active (", cm_act[3]," beads): ", ra.x,ra.y,ra.z)
    print ("relative center of mass inactive (", cm_inact[3]," beads): ", rb.x,rb.y,rb.z)
    #-------------------
    
    maxdist=0
    arr_dist=[]
    for i in range(len(polymer.bd)):
        ndist=distance(polymer.bd[i], c)
        arr_dist.append(ndist)
        if maxdist<ndist:
            maxdist=ndist
    arr=np.array(arr_dist)
    arr.sort()
    print("First and last distances: ",arr[0], arr[-1])
    avermax=0
    arr_desc=arr[::-1]
    for i in range(100):
        avermax+=arr_desc[i]
    #print("Avermax value: ", avermax/100)
    f.write('%f\n'%float(avermax/100))
    print("maximum distance to CoM all: ", maxdist)
    f.close()

def main():
    cell=['a2', 'a3', 'a5', 'a6', 'a8', 'a9' ,'b3','b6','b15','b16','b19','b26','b31','sc1','sc16','sc19','sc21','sc23','sc24','sc29']
    #cell=['a2','b16','b19','b31']
    bpath = 'path/to/testarts/'
    for i in range(len(cell)):
        chr_act=[]
        chr_ina=[]
        chr_neu=[]
        chr_pol=[]
        chr_all=[]
        path=bpath+cell[i] #+"//"+'restart3.dat'
        f0=open(path+'/restart.dat')
        polymer0=chain()
        read_rst(f0,polymer0)
        polymer0.bd.sort()
        removePBC (polymer0)
        apply_epi(polymer0)
        for j in range(len(polymer0.bd)):
            if polymer0.bd[j].typep==1:
                chr_act.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==2:
                chr_ina.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==3:
                chr_neu.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==4:
                chr_pol.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
        chr_all.extend(chr_act)
        chr_all.extend(chr_ina)
        chr_all.extend(chr_neu)
        chr_all.extend(chr_pol)
        
        nchr_act=np.array(chr_act)
        nchr_ina=np.array(chr_ina)
        nchr_neu=np.array(chr_neu)
        nchr_pol=np.array(chr_pol)
        nchr_all=np.array(chr_all)
        hull = ConvexHull(nchr_all)
        print('active ',len(nchr_act))
        da=dist(hull, nchr_act)
        di=dist(hull, nchr_ina)
        dn=dist(hull, nchr_neu)
        dp=dist(hull, nchr_pol)
        #np.savetxt(bpath+cell[i]+'a.dat',da)
        #np.savetxt(bpath+cell[i]+'i.dat',di)
        #np.savetxt(bpath+cell[i]+'n.dat',dn)
        np.savetxt(bpath+cell[i]+'p.dat',dp)

        data=np.array((da,di,dn,dp))
        plt.boxplot(-1*data)
        plt.xticks([1, 2, 3, 4], ['active', 'inactive', 'neutral', 'polycomb'])
        plt.savefig(bpath+cell[i]+'_bp.png', dpi = 300)
        plt.show()


def main2():
    cell=['a2', 'a3', 'a5', 'a6', 'a8', 'a9' ,'b3','b6','b15','b16','b19','b26','b31','sc1','sc16','sc19','sc21','sc23','sc24','sc29']
    #cell=['a6','b16','b19','b31']
    bpath = 'path/to/restarts/'
    #f=open('C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/all_data_z0_z4_others_z0-z4_dst2cm.dat', 'w')
    for i in range(len(cell)):
        path=bpath+cell[i]+'/restart3.dat'
        f0=open(path)
        polymer0=chain()
        read_rst(f0,polymer0)
        polymer0.bd.sort()
        removePBC (polymer0)
        apply_epi_expression(polymer0,'C:/Kos/Analysis/Bio/single_cell/best_res/epigenetics/table_expression.tsv')
        cm_epi(polymer0)
        zone0=[]
        zone4=[]
        others=[]
        chr_all=[]
        for j in range(len(polymer0.bd)):
            if polymer0.bd[j].typep==1:
                zone0.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==2:
                zone4.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==3:
                others.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
        chr_all.extend(zone0)
        chr_all.extend(zone4)
        chr_all.extend(others)
        
        nchr_zone0=np.array(zone0)
        nchr_zone4=np.array(zone4)
        nchr_others=np.array(others)

        nchr_all=np.array(chr_all)
        #f.write("%f\t%f\t%f\t%f\n"%(calc_cm_dif(nchr_all,nchr_zone0),calc_cm_dif(nchr_all,nchr_zone4),calc_cm_dif(nchr_all,nchr_others),calc_cm_dif(nchr_zone0,nchr_zone4)))
    #f.close()
        #calc_cm_dif(nchr_zone0,nchr_zone4,'c')
        #calc_cm_dif(nchr_all,nchr_zone0,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_z0z4others_dst2cm.dat')
        #calc_cm_dif(nchr_all,nchr_zone4,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_z0z4others_dst2cm.dat')
        #calc_cm_dif(nchr_all,nchr_others,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_z0z4others_dst2cm.dat')
        #calc_cm(nchr_all,zone0,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_z0_dst2cm.dat')
        #calc_cm(nchr_all,zone4,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_z4_dst2cm.dat')
        #calc_cm(nchr_all,others,'C:/Kos/Analysis/Bio/single_cell/best_res/dpd/10k/'+cell[i]+'_others_dst2cm.dat')
        hull = ConvexHull(nchr_all)
        print('active ',len(nchr_zone0))
        da=dist(hull, nchr_zone0)
        di=dist(hull, nchr_zone4)
        dn=dist(hull, nchr_others)
        np.savetxt(bpath+cell[i]+'z0_expr.dat',-1*da)
        np.savetxt(bpath+cell[i]+'z4_expr.dat',-1*di)
        np.savetxt(bpath+cell[i]+'others_expr.dat',-1*dn)
        data=np.array((da,di,dn))
        plt.boxplot(-1*data)
        plt.xticks([1, 2, 3], ['1', '4', 'etc'])
        plt.show()
#main2()

def main3():
    cell=['a2', 'a3', 'a5', 'a6', 'a8', 'a9' ,'b3','b6','b15','b16','b19','b26','b31','sc1','sc16','sc19','sc21','sc23','sc24','sc29']
    #cell=['a6','b16','b19','b31']
    bpath = 'path/to/restarts/'

    for i in range(len(cell)):
        path=bpath+cell[i]+'/restart.dat'
        f0=open(path)
        polymer0=chain()
        read_rst(f0,polymer0)
        polymer0.bd.sort()
        removePBC (polymer0)
        apply_epi_expression(polymer0,'path/to/epigenetics/table_expression.tsv')
        cm_epi(polymer0)
        zone0=[]
        zone4=[]
        others=[]
        chr_all=[]
        for j in range(len(polymer0.bd)):
            if polymer0.bd[j].typep==1:
                zone0.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==2:
                zone4.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
            elif polymer0.bd[j].typep==3:
                others.append((polymer0.bd[j].x,polymer0.bd[j].y,polymer0.bd[j].z))
        chr_all.extend(zone0)
        chr_all.extend(zone4)
        chr_all.extend(others)
        
        nchr_zone0=np.array(zone0)
        nchr_zone4=np.array(zone4)
        nchr_others=np.array(others)

        nchr_all=np.array(chr_all)
        hull = ConvexHull(nchr_all)
        #f.write("%f\t%f\t%f\t%f\n"%(calc_cm_dif(nchr_all,nchr_zone0),calc_cm_dif(nchr_all,nchr_zone4),calc_cm_dif(nchr_all,nchr_others),calc_cm_dif(nchr_zone0,nchr_zone4)))
        da=dist(hull, nchr_zone0)
        di=dist(hull, nchr_zone4)
        dn=dist(hull, nchr_others)
        np.savetxt(bpath+cell[i]+'z0_expr.dat',-1*da)
        np.savetxt(bpath+cell[i]+'z4_expr.dat',-1*di)
        np.savetxt(bpath+cell[i]+'others_expr.dat',-1*dn)
    #f.close()