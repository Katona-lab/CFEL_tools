#!/usr/bin/python
#
# (c) 2019 Gergely Katona <gergely.katona@gu.se>

import io
import pandas as pd
from cctbx.array_family import flex
from cctbx import crystal
from cctbx import uctbx
from cctbx import sgtbx
from cctbx import miller
import scipy as sp
import scipy.stats
from iotbx import reflection_file_reader
import numpy as np
import gc
from multiprocessing import Pool

unit_cell=(226.500,  226.500,  113.900, 90.0, 90.0, 90.0)
uc = uctbx.unit_cell(unit_cell)

wavelength = 1.89

xtal_symm = crystal.symmetry(
        unit_cell=unit_cell,
        space_group_symbol="P 43 21 2")
msnam = miller.build_set(crystal_symmetry=xtal_symm,anomalous_flag=False, d_min=2.3)

j=0
drows=[]
chunkbatch=5000

def wr_pkl(drows,no):
    global state
    df=pd.DataFrame(drows)                            
    df['I']=df['I'].astype('float32')
    df['sigI']=df['sigI'].astype('float32')

    miller_indices = flex.miller_index()
    data = flex.double()
    sigmas = flex.double()
    batch = flex.int()

    miller_indices=flex.miller_index(list(df['hkl'].values))
    data=flex.double(list(df['I'].values))
    sigmas=flex.double(list(df['sigI'].values))

    miller_set=miller.set(crystal_symmetry=xtal_symm,indices=miller_indices,anomalous_flag=False)

    miller_array = miller_set.array(data=data,sigmas=sigmas).set_observation_type_xray_intensity()
    mapped=miller_array.map_to_asu()

    df['hkl_asu']=list(mapped.indices())
    df['hkl_asu']=df['hkl_asu'].astype('category')
    df=df.drop(columns='hkl')
    df.to_pickle(state+"_%03i.pkl" % (no))
    return


with open(root+stream) as f:
    while True:
        line=f.readline()
        if line=="----- Begin chunk -----\n" :
            while True:
                line=f.readline()
                if "End chunk" in line:
                    j=j+1
                    if j % chunkbatch == 0:
                        wr_pkl(drows,j/chunkbatch)
                        print (j)
                        last=j/chunkbatch
                        drows=[]
                    break
                if "Begin crystal" in line:
                    while True:
                        line=f.readline()
                        if "End crystal" in line:
                            break
                        if "sigma(I)" in line:
                            while True:
                                drow={}
                                line=f.readline()
                                if "End" in line:
                                    break
                                sline=line.split()
                                drow['hkl']=(int(sline[0]),int(sline[1]),int(sline[2]))
                                drow['I']=float(sline[3])
                                drow['sigI']=float(sline[4])

                                drows.append(drow)                         
        if not line:
            wr_pkl(drows,last+1)
            break

                
count=0
def looker(etind):
    global count
    drows=[]
    sub=cf[cf.hkl_asu==etind]
    print (str(etind)+" "+str(count))
    count=count+1
    if len(sub['I'])>2:
        Is=np.array(sub['I'])
        lenIs=len(sub['I'])
        mean=np.mean(Is)
        sem=sp.stats.sem(Is)
        drow={}
        drow['I']=mean
        drow['SIGI']=sem
        drow['batch']=0
        drow['hkl']=etind
        drows.append(drow)
        for i in range(1,100):
            drow={}
            resam=np.random.choice(Is, size=lenIs)
            mean=np.mean(resam)
            sem=sp.stats.sem(resam)
            drow['I']=mean
            drow['SIGI']=sem
            drow['batch']=i
            drow['hkl']=etind
            drows.append(drow)
        gf=pd.DataFrame(drows)
        gf['I']=gf['I'].astype('float32')
        gf['SIGI']=gf['SIGI'].astype('float32')
        gf['batch']=gf['batch'].astype('int16')
        return gf
    return pd.DataFrame()
        

def resample(i):
    global Is,lenIs
    drow={}
    resam=np.random.choice(Is, size=lenIs)
    mean=np.mean(resam)
    sem=sp.stats.sem(resam)
    drow['I']=mean
    drow['SIGI']=sem
    drow['batch']=i
    drow['hkl']=etind
    return drow

p=Pool(processes=24)
resa=p.map(looker,list(msnam.indices())[:])


rf=pd.concat(resa)
rf.to_pickle('LCLS2015_'+state+'_bootstrapped.pkl')
