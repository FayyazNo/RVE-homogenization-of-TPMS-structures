# Homogenization thermal properties of the voxelized RVE
# This Code is Written By Fayyaz Nosouhi (dehnavifn@gmail.com)
# Tehran University March-2020

from abaqus import *
from abaqusConstants import *
import __main__
import section
import regionToolset
import displayGroupMdbToolset as dgm
import part
import material
import assembly
import step
import interaction
import load
import mesh
import optimization
import job
import sketch
import visualization
import xyPlot
import displayGroupOdbToolset as dgo
import numpy as np
import time
from operator import itemgetter 
time_strat=time.asctime( time.localtime(time.time()) )

#Parameters---------------------------------------------------------------------
Nx=60
Ny=60
Nz=60
ax=1                # factor of the cell lenhth in x direction
ay=1
az=1
E0=.0000001         #thermal cond. of phase 0 (void)
v0=.1
E1=72000              # ther. cond. of phase 1 (solid)
v1=.33

CPU_No=2
ModelName="Schwarz0_1"
RVE_Type='Schwarz'     # RVE_Type:  IWP,   Schwarz,  Gyroid,  Diamond
Vol=30
cx=cy=cz=1         # cell number
jobName=ModelName
#---------------------------------------------------------
execfile('vol.py')
xyz=Nx*Ny*Nz
Lx=2*cx*pi
Ly=2*cy*pi
Lz=2*cz*pi

C=np.zeros([xyz])
c1=(Lx)/Nx
c2=(Ly)/Ny
c3=(Lz)/Nz
if   (RVE_Type== 'IWP'):
	vol_t=iwp_vol
elif (RVE_Type== 'Schwarz'):
	vol_t=sch_vol
elif (RVE_Type== 'Gyroid'):
	vol_t=gyr_vol
else  :  
	vol_t=dim_vol



	
t=0
Q=0
#
for q in range (len(vol_t)/2):
	if vol_t[q*2+1]>Vol:
		t=vol_t[2*q-2]
		Q=q-1
		break



#
while (True):
	C=np.zeros([xyz])
	for i in range (Nx):
		for j in range (Ny):
			for k in range (Nz):
				if (RVE_Type== 'IWP' or RVE_Type== 'Schwarz' or RVE_Type== 'Diamond'):
					x=(((i-Nx/2.)*c1)-.0*pi)
					y=(((j-Ny/2.)*c2)-.0*pi)
					z=(((k-Nz/2.)*c3)-.0*pi)
				elif (RVE_Type== 'Gyroid'):
					x=(((i-Nx/2.)*c1)-.28*pi)
					y=(((j-Ny/2.)*c2)-.28*pi)
					z=(((k-Nz/2.)*c3)-.28*pi)
				else:
					print("RVE_Type is not valid! select IWP,   Schwarz,  Gyroid or Diamond")
					break
				#
				if  (RVE_Type== 'Diamond'):
					a=sin(x)*sin(y)*sin(z)+sin(x)*cos(y)*cos(z)+cos(x)*sin(y)*cos(z)+cos(x)*cos(y)*sin(z) #DIAMOND
				elif(RVE_Type== 'Gyroid'):
					a=cos(x)*sin(y)+cos(y)*sin(z)+cos(z)*sin(x)  # Gyroide <1.2 & yzx-0.28*pi
				elif (RVE_Type== 'Schwarz'):
					a=cos(x)+cos(y)+cos(z) #schwartz
				else  :  
					a=(cos(2*x)+cos(2*y)+cos(2*z)-2*(cos(x)*cos(y)+cos(y)*cos(z)+cos(z)*cos(x))) #IWP
				#-------------------------------------------------------------------------------------
				if a<t :
					C[i*Ny*Nz+j*Nz+k]=1
	VF=100.0*np.sum(C)/len(C)
	if (abs(VF-Vol)<.1 or VF>Vol):
		if (RVE_Type== 'Schwarz'):
			C=np.fft.ifftshift((C.reshape(Nx, Ny, Nz)))
			C=C.reshape(Nx*Ny*Nz)
		break
	else:
		t=t+0.001
	print(VF, t)

#Generating a RVE size Nx*Ny*Nz(pixel size 1)-----------------------------------
Nx=Nx*ax 
Ny=Ny*ay 
Nz=Nz*az
mdb.Model(name=ModelName, modelType=STANDARD_EXPLICIT)
#os.chdir(os.getcwd()+dirName)
s = mdb.models[ModelName].ConstrainedSketch(name='__profile__', 
	sheetSize=200.0)
g, v, d, c = s.geometry, s.vertices, s.dimensions, s.constraints
s.setPrimaryObject(option=STANDALONE)
s.rectangle(point1=(0.0, 0.0), point2=(Nx, Ny))
p = mdb.models[ModelName].Part(name='Part-1', dimensionality=THREE_D, 
	type=DEFORMABLE_BODY)
p = mdb.models[ModelName].parts['Part-1']
p.BaseSolidExtrude(sketch=s, depth=Nz)
s.unsetPrimaryObject()
p = mdb.models[ModelName].parts['Part-1']
session.viewports['Viewport: 1'].setValues(displayedObject=p)
del mdb.models[ModelName].sketches['__profile__']
#Material---------------------------------------------------------------------
mdb.models[ModelName].Material(name='Material-0')
mdb.models[ModelName].materials['Material-0'].Elastic(table=((E0, v0), ))
mdb.models[ModelName].Material(name='Material-1')
mdb.models[ModelName].materials['Material-1'].Elastic(table=((E1, v1), ))

mdb.models[ModelName].HomogeneousSolidSection(name='Section-0', material='Material-0', thickness=None)
mdb.models[ModelName].HomogeneousSolidSection(name='Section-1', material='Material-1', thickness=None)
#Assembly--------------------------------------------------------------------
a = mdb.models[ModelName].rootAssembly
session.viewports['Viewport: 1'].setValues(displayedObject=a)
session.viewports['Viewport: 1'].assemblyDisplay.setValues(
	optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
a = mdb.models[ModelName].rootAssembly
a.DatumCsysByDefault(CARTESIAN)
p = mdb.models[ModelName].parts['Part-1']
a.Instance(name='Part-1-1', part=p, dependent=ON)
#Steps-----------------------------------------------------------------------
mdb.models[ModelName].StaticStep(name='Step-1', previous='Initial')
mdb.models[ModelName].steps['Step-1'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)
mdb.models[ModelName].StaticStep(name='Step-2', previous='Step-1')
mdb.models[ModelName].steps['Step-2'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)
mdb.models[ModelName].StaticStep(name='Step-3', previous='Step-2')
mdb.models[ModelName].steps['Step-3'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)
mdb.models[ModelName].StaticStep(name='Step-4', previous='Step-3')
mdb.models[ModelName].steps['Step-4'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)
mdb.models[ModelName].StaticStep(name='Step-5', previous='Step-4')
mdb.models[ModelName].steps['Step-5'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)
mdb.models[ModelName].StaticStep(name='Step-6', previous='Step-5')
mdb.models[ModelName].steps['Step-6'].setValues(matrixSolver=ITERATIVE,  matrixStorage=SYMMETRIC)

mdb.models[ModelName].fieldOutputRequests['F-Output-1'].setValues(variables=('S', 'LE'))

#boundary condition-------------------------------------------------------
a = mdb.models[ModelName].rootAssembly
f1 = a.instances['Part-1-1'].faces

#1- tension  direction x----------
faces1 = f1.findAt(((Nx, .001*Ny, 0.001*Nz),))
region = a.Set(faces=faces1, name='Set-2')
mdb.models[ModelName].DisplacementBC(name='BC-1', createStepName='Step-1', 
    region=region, u1=0.1, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
#
faces1 = f1.findAt(((0, .001*Ny, 0.001*Nz),))
region = a.Set(faces=faces1, name='Set-3')
mdb.models[ModelName].DisplacementBC(name='BC-2', createStepName='Step-1', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

mdb.models[ModelName].boundaryConditions['BC-1'].deactivate('Step-2')
mdb.models[ModelName].boundaryConditions['BC-2'].deactivate('Step-2')

# tension direction y-----------
faces1 = f1.findAt(((.001*Nx, Ny, 0.001*Nz),))
region = a.Set(faces=faces1, name='Set-4')
mdb.models[ModelName].DisplacementBC(name='BC-3', createStepName='Step-2', 
    region=region, u1=0.0, u2=0.1, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

faces1 = f1.findAt(((.001*Nx, 0, 0.001*Nz),))
region = a.Set(faces=faces1, name='Set-5')
mdb.models[ModelName].DisplacementBC(name='BC-4', createStepName='Step-2', 
    region=region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

mdb.models[ModelName].boundaryConditions['BC-3'].deactivate('Step-3')
mdb.models[ModelName].boundaryConditions['BC-4'].deactivate('Step-3')

# tension direction z-----------
faces1 = f1.findAt(((.001*Nx, 0.001*Ny, Nz),))
region = a.Set(faces=faces1, name='Set-6')
mdb.models[ModelName].DisplacementBC(name='BC-5', createStepName='Step-3', 
    region=region, u1=0.0, u2=0.0, u3=0.1, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)


faces1 = f1.findAt(((.001*Nx, 0.001*Ny, 0),))
region = a.Set(faces=faces1, name='Set-7')
mdb.models[ModelName].DisplacementBC(name='BC-6', createStepName='Step-3', 
    region=region, u1=0.0, u2=0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

mdb.models[ModelName].boundaryConditions['BC-5'].deactivate('Step-4')
mdb.models[ModelName].boundaryConditions['BC-6'].deactivate('Step-4')

# shear direction xy----------

faces1 = f1.findAt(((Nx, 0.001*Ny, .001*Nz),))
region = a.Set(faces=faces1, name='Set-8')
mdb.models[ModelName].DisplacementBC(name='BC-7', createStepName='Step-4', 
    region=region, u1=0.0, u2=0.05, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

faces1 = f1.findAt(((0, 0.001*Ny, .001*Nz),))
region = a.Set(faces=faces1, name='Set-9')
mdb.models[ModelName].DisplacementBC(name='BC-8', createStepName='Step-4', 
    region=region, u1=0.0, u2=-0.05, u3=0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

mdb.models[ModelName].boundaryConditions['BC-7'].deactivate('Step-5')
mdb.models[ModelName].boundaryConditions['BC-8'].deactivate('Step-5')

# shear direction xz----------
faces1 =f1.findAt(((0.001*Nx, 0.001*Ny, Nz),))
region = a.Set(faces=faces1, name='Set-10')
mdb.models[ModelName].DisplacementBC(name='BC-9', createStepName='Step-5', 
    region=region, u1=0.05, u2=0, u3=0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

faces1 = f1.findAt(((0.001*Nx, 0.001*Ny, 0),))
region = a.Set(faces=faces1, name='Set-11')
mdb.models[ModelName].DisplacementBC(name='BC-10', createStepName='Step-5', 
    region=region, u1=-0.05, u2=0, u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

mdb.models[ModelName].boundaryConditions['BC-9'].deactivate('Step-6')
mdb.models[ModelName].boundaryConditions['BC-10'].deactivate('Step-6')
# shear direction xz----------
faces1 =f1.findAt(((0.001*Nx, Ny, .001*Nz),))
region = a.Set(faces=faces1, name='Set-12')
mdb.models[ModelName].DisplacementBC(name='BC-11', createStepName='Step-6', 
    region=region, u1=0.0, u2=0.0, u3=0.05, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)

faces1 = f1.findAt(((0.001*Nx, 0, .001*Nz),))
region = a.Set(faces=faces1, name='Set-13')
mdb.models[ModelName].DisplacementBC(name='BC-12', createStepName='Step-6', 
    region=region, u1=0, u2=0.0, u3=-0.05, ur1=0.0, ur2=0.0, ur3=0.0, 
    amplitude=UNSET, fixed=OFF, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
#--------------------------------------------------------------------------------
#Meshing-------------------------------------------------------------------------
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=STANDARD)
p = mdb.models[ModelName].parts['Part-1']
e = p.edges
#
e1 = e.findAt(((Nx, 0,  .001*Nz),),)
e2 = e.findAt(((0,  Ny, .001*Nz),),)
e3 = e.findAt(((0,  0,  .001*Nz),),)
e4 = e.findAt(((Nx, Ny, .001*Nz),),)

e5 = e.findAt(((Nx, .001*Ny, 0),),)
e6 = e.findAt(((0,  .001*Ny, Nz),),)
e7 = e.findAt(((Nx, .001*Ny, Nz),),)
e8 = e.findAt(((0,  .001*Ny, 0),),)

e5 = e.findAt(((.001*Nx,  0, .0),),)
e6 = e.findAt(((.001*Nx, Ny,  0),),)
#
p.seedEdgeByNumber(edges=(e5+e6), number=Nx/ax, constraint=FINER)
p.seedEdgeByNumber(edges=(e3+e4), number=Ny/ay, constraint=FINER)
p.seedEdgeByNumber(edges=(e1+e2), number=Nz/az, constraint=FINER)
#
p.generateMesh()
c = p.cells
#cells = c.getSequenceFromMask(mask=('[#1 ]', ), )
cells = c.findAt(((.001*Nx, .001*Ny, .001*Nz),),)
p.setElementType(regions=(cells, ), elemTypes=(elemType1, ))
#---------------------------------------------
#Creat two set for two diffrent materials (0 and 1)
setMat0=[]
setMat1=[]
el=p.elements
for i in range (len(C)):
	if C[i]==0:
		setMat0.append(el[i].label)
	else:

		setMat1.append(el[i].label)


#-
elset0=el.sequenceFromLabels(setMat0)
elset1=el.sequenceFromLabels(setMat1)
region1=p.Set(elements=elset0, name='SetMat0')
region2=p.Set(elements=elset1, name='SetMat1')
p.SectionAssignment(region=region1, sectionName='Section-0', offset=0.0, 
	offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
p.SectionAssignment(region=region2, sectionName='Section-1', offset=0.0, 
	offsetType=MIDDLE_SURFACE, offsetField='', thicknessAssignment=FROM_SECTION)
#Job-------------------------------------------------------------------------------

mdb.Job(name=jobName, model=ModelName, description='', type=ANALYSIS, 
    atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
    memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
    explicitPrecision=SINGLE, nodalOutputPrecision=SINGLE, echoPrint=OFF, 
    modelPrint=OFF, contactPrint=OFF, historyPrint=OFF, userSubroutine='', 
    scratch='', multiprocessingMode=DEFAULT, numCpus=CPU_No, numDomains=CPU_No, 
    numGPUs=0)
mdb.jobs[jobName].submit(consistencyChecking=OFF)
mdb.jobs[jobName].waitForCompletion()
#Homogenization------------------------------------------------------------------------

modb = session.openOdb(name=jobName+'.odb')
n1=modb.rootAssembly.instances['PART-1-1']

epMatrix = np.zeros(shape=(6,6), dtype=float)
sigMatrix= np.zeros(shape=(6,6), dtype=float)
Cmatrix  = np.zeros(shape=(6,6), dtype=float)

steps=['Step-1','Step-2','Step-3','Step-4','Step-5','Step-6']

for i in range (6):
	S=modb.steps[steps[i]].frames[-1].fieldOutputs['S'].getSubset(region=n1)
	E=modb.steps[steps[i]].frames[-1].fieldOutputs['E'].getSubset(region=n1)
	for j in range (6):
		for k in range (len(S.values)):
			epMatrix [j][i]= epMatrix [j][i]+E.values[k].data[j]				
			sigMatrix[i][j]= sigMatrix[i][j]+S.values[k].data[j]



ep_inv= np.linalg.inv(epMatrix)
for i in range (6):
	Cmatrix[:,i]=np.dot(ep_inv, sigMatrix[i,:])


time_fin=time.asctime( time.localtime(time.time()) )

fid=open(ModelName+'.txt', 'w')
fid.write("E l a s t i c   H o m o g e n i z a t i o n  U s i n g  A b a q u s\n\n" )
fid.write("-------------------------------------------------\n"  )


vol_phase1=C.sum()/(1.0*Nx*Ny*Nz)

fid.write("Volume fraction Phase 1 = %5.4f\n" % (vol_phase1) )
fid.write("E Phase 1 = %5.4f\n" % (E1) )
fid.write("v Phase 1 = %5.4f\n\n" % (v1) )

fid.write("-------------------------------------------------\n"  )
fid.write("Volume fraction Phase 0 = %5.4f\n" % (1-vol_phase1) )
fid.write("E Phase 0 = %5.4f\n" % (E0) )
fid.write("v Phase 0 = %5.4f\n\n" % (v0) )
fid.write("-------------------------------------------------\n"  )
fid.write("C_homogenized=\n" )
np.savetxt(fid, (Cmatrix[0,:],Cmatrix[1,:],Cmatrix[2,:],Cmatrix[3,:],Cmatrix[4,:],Cmatrix[5,:], ), fmt="%f %f %f %f %f %f")
Smatrix=np.linalg.inv(Cmatrix)

fid.write(" \n\n" )
fid.write("-------------------------------------------------\n"  )
fid.write("E11= %5.4f\n" % (1./Smatrix[0,0]) )
fid.write("E22= %5.4f\n" % (1./Smatrix[1,1]) )
fid.write("E33= %5.4f\n" % (1./Smatrix[2,2]) )


fid.write("-------------------------------------------------\n"  )
fid.write("Start and finish time\n"  )
fid.write(time_strat)
fid.write("\n")
fid.write(time_fin)

fid.close()

dir_name = os.getcwd()
test = os.listdir(dir_name)
for item in test:
    if ( item.endswith(".com") or item.endswith(".inp") or item.endswith(".ipm")or \
		item.endswith(".log") or item.endswith(".prt") or item.endswith(".sim")or \
		item.endswith(".sta") or item.endswith(".rpy") or item.endswith(".1") or item.endswith(".rec")or item.endswith(".C")):
        os.remove(os.path.join(dir_name, item))




























