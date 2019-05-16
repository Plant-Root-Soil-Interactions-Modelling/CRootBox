import py_rootbox as rb
from rb_tools import *
import matplotlib.pyplot as plt

rootsystem = rb.RootSystem()

# Open plant and root parameter from a file
#name = "Zea_mays_1_Leitner_2010" # "Anagallis_femina_Leitner_2010"
#name = "Crypsis_aculeata_Clausnitzer_1994" # "Anagallis_femina_Leitner_2010"
#name = "Juncus_squarrosus_Clausnitzer_1994"
#name = Brassica_oleracea_Vansteenkiste_2014 #Brassica_napus_a_Leitner_2010
name = "Mai2019/Rice_NERICA4_NoP_Drying_Mai_2019V5_1"
outputname = name
ColumnHeight = 55
ColumnRadius = 8
IsBox = False
simtime = 51
#massPerLength = 1/59/100
massPerVolume = 0.15 # g cm-3
realData = [43.75, 25.00, 31.25]
realMass = [1.6]
realVolume = realMass[0]/massPerVolume
#realMass = realVolume[0]*massPerVolume #/5
#realSurface = [2534.11533156*1.8]
allRS = [ ]
Column = rb.SDF_PlantContainer(ColumnRadius, ColumnRadius, ColumnHeight, IsBox)
###########
#
# Creates N root systems
#
N=1
for i in range(0,N):
     rs = rb.RootSystem()
     rs.openFile(name)
     #rparams = rs.getRootSystemParameter()
     #rparams.delayB = 1
     #rparams.interBranchDistanceScaleFactor = 5.5
     #rparams.interBranchDistanceTransitionDepth = -18
     #rparams.laterTypeTransitionDepth = -25
     p1 = rs.getRootTypeParameter(1)
     #p1.r = 1.5
     #p1.lb = 1000
     #p1.theta = 1
     #p1.la = 0.5
     #p2 = rs.getRootTypeParameter(2)
     #p2.lmax = 2
     #p5 = rs.getRootTypeParameter(5)
     #p5.lmax = 4
     #p1.r = 1
     rs.setGeometry(Column)
     allRS.append(rs)

#
# Simulate
#
for rs in allRS:
    rs.initialize()
    #for i in range(50):
    # rs.simulate(1)
    # rs.simulate(25)
    rs.simulate(simtime)
#    rs.simulate(1)

##
## Export results as single vtp files
##
c = 0
for rs in allRS:
      c += 1 # root system number
      #vtpname = "results/"+outputname+str(c)+".vtp"
      #rs.write(vtpname)
      ana = rb.SegmentAnalyser(rs)
      ana.filter(3, 0, simtime) #creationTimeId = 3 #st_time
      ana.write("results/"+outputname+str(c)+"_.dgf")
      ana.write("results/"+outputname+str(c)+"_.vtp")

#
# Compute vertical RLD distribution in layers
#
nl = 55; # number of layers
vRLD=np.zeros((N,nl)); # N rows, nl columns
vRootSurfacePerLayer=np.zeros((N,nl)); # N rows, nl columns
vRootLengthPerLayer=np.zeros((N,nl)); # N rows, nl columns
vRL=np.zeros((N,3));
vRootSurface=np.zeros((N,3));
vRootLength=np.zeros((N,3));
TotalRootVolume = np.zeros((N,1));
TotalRootSurface = np.zeros((N,1));
TotalRootLength = np.zeros((N,1));
depth = ColumnHeight
c=0
for rs in allRS:
      analysis = rb.SegmentAnalyser(rs)
      RootLengthPerLayer = analysis.distribution(rb.ScalarType.length,0,depth,nl,True)
      RLD = analysis.distribution(rb.ScalarType.volume,0,depth,nl,True)
      RootSurfacePerLayer = analysis.distribution(rb.ScalarType.surface,0,depth,nl,True)
      vRLD[c,:] = RLD
      vRootSurfacePerLayer[c,:] = RootSurfacePerLayer
      vRootLengthPerLayer[c,:] = RootLengthPerLayer
      #vRLD[c,:] /= (depth/nl)
      for i in range(0,15):
            vRL[c,0] += vRLD[c,i]
            vRootSurface[c,0] += vRootSurfacePerLayer[c,i]
            vRootLength[c,0] += vRootLengthPerLayer[c,i]
            TotalRootVolume[c] += vRLD[c,i]
            TotalRootSurface[c] += vRootSurfacePerLayer[c,i]
            TotalRootLength[c] += vRootLengthPerLayer[c,i]
      for i in range(15,30):
            vRL[c,1] += vRLD[c,i]
            vRootSurface[c,1] += vRootSurfacePerLayer[c,i]
            vRootLength[c,1] += vRootLengthPerLayer[c,i]
            TotalRootVolume[c] += vRLD[c,i]
            TotalRootSurface[c] += vRootSurfacePerLayer[c,i]
            TotalRootLength[c] += vRootLengthPerLayer[c,i]
      for i in range(30,55):
            vRL[c,2] += vRLD[c,i]
            vRootSurface[c,2] += vRootSurfacePerLayer[c,i]
            vRootLength[c,2] += vRootLengthPerLayer[c,i]
            TotalRootVolume[c] += vRLD[c,i]
            TotalRootSurface[c] += vRootSurfacePerLayer[c,i]
            TotalRootLength[c] += vRootLengthPerLayer[c,i]
      c += 1 # root system number

#
# Ploting
#
z=[-15,-30,-55]
total = np.sum(vRL,axis=0)
mean=np.mean(vRL,axis=0)
total = mean[0]+mean[1]+mean[2];
#std=np.std(vRL,axis=0)
#plt.plot(mean*massPerLength,z,'k-', color="blue",linewidth=2)
##plt.plot(mean/total,z,'k-', color="blue",linewidth=2)
#plt.fill_betweenx(z,(mean+std)*massPerLength,(mean-std)*massPerLength,color="blue",edgecolor='',alpha=0.5)
#x1,x2,y1,y2 = plt.axis()
#plt.axis((0,x2,y1,y2)) # set min of x axis to 0, because of some negative values of (mean-std)
#plt.xlabel('g')
#plt.ylabel('Depth (cm)')
#plt.show()
print ("Simulated Volume[cm3] - Data volume[cm3]")
print (TotalRootVolume[0], realVolume)
print ("Simulated Mass[g] - Data mass[g]")
print (TotalRootVolume[0]*massPerVolume, realMass)
#print ("Simulated Surface - Data surface")
#print (TotalRootSurface[0], realSurface)
print ("Simulated Volume in each soil segments")
print (vRL[0][0]/TotalRootVolume[0]*100, vRL[0][1]/TotalRootVolume[0]*100, vRL[0][2]/TotalRootVolume[0]*100)
print ("Simulated Surface in each soil segments")
print (vRootSurface[0][0]/TotalRootSurface[0]*100, vRootSurface[0][1]/TotalRootSurface[0]*100, vRootSurface[0][2]/TotalRootSurface[0]*100, TotalRootSurface[0])
print ("Simulated Length in each soil segments")
print (vRootLength[0][0]/TotalRootLength[0]*100, vRootLength[0][1]/TotalRootLength[0]*100, vRootLength[0][2]/TotalRootLength[0]*100, TotalRootLength[0])

basalZone = plt.scatter(realData[0],vRL[0][0]/TotalRootVolume[0]*100, color="blue", marker='x')
intermediateZone = plt.scatter(realData[1],vRL[0][1]/TotalRootVolume[0]*100, color="red", marker='o')
deepZone = plt.scatter(realData[2],vRL[0][2]/TotalRootVolume[0]*100, color="green", marker='<')

refLine = np.arange(np.min(realData), np.max(realData),(np.max(realData) - np.min(realData))/5)
plt.plot([np.min(realData), np.max(realData)], [np.min(realData), np.max(realData)], 'b--')
plt.legend((basalZone, intermediateZone, deepZone),
           ('basalZone', 'intermediateZone', 'deepZone'))
plt.xlabel('Root fraction from experiment data (%)')
plt.ylabel('Root fraction from from CRootBox (%)')
#x1,x2,y1,y2 = plt.axis()
#plt.axis((x1-10,x2+10,x1-10,x2+10)) # set min of x axis to 0, because of some negative values of (mean-std)
plt.show()

#
# Ploting
##
#nameOutputFile = name+"_RVD"
#z=np.linspace(0,depth*(-1),nl)   # depth*-1 is the (negativ) z coordinate
#mean=np.mean(vRLD,axis=0)
#std=np.std(vRLD,axis=0)
##plt.figure(figsize=(3.8,3))
#plt.plot(mean,z,'k-', color="blue",linewidth=2)
#plt.fill_betweenx(z,mean+std,mean-std,color="blue",edgecolor='',alpha=0.5)
#x1,x2,y1,y2 = plt.axis()
#plt.axis((0,x2,y1,y2)) # set min of x axis to 0, because of some negative values of (mean-std)
#plt.xlabel('RVD (cm3/cm)')
#plt.ylabel('Depth (cm)')
##plt.savefig(nameOutputFile+".png")
#plt.show()
## save data: z, mean, std in npz format (for multi array)
##https://docs.scipy.org/doc/numpy/reference/generated/numpy.savez.html
##np.savez(nameOutputFile, z=z, mean=mean, std=std)

