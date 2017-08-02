import scipy.sparse.linalg as LA
from scipy import sparse

import py_rootbox as rb    
from rb_tools import *

import xylem_flux 

# Simulate a root system
name = "Anagallis_femina_Leitner_2010" #"Sorghum_bicolor_NA_NA" "Zea_mays_4_Leitner_2014" # 
rs = rb.RootSystem()
rs.openFile(name)
rs.initialize() 

# example
simtime = 60 # days
dt = 1
N = round(simtime/dt)


# Incrementally build nodes and segments
nodes = vv2a(rs.getNodes()) # contains the initial nodes of tap, basal and shootborne roots
seg = np.array([], dtype=np.int64).reshape(0,2)

print(nodes)

print()
for i in range(0,N):        
  
     rs.simulate(dt)           
     print("Number of nodes", rs.getNumberOfNodes()) # equals the number of new segments
     print("Number of roots", len(rs.getRootTips()))
            
     print("Number of new nodes", rs.getNumberOfNewNodes()) # equals the number of new segments
     print("Number of new roots", rs.getNumberOfNewRoots())

     uni = v2ai(rs.getUpdatedNodeIndices())    
     unodes =  vv2a(rs.getUpdatedNodes())
     print("Number of node updates", len(unodes), len(uni))
     # do the update
          
     nodes[uni] = unodes     
             
     newnodes = vv2a(rs.getNewNodes())
     newsegs = seg2a(rs.getNewSegments())
     
     if rs.getNumberOfNewNodes()!=newnodes.shape[0]:
         print("oh noooooooo") 
     
#      print("New nodes: ")
#      print(newnodes,"\n")
#      print("New segments: ")
#      print(newsegs)
     if len(newnodes)!=0:
         nodes = np.vstack((nodes,newnodes))
     if len(newsegs)!=0:
         seg = np.vstack((seg,newsegs))
     
    # test is everything right?
     nodes_ = vv2a(rs.getNodes());
     seg_ = seg2a(rs.getSegments());
    
     uneq = np.sum(nodes_!=nodes)/3
    
     if uneq>0:    
#         print("Incremental nodes", len(nodes), " real", len(nodes_))
#         print("Incremental segments", len(seg), " real", len(seg_))   
        print("unequal nodes: ", uneq) 

        i = np.nonzero(nodes_[:,0]!=nodes[:,0])
        print(i)
        print()

        print(nodes[i,:])
        print(nodes_[i,:])     
         
     
     
     print()



# node indices have meaning, the ordering should be the same




# # segment indices have no special meaning, and the ordering is different
# seg = np.sort(seg,axis=0) # per default along the last axis 
# seg_ = np.sort(seg_,axis=0) 
# print("unequal segs: ", np.sum(seg_!=seg)/2) 
# 
# 
# rs.write("results/example_5b.vtp")

