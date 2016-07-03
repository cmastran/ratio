# -*- coding: utf-8 -*-
"""
Created on Sun Jul 03 00:11:30 2016

@author: cmast
"""

def kill_proc(proc, timeout):
  timeout["value"] = True
  proc.kill()

def process_run(cmd, timeout_sec):
#  proc = subprocess.Popen(shlex.split(cmd),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  proc = subprocess.Popen(shlex.split(cmd),shell=True,stdout=subprocess.PIPE, stderr=subprocess.PIPE)
  timeout = {"value": False}
  timer = Timer(timeout_sec, kill_proc, [proc, timeout])
  try:
      timer.start()
      stdout1, stderr1 = proc.communicate()
  finally:
      timer.cancel()
  return proc, proc.returncode, stdout1, stderr1, timeout["value"]
#  return proc.returncode, stdout.decode("utf-8"), stderr.decode("utf-8"), timeout["value"]



def  plot_pressure(pm,size):
    x =[]
    y = []
    z = []
    plt.figure(figsize=(size,size/1.25))
    # define grid.
    npts = len(pm)
    for jj in range(0,len(pm)):
        x.append(pm[jj][0])
        y.append(pm[jj][1])
        z.append(pm[jj][3])
        
    xi = np.linspace(-17,17,200)
    yi = np.linspace(-17,17,200)
    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,zi,21,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,21,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    plt.title('griddata pressure (%d points)' % npts)
    plt.show()
    return()
    
def  plot_cload(pm,size):
    x =[]
    y = []
    z = []
    plt.figure(figsize=(size,size/1.25))
    # define grid.
    npts = len(pm)
    for jj in range(0,len(pm)):
        x.append(pm[jj][0])
        y.append(pm[jj][1])
        z.append(pm[jj][3])
        
    xi = np.linspace(-17,17,200)
    yi = np.linspace(-17,17,200)
    # grid the data.
    zi = griddata((x, y), z, (xi[None,:], yi[:,None]), method='cubic')
    # contour the gridded data, plotting dots at the randomly spaced data points.
    CS = plt.contour(xi,yi,zi,21,linewidths=0.5,colors='k')
    CS = plt.contourf(xi,yi,zi,21,cmap=plt.cm.jet)
    plt.colorbar() # draw colorbar
    # plot data points.
    plt.scatter(x,y,marker='o',c='b',s=5)
    plt.xlim(-20,20)
    plt.ylim(-20,20)
    plt.title('griddata cload (%d points)' % npts)
    plt.show()
    return()  
    
#    function to make geo file
#
def make_geo_file(filename):
   f = open(filename,"w")
   f.write("Merge \"rplate.stl\";\n")
   f.write("Surface Loop (1)= {1};\n")
   f.write("Volume (1)= {1};\n")
   f.write("Recombine Surface {1}; // recombine triangs into quads\n")
   f.write("Coherence Mesh;\n")
   f.write("Mesh.Algorithm = 1;\n") 
   f.write("Mesh.RemeshAlgorithm = 1;\n // automatic\n")
   f.write("Mesh 3;\n")
   f.write("OptimizeMesh \"Netgen\";\n")
   # comment this remeshing one
#   f.write("RefineMesh;\n")
   #
   f.write("Mesh.Smoothing = 200;\n")
   f.write("OptimizeMesh \"Netgen\";\n")
   f.write("Mesh.Smoothing = 200;\n")
   f.write("OptimizeMesh \"Netgen\";\n")

   f.write("Mesh.Format = 1;\n")
   f.write("Mesh.SaveAll = 1;\n")
 
#  f.write("Mesh.Smoothing = 100;\n")
   f.write("Save \"mesh3dproc.msh\";\n")
   f.close()
   return
#
   
def   fnodewriter(f,node_list, spacer):
    n = len(node_list)
    nc = 0
    f.write("\t")
    for i in range(0,n):
        nc = nc + 1  # add 1 as nodes in getfem start from 0, not 1
        node = node_list[i]+1
        f.write(str(node))
        if (i < n-1 ):
            f.write(",\t")
        if ( nc == spacer ):
            f.write("\n\t")
            nc = 0
    f.write("\n")        
    return()
    
def   fcloadlabelwriter(f,bottom):
    for i in range(0,len(bottom)):
        f.write("*NSET, NSET=LOAD"+str(bottom[i]+1)+"\n")
        f.write(str(bottom[i]+1)+"\n")
    return()
    
def   fcloadwriter(f,bottom,cload):
    for i in range(0,len(bottom)):
        f.write("LOAD"+str(bottom[i]+1)+", 3, "+str(cload[i][2])+"\n")
    return()
   
def   assemble_input_deck(numel,material,E,nu,deckfile,meshfile,side_pts,bottom_pts,bottom_pts_cload):
    f = open(deckfile,"w")
#    f.write("*HEADING\n")
#    f.write("Model: diaphragm     Date: whoknows\n")
    f.write("*INCLUDE, INPUT="+meshfile+"\n")
    f.write("*NSET, NSET=FIX"+"\n")
    fnodewriter(f,side_pts.tolist(),10)
    f.write("*BOUNDARY \nFIX, 1"+"\n")
    f.write("*BOUNDARY \nFIX, 2"+"\n")
    f.write("*BOUNDARY \nFIX, 3"+"\n")
    f.write("*NSET, NSET=Nall, GENERATE"+"\n")
    f.write("1,"+str(numel)+"\n")
    f.write("*MATERIAL, NAME="+material+"\n")
    f.write("*ELASTIC \n")
    f.write(str(E)+",\t"+str(nu)+"\n")
    f.write("*SOLID SECTION, ELSET=Volume0, MATERIAL="+material+"\n")
    fcloadlabelwriter(f,bottom_pts)
    f.write("**\n")
    f.write("*STEP\n")
    f.write("*STATIC\n")
    f.write("*CLOAD\n")
    fcloadwriter(f,bottom_pts,bottom_pts_cload)
    f.write("*NODE PRINT, NSET=Nall\n")
    f.write("U\n")
    f.write("*EL PRINT, ELSET=Volume0\n")
    f.write("S\n")
    f.write("*NODE FILE\n")
    f.write("U\n")
    f.write("*EL FILE\n")
    f.write("S\n")
    f.write("*END STEP\n")
    f.close()
    return()
   
   
def   add_segment_number_to_oscad(file,nseg):
   ii = 0
   line_list = []
   f = open(file)
   line_list.append(f.readline())
   while len(line_list[ii])!=0:
        pp = "ii="+str(ii)+":"+line_list[ii]
##        sys.stdout.write(pp)
        line_list.append(f.readline())
        ii = ii+1
   f.close()
   num_lines = ii 
   # now rewrite file with first line fixed 
   f = open(file,"w")
   f.write("$fn="+str(nseg)+";\n")
   for ii in range(0,len(line_list)):
       f.write(line_list[ii])
   f.close()
   return
#  
#       Python function to strip 2D Triangular elements from mesh 
#
def strip_2D_elements_from_mesh(infile,outfile):
   #
   #    read the input mesh file
   #
   ii = 0
   line_list = []
   f = open(infile)
   line_list.append(f.readline())
   while len(line_list[ii])!=0:
#        pp = "ii="+str(ii)+":"+line_list[ii]
##        sys.stdout.write(pp)
        line_list.append(f.readline())
        ii = ii+1
   f.close()
   num_lines = ii # this is the end of the file
   #
   #  now scan the lines and read nodes
   #
   nodeline = line_list.index("$Nodes\n")
   endnodeline = line_list.index("$EndNodes\n")
   numnodes = int(line_list[nodeline+1])
   for ii in range(nodeline+2,endnodeline):
       ml =  line_list[ii].split(" ")
       stlin = "%d %1.9f %1.9f %1.9f\n" % (int(ml[0]), float(ml[1]),  float(ml[2]),  float(ml[3]))
       line_list[ii] = stlin
        
   elementline = line_list.index("$Elements\n")
   endelementline = line_list.index("$EndElements\n")
   numelements = int(line_list[elementline+1])
   #
   #  now scan the lines to see beginning of elements
   #
   elementline = line_list.index("$Elements\n")
   endelementline = line_list.index("$EndElements\n")
   numelements = int(line_list[elementline+1])
   
##   print "Elements = ", elementline
##   print "Elements = ", line_list[elementline]
##   print "Elements = ", line_list[elementline+1]
   #
   #   here count the number of 3D elements
   solid_count = 0
   elementcountline = elementline+1
   firstelementline = elementline+2
 
   for ii in range(firstelementline,endelementline):
      mylist = line_list[ii].split(" ")
#      mylist[0] = str(int(mylist[0])-firstel_num + 1)  # reset first element num to 1
      if (mylist[1] == "4" ):
         solid_count = solid_count+1
#      linestr = mylist[0] + " " + mylist[1] + " " + mylist[2] + " "
#      linestr = linestr + mylist[3] + " " + mylist[4] + " " + mylist[5] + " "
#      linestr = linestr + mylist[6] + " " + mylist[7] + " " + mylist[8]
#      line_list[ii] = linestr
##   print "mesh has ",solid_count, "tetrahedrons"
   sys.stdout.write("("+str(solid_count)+" tetrahedrons) ... ")
##   print "triangles =", numelements-solid_count
   #
   #  now print the file without triangles
   #
   fo = open(outfile,"w")
   for ii in range(0,elementcountline):
      fo.write(line_list[ii])
   #
   #  now write the number of elements
   #
   fo.write(str(solid_count)+"\n")
   #
   #  now print solids with renamed element numbers
   #
   first_tet_found = False
   count = 0
   for ii in range(firstelementline,endelementline):
      mylist = line_list[ii].split(" ")
      if (mylist[1] == "4" ):
          count = count + 1
          if (first_tet_found == False ):
              firstet =  line_list[ii].split(" ")
              #print "firstet = ", firstet
              firstet_num = int(firstet[0])  # get first number 
              #print "firstet_num = ", firstet_num
              first_tet_found = True
              
          mylist = line_list[ii].split(" ")
          long_str = str(count) + " " + mylist[1] + " " + mylist[2] + " "
          long_str = long_str + mylist[3] + " " + mylist[4] + " " + mylist[5] + " "
          long_str = long_str + mylist[6] + " " + mylist[7] + " " + mylist[8]
          line_list[ii] = long_str
          fo.write(line_list[ii]) 
   #
   #   write last line
   #
   fo.write("$EndElements")
   fo.close()
   return
   
def clean_inp_file(infile,outfile):
   #
   #    read the input inp file
   #
   ii = 0
   line_list = []
   f = open(infile)
   line_list.append(f.readline())
   while len(line_list[ii])!=0 :
        line_list.append(f.readline())
        ii = ii+1
   f.close()
   num_lines = ii # this is the end of the file
   #
   #  now scan the lines and read nodes
   #
   nodeline = line_list.index("*NODE\n")
   endnodeline = line_list.index("******* E L E M E N T S *************\n")
   numnodes = endnodeline - nodeline - 1
   for ii in range(nodeline+1,endnodeline):
       ml =  line_list[ii]
       ml = ml.replace(',',' ')
    #   print "ml = ", ml
       ml2 =  ml.split()
   #    print "ml2 = ", ml2
       ml =  ml.split()
       stlin = "%d, %1.9f, %1.9f, %1.9f\n" % (int(ml[0]), float(ml[1]),  float(ml[2]),  float(ml[3]))
       line_list[ii] = stlin
        
   #
   #  now print the nodes with floats
   #
   fo = open(outfile,"w")
   for ii in range(0,num_lines):
      fo.write(line_list[ii])
   #
   fo.close()
   return 
   
   
def calculix_extreme_z(infile):
   #
   #    read the result calculix data file
   #
   ii = 0
   line_list = []
   nn = []
   dx = []
   dy = []
   dz = []
   f = open(infile)
   line_list.append(f.readline())
   while len(line_list[ii])!=0 :
        line_list.append(f.readline())
        ii = ii+1
   f.close()
   num_lines = ii # this is the end of the file
   #
   #  now scan the lines and read nodes
   #
   
   sub = "displacements (vx,vy,vz)"
   for ii in range(0,len(line_list)):
       text = line_list[ii]
       if sub in text:
           nodeline = ii
    
      
   sub = "stresses"
   for ii in range(0,len(line_list)):
       text = line_list[ii]
       if sub in text:
           endnodeline = ii
   
   numnodes = endnodeline - nodeline - 3
   
   j = 0
   for ii in range(nodeline+2,endnodeline-1):
       ml =  line_list[ii]
   #    ml = ml.replace(',',' ')
    #   print "ml = ", ml
    #   ml2 =  ml.split()
   #    print "ml2 = ", ml2
       ml =  ml.split()
       nn.append(int(ml[0]))
       dx.append(float(ml[1]))
       dy.append(float(ml[2]))
       dz.append(float(ml[3]))
       j = j+1
   #
   #  now find the extremes
   #
   count = j
   mindz = dz[0]
   maxdz = dz[0]
 
   for ii in range(0,len(nn)):
       if ( mindz > dz[ii]):
           mindz = dz[ii]
       if ( maxdz < dz[ii]):
           maxdz = dz[ii]
   #
   return(mindz,maxdz)    
   
   
   
def find_solid_centroids(m1):
    ct = []
    for jj in range(0,m1.nbcvs()):
        # get points from the mesh solid tetrahedrons
        (pts,ID1) = m1.pts_from_cvid(jj)
        # calculate the centroid of the four points
        ct1 = [sum(pts[0,:])/4.0,sum(pts[1,:])/4.0 ,sum(pts[2,:])/4.0]
        # append the centroid coordinates to the list
        ct.append(ct1)
    return(ct) 
    

    
def force_at_bottom_points(m2,bottom_pts_ids,bottom_pts_cload):
    fmsh = []
    pts = m2.pts(bottom_pts_ids)
    ptsl =  m2.pts().tolist()
    for i in range(0,len(bottom_pts_ids)):
        num = bottom_pts_ids[i]
#        print "num = ", num
#        print "len pts = ",len(ptsl[0])
        p1 = m2.pts(num)
        x = ptsl[0][num]
        y = ptsl[1][num]
        z = ptsl[2][num]
        dat = bottom_pts_cload[i]
 #       print "dat = ", dat
        dat1 = dat[2]
 #       print dat1
 #       dat1 = x
        fm = [ x, y, z, dat1 ]

        fmsh.append(fm)
    return(fmsh)
        
    return(fmsh)
#unit normal vector of plane defined by points a, b, and c
def unit_normal(a, b, c):
    x = np.linalg.det([[1,a[1],a[2]],
         [1,b[1],b[2]],
         [1,c[1],c[2]]])
    y = np.linalg.det([[a[0],1,a[2]],
         [b[0],1,b[2]],
         [c[0],1,c[2]]])
    z = np.linalg.det([[a[0],a[1],1],
         [b[0],b[1],1],
         [c[0],c[1],1]])
    magnitude = (x**2 + y**2 + z**2)**.5
    return (x/magnitude, y/magnitude, z/magnitude)

#area of polygon poly
def poly_area(poly):
    if len(poly) < 3: # not a plane - no area
        return 0
    total = [0, 0, 0]
    N = len(poly)
    for i in range(N):
        vi1 = poly[i]
        vi2 = poly[(i+1) % N]
        prod = np.cross(vi1, vi2)
        total[0] += prod[0]
        total[1] += prod[1]
        total[2] += prod[2]
    result = np.dot(total, unit_normal(poly[0], poly[1], poly[2]))
    return abs(result/2.0)   

def find_solid_bottom_facet_cloads(m3,bottom,centroid_force):
   
    # find the IDs of all points at bottom
    bottom_pts_ida= m3.pid_in_faces(bottom)
    bottom_pts_ids= m3.pid_in_faces(bottom).tolist()
    cloads = []
    pids = []
    pts = m3.pts(bottom_pts_ida) # it does need argument equal to pid list
    num_bot_pts = len(bottom_pts_ids)
    # initialize loads to [0.0, 0.0, 0.0]
    for index in range(0,num_bot_pts):
        pids.append(bottom_pts_ids[index])
        cloads.append([0.0, 0.0, 0.0])
        
#    return(pids,cloads)
    #for all faces
    for face_index in range(0,len(bottom[0])):
        # here find the point IDs
        element_point_ids = m3.pid_in_faces(bottom[:,face_index])
#        print "face ids =" , bottom[:,face_index]
#        print "element point ids =" , element_point_ids
#        print
         # set force on corresponding node for all 3 nodes 
        for j in range(0,3):
            # set force on corresponding node for all 3 nodes 
            elind = pids.index(element_point_ids[j])
            cload_node = cloads[elind]
#            print "for point =", element_point_ids[j]
#            print "orig cload_node =", cload_node
            # only z-component matters
            cload_node[2] = cload_node[2] + 1/3.0*centroid_force[face_index]
#            cload_node[2] = 1/3.0*centroid_force[face_index]
            coord = pts[:,elind]
            x = coord[0]
 #           cload_node[2] = x
#            print "coord = ",coord
#            print "fin cload_node =", cload_node
            cloads[elind] = cload_node
#            print "cload_node =", cload_node
 #           wait = raw_input("PRESS ENTER TO CONTINUE.")
#        print
    
   # sys.exit("Stopping here")
    return(pids,cloads)
  
    
def find_solid_bottom_facet_centroids1(m2,fbot,hb,epsilon):
    ct = []
    ta = []
    fac_pts = []
    pc = []
    tot_pt_IDs = m2.pid_in_faces(fbot)
    p = m2.pts(tot_pt_IDs)
#    print "pts = ", p
#    print  "points in bottom = ", len(p[0])
    for jj in range(0,len(fbot[0])):
        # get points from the mesh solid tetrahedrons
        ftlist =  fbot.tolist()[0]
        fnlist =  fbot.tolist()[1]
        facet = [ftlist[jj],fnlist[jj]]
        facet_arr = array([facet])
        face = facet_arr.transpose()
        pt_IDs = m2.pid_in_faces(face)  #here are the points
        
#        print "pt_IDs = ", pt_IDs
        #
        # calculate the centroid of the four points
        #  check that at leas 3 points are on hb plane
        #
        pt_IDs_list = pt_IDs.tolist()
        np = 0
        sx = 0.0
        sy = 0.0
        sz = 0.0
        pid = []
        pp = []
        for kk in range(0,3):  # for all 3 points in face
 #           print "select point", pt_IDs_list[kk]
            num = pt_IDs_list[kk]
 #           pl = p[:,num]
            pl0 = m2.pts(num)
            pl = [ float(pl0[0]), float(pl0[1]), float(pl0[2]) ]
           
#            print "pl = ",pl
            if ( abs(pl[2] -hb) < epsilon*abs(hb) ):
                sx = sx + pl[0]
                sy = sy + pl[1]
                sz = sz + pl[2]
                
                pt1 = [pl[0],pl[1],pl[2]]
                pp.append(pt1)
                pid.append(pt_IDs_list[kk])
                np = np + 1
                
        sx = sx/3.0
        sy = sy/3.0
        sz = sz/3.0
        
        # here print centroid and points
#        print "face",face
#        print "pts = ",pt_IDs_list 
        for kk in range(0,3):  # for all 3 points in face
            num = pt_IDs_list[kk]
            pl0 = m2.pts(num)
            pl = [ float(pl0[0]), float(pl0[1]), float(pl0[2]) ]
#            print "pt: ", pt_IDs_list[kk], pl
#        print "av: ","[",sx, " ", sy, " ",sz,"]"
 #       wait = raw_input("PRESS ENTER TO CONTINUE.")
        #
        #     find the convex and facet number
        #
        if ( np >= 3 ): # then it is a plane facet
        # insert points which are not already in the list
            for ii  in range (0,2):
                num =  pt_IDs_list[ii]
                try:
                    fac_pts.index(num)
                except ValueError:
                    # if not found, then add it
 #                   pl = p[:,npt]
                    pl = m2.pts(num)
 #                   pl1 =  sum(pl, [])
#                    print "npt=",npt,"pl =", pl
#                    wait = raw_input("PRESS ENTER TO CONTINUE.")
                    fac_pts.append(num)
                    x = float(pl[0])
                    y = float(pl[1])
                    z = float(pl[2])
                    vc = [ x, y, z ]
#                    print "vc = ", vc
                    pc.append(vc)
         
            ct1 = [sx, sy ,sz]  # here add avg
            ta.append(poly_area(pp))
            ct.append(ct1)
                
                # now sort points according to their x-coordinate
 #   wait = raw_input("PRESS ENTER TO CONTINUE.")    
    sort_done = False
    while ( sort_done == False ):
        sort_done = True
        for ii in range(0,len(fac_pts)-1):
            p1 = pc[ii]
            p2 = pc[ii+1]
            x1 = p1[0]
            x2 = p2[0]
            if  ( x1 > x2 ):
                # swap them
                sort_done = False
                temp_id = fac_pts[ii]
                fac_pts[ii] = fac_pts[ii+1]
                fac_pts[ii+1] = temp_id
                temp_vc = pc[ii]
                pc[ii] = pc[ii+1]
                pc[ii+1] = temp_vc
       # done sorting 
#    print "done sorting "
    return(ct,ta,fac_pts,pc)
      
        
def sum_of_centroids_on_bot(centr,m2,cbot):
    av = [0.0, 0.0, 0.0 ]
    kk = 0
    for jj in range(0,len(centr)):
        # get points from the mesh solid tetrahedrons
      if ( centr[jj][0]*centr[jj][0]+centr[jj][1]*centr[jj][1]+centr[jj][2]*centr[jj][2] > 1.0e-5):
          kk = kk + 1
          av[0] = av[0] + centr[jj][0]
          av[1] = av[1] + centr[jj][1]
          av[2] = av[2] + centr[jj][2]
     
    av[0] = av[0]/kk
    av[1] = av[1]/kk
    av[2] = av[2]/kk 
    sr = "["+str(av[0])+","+str(av[1])+","+str(av[2])+"]"
    return(sr)
    
def sum_of_xy_on_bot(pp,cb):
    av = [0.0, 0.0, 0.0]
    cnt = 0
    print "len pp = ", len(pp[0])
    for jj in range(0,len(pp[0])):
        if (cb[jj] == True):
            cnt = cnt +1
    print "point at bot count = ",cnt
    for jj in range(0,len(pp[0])):
        # get points from the mesh solid tetrahedrons
        if ( cb[jj] == True ):
            av[0] = av[0] + pp[0][jj]
            av[1] = av[1] + pp[1][jj]
            av[2] = av[2] + pp[2][jj]
            
    sr = "["+str(av[0])+","+str(av[1])+","+str(av[2])+"]"
    return(sr)
    
def pressure_at_centroids(ct1,pu,pcoeff):
    pc = []
    pm = []
    for jj in range(0,len(ct1)):
  #      pr = 0.5*(ct1[jj][0]/15.0)
        # calculate the face presure at the mid point
        pr = pcoeff*ct1[jj][0] + pu
#
        x = ct1[jj][0]
        y = ct1[jj][1]
      
        pm1 = [ct1[jj][0],ct1[jj][1],ct1[jj][2],pr]
       
        pc.append(pr)
        pm.append(pm1)
    return(pc,pm) 
    
def force_at_centroids(pc,bfc,af):
    fc = []
    for jj in range(0,len(bfc)):    
 #       fc1 = pc[jj]*af[jj]
        fc1 = pc[jj]*af[jj] 
#        fc1 = pc[jj]*af[0] 
        fc.append(fc1)
    return(fc)     
    
def find_pressure(ffbot,bot_facets,press_centroid):
    
    for i in range(0,len(press_centroid)):
        if (ffbot[0] == bot_facets[i][0] and ffbot[1] == bot_facets[i][1]):
            return(press_centroid[i])         
    return (-1.0)

def find_force(ffbot,bot_facets,force_centroid):
    
    for i in range(0,len(press_centroid)):
        if (ffbot[0] == bot_facets[i][0] and ffbot[1] == bot_facets[i][1]):
            return(force_centroid[i])
            
    return (-1.0)   
    
#  
def tensors_from_pressure(pr1):
    pr2 = []
    for kk in range(0,len(pr1)):
        # put uniform pressure at the centroid of the face 
        pr3 = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,pr1[kk]]]
     #   pr3 = [[0.0,0.0,0.0],[0.0,0.0,0.0],[pr1[kk],0.0,0.0]]
        # the line below is for uniform loading
        # pr3 = [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,1.0]]
        pr2.append(pr3)
    pr4 = np.transpose(array(pr2))
    return(pr4)       
  
def pointnotinlist(xpt,ypt,zpt,xpa,ypa,zpa,dip):
    dt = 0.0
    for gg in range(0,len(xpa)):
        dt = (xpt-xpa[gg])*(xpt-xpa[gg])+(ypt-ypa[gg])*(ypt-ypa[gg])
        if (math.sqrt(dt) < dip):
            return(False) # point in list already
    return(True)
    
def backplate_borderpoints(mesh,zplate):
    xp = []
    yp = []
    zp = []
    total_points = mesh.points
    for mm in range(0,len(mesh)):
        zf = True
        for oo in range(0,3):   # cycle for all points in triangle
            zp1 = total_points[mm][2+3*oo]
            zf = zf and ( math.fabs(zp1-zplate) < 0.01 )
    # na=make sure all points are on the plane
        if ( zf == True ):
    # if point on back plate then add it to list
            for oo in range(0,3):
                xp1 = total_points[mm][0+3*oo]
                yp1 = total_points[mm][1+3*oo]
                zp1 = total_points[mm][2+3*oo]
                if ( pointnotinlist(xp1,yp1,zp1,xp,yp,zp,0.01) == True):
            # attach it if not repeated
                    xp.append(total_points[mm][0+3*oo])
                    yp.append(total_points[mm][1+3*oo])
                    zp.append(total_points[mm][2+3*oo])
    return(xp,yp,zp)
    
#
#      read slt file 
#
def ReadSTLTriangleFile(infile):
   #
   #    read the input STL file
   #
   ii = 0
   line_list = []
   triangles = []
   f = open(infile)
   line_list.append(f.readline())
   while len(line_list[ii])!=0 :
#        pp = "ii="+str(ii)+":"+line_list[ii]
##        sys.stdout.write(pp)
        line_list.append(f.readline())
        ii = ii+1
   f.close()
   num_lines = ii # this is the end of the file
   #
   # now parse list and extract triangles
   #
   last_line = line_list[num_lines-1] # This is the first line
   name = line_list[0]  # The first line is the solid name
   triangles = ReadTriangles(line_list)
   return(name,last_line,triangles)
   
#
#   write stl file
# 
def WriteSTLTriangleFile(tfirst,tlast,tri,outfile):
    f = open(outfile,"w")
    f.write(tfirst)
    for ii in range(0,len(tri)):
        tri1 = tri[ii]
        for jj in range(0,len(tri1)):
            f.write(tri1[jj])
    f.write(tlast)
    f.close
#
#  read individual triangles
#
def ReadTriangles(lines):
    tr = []
    numl = len(lines)
 #   print "numl = ", numl
    
    ii = 1
    while ii < numl-2:
        pr = []
      #  print "ii =", ii, "+",lines[ii]
        l1 = lines[ii+0]  # facet statement
        l2 = lines[ii+1] # loop statement
        v1 = lines[ii+2] # first vertex
        v2 = lines[ii+3] # second vertex
        v3 = lines[ii+4] # third vertex
        l3 = lines[ii+5] # endloop statement
        l4 = lines[ii+6] # end facet statement
        # now assemble triangle array
        pr = [l1,l2,v1,v2,v3,l3,l4]
        tr.append(pr)
        ii = ii+7        # advance counter
    return(tr)
#
#     find triangles which are coplanar with z = zp within dist1
#
def TriangleAtHeight(tri,zplate,dist1):
    #
    #     test if triangle contined at z= zp
    # 
    zf = True
    for oo in range(0,3):   # cycle for all points in triangle
    #   print "tri =", tri
        zps = tri[2+oo]
    #   print "zps[",oo,"] =", zps
        zpc = zps.lstrip()
        zp1 = re.split(" +",zpc)
    #   print "zp1 = ", zp1
        zpn = float(zp1[3])
    #   print "zpn = ", zpn
        zf = zf and ( math.fabs(zpn-zplate) < dist1 )
    if zf == True :
        return(zf)
    return(False)

# 
#   remove triangles coplanar at z = zp and return new list
#
def RemoveTrianglesAtHeight(trlist,zp,dist1):
    trlist1 = []
    for ii in range(0,len(trlist)):
        if TriangleAtHeight(trlist[ii],zp,dist1) == False:
            trlist1.append(trlist[ii])
    return(trlist1)
    
    
def inlist(a,b,index):
    #print "a[index]=", a[index]
    if a[index] in b:
        return(True)
    else:
        return(False)
