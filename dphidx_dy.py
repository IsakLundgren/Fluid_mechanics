def dphidx_dy(x_2d,y_2d,phi_2d):
   import numpy as np

   dphidx_p=phi_2d
   dphidy_p=phi_2d

#%%%%%%%%%%% west face
   sx_w=np.diff(x_2d,axis=1)
   sy_w=np.diff(y_2d,axis=1)
# duplicate last column and put it at the end
   sx_w=np.insert(sx_w,-1,sx_w[:,-1],axis=1)
   sy_w=np.insert(sy_w,-1,sy_w[:,-1],axis=1)

# normalize
   d=np.sqrt(sx_w**2+sy_w**2)
   nx_w=-sy_w/d
   ny_w=sx_w/d
   d_w=d
   phi_w_2d=0.5*(np.add(phi_2d,np.roll(phi_2d,-1,axis=0)))
# delete last row
   phi_w_2d = np.delete(phi_w_2d, -1, 0)
# delete last columns
   phi_w_2d = np.delete(phi_w_2d, -1, 1)

# fix west boundary
   phi_w_2d[0,:]=phi_2d[0,0:-1]
# fix east boundary
   phi_w_2d[-1,:]=phi_2d[-1,0:-1]

#%%%%%%%%%%% south face
   sx_s=np.diff(x_2d,axis=0)
   sy_s=np.diff(y_2d,axis=0)
# duplicate last row and put it at the end
   sx_s=np.insert(sx_s,-1,sx_s[-1],axis=0)
   sy_s=np.insert(sy_s,-1,sy_s[-1],axis=0)

# normalize
   d=np.sqrt(sx_s**2+sy_s**2)
   nx_s=sy_s/d
   ny_s=-sx_s/d
   d_s=d
   phi_s_2d=0.5*(np.add(phi_2d,np.roll(phi_2d,-1,axis=1)))
# delete last row
   phi_s_2d = np.delete(phi_s_2d, -1, 0)
# delete last columns
   phi_s_2d = np.delete(phi_s_2d, -1, 1)
# fix south boundary
   phi_s_2d[:,0]=phi_2d[0:-1,0]
# fix north boundary
   phi_s_2d[:,-1]=phi_2d[0:-1,-1]

# area approaximated as the vector product of two triangles for cells
   ax=sx_w
   ay=sy_w
   bx=sx_s
   by=sy_s
   area_p1=0.5*np.absolute(ax[0:-1,0:-1]*by[0:-1,0:-1]-ay[0:-1,0:-1]*bx[0:-1,0:-1])
   
   ax=np.roll(ax,1,axis=0)
   ay=np.roll(ay,1,axis=0)
   bx=np.roll(bx,1,axis=1)
   by=np.roll(by,1,axis=1)
   area_p2=0.5*np.absolute(ax[0:-1,0:-1]*by[0:-1,0:-1]-ay[0:-1,0:-1]*bx[0:-1,0:-1])
   
   area_p=area_p1+area_p2


# west face's contribution to Gauss integral
   wx=phi_w_2d*nx_w*d_w
   wy=phi_w_2d*ny_w*d_w
# south face's contribution to Gauss integral
   sx=phi_s_2d*nx_s*d_s
   sy=phi_s_2d*ny_s*d_s

# gauss integration west + east face (neg. normal vector) for each cell
   ewx=-wx[1:,0:-1]+wx[0:-1,0:-1]
   ewy=-wy[1:,0:-1]+wy[0:-1,0:-1]

# gauss integration south + north face (neg. normal vector)
   nsx=-sx[0:-1,1:]+sx[0:-1,0:-1]
   nsy=-sy[0:-1,1:]+sy[0:-1,0:-1]

# gradient
   gradx=(ewx+nsx)/area_p
   grady=(ewy+nsy)/area_p

   dphidx_p=gradx
   dphidy_p=grady

####################### Neumann b.c. at boundaries
# duplicate row 0 and put it before row 0 (west boundary)
   dphidx_p=np.insert(dphidx_p,0,dphidx_p[0,:],axis=0)
   dphidy_p=np.insert(dphidy_p,0,dphidy_p[0,:],axis=0)
# duplicate last row and put it at the end (east boundary)
   dphidx_p=np.insert(dphidx_p,-1,dphidx_p[-1,:],axis=0)
   dphidy_p=np.insert(dphidy_p,-1,dphidy_p[-1,:],axis=0)
# duplicate first column and put it at the end
   dphidx_p=np.insert(dphidx_p,0,dphidx_p[:,0],axis=1)
   dphidy_p=np.insert(dphidy_p,0,dphidy_p[:,0],axis=1)
# duplicate last column and put it at the end
   dphidx_p=np.insert(dphidx_p,-1,dphidx_p[:,-1],axis=1)
   dphidy_p=np.insert(dphidy_p,-1,dphidy_p[:,-1],axis=1)

   return dphidx_p,dphidy_p




