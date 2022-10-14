import scipy.io as sio
import numpy as np
import csv
import matplotlib.pyplot as plt
from dphidx_dy import dphidx_dy
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
plt.rcParams.update({'font.size': 22})


#--------------------------------------------------------------------------
ni=200 #Do not change it.
nj=200 #Do not change it.
#--------------------------------------------------------------------------
tt=np.arange(nj+1)
tt[1]=int(0)
for j  in range (2,nj+1):
    tt[j]=tt[j-1]+j

ss=np.zeros((ni+1,nj+1),dtype=int)
ss[:,1]= tt
for j  in range (1,nj+1):
   for i  in range (2,ni+1):
        ss[j,i]=ss[j,i-1]+i-1+(j-1)

pp=np.zeros((ni),dtype=int)
for i  in range (1,ni):
    pp[i] = int((i)**2)

subtx=np.ones((ni+1,nj+1),dtype=int)
count=0
for n in range (nj,0,-1):
    count = nj-n+1
    subtx[1:count+1,n]=0
    subtx[count+1:,n]=pp[1:ni-count+1]

gridIndex=ss-subtx
#--------------------------------------------------------------------------
nmax=ni*nj
u=np.zeros(nmax)
v=np.zeros(nmax)
p=np.zeros(nmax)
te=np.zeros(nmax)
diss=np.zeros(nmax)
vist=np.zeros(nmax)
x=np.zeros(nmax)
y=np.zeros(nmax)
with open('./CSVOutputs/output_standard-keps-low-re.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = 0  
    for row in csv_reader:
        if n == 0:
            print(f'Column names are {", ".join(row)}')
            n += 1
        else:
#           if n < 10:  # print the 10 first lines
                #print(f'\t{row[0]}: U={row[1]}, V={row[2]},  P={row[3]},  k={row[4]},  eps={row[5]}, vist={row[6]}, x={row[7]}, y={row[8]}')

            u[n-1]=row[1]
            v[n-1]=row[2]
            p[n-1]=row[3]
            te[n-1]=row[4]
            diss[n-1]=row[5]
            vist[n-1]=row[6]
            x[n-1]=row[7]
            y[n-1]=row[8]
            n += 1
        #if n == 10:
           #break

print(f'Processed {n} lines.')


x1_2d=np.zeros((ni,nj))
x2_2d=np.zeros((ni,nj))
v1_2d=np.zeros((ni,nj))
p_2d=np.zeros((ni,nj))
v2_2d=np.zeros((ni,nj))
te_2d=np.zeros((ni,nj))
vist_2d=np.zeros((ni,nj))
diss_2d=np.zeros((ni,nj))
vist_2d=np.zeros((ni,nj))

for j  in range (1,nj+1):
   for i  in range (1,ni+1):
       n=gridIndex[i,j]
       x1_2d[i-1,j-1]=x[n]
       x2_2d[i-1,j-1]=y[n]
       p_2d[i-1,j-1]=p[n]
       v1_2d[i-1,j-1]=u[n]
       v2_2d[i-1,j-1]=v[n]
       te_2d[i-1,j-1]=te[n]
       vist_2d[i-1,j-1]=vist[n]
       diss_2d[i-1,j-1]=diss[n]

v1_2d_org=v1_2d
x2_2d_org=x2_2d

x1_2d=np.flipud(np.transpose(x1_2d))
x2_2d=np.flipud(np.transpose(x2_2d))
v1_2d=np.flipud(np.transpose(v1_2d))
v2_2d=np.flipud(np.transpose(v2_2d))
p_2d=np.flipud(np.transpose(p_2d))
te_2d=np.flipud(np.transpose(te_2d))
vist_2d=np.flipud(np.transpose(vist_2d))
diss_2d=np.flipud(np.transpose(diss_2d))


# The STAR-CCM data do no include the boundaries. 
# below we add wall at the bottom (low x_2, south) and top (high x_2, north)
deltaYBottom=(x2_2d[:,1]-x2_2d[:,0])/2
deltaYTop=(x2_2d[:,-1]-x2_2d[:,-1-1])/2


# duplicate first column  (south boundary)
x1_2d=np.insert(x1_2d,0,x1_2d[:,0],axis=1)
x2_2d=np.insert(x2_2d,0,x2_2d[:,0],axis=1)
v1_2d=np.insert(v1_2d,0,v1_2d[:,0],axis=1)
v2_2d=np.insert(v2_2d,0,v2_2d[:,0],axis=1)
te_2d=np.insert(te_2d,0,te_2d[:,0],axis=1)
diss_2d=np.insert(diss_2d,0,diss_2d[:,0],axis=1)
vist_2d=np.insert(vist_2d,0,vist_2d[:,0],axis=1)
p_2d=np.insert(p_2d,0,p_2d[:,0],axis=1)
diss_2d=np.insert(diss_2d,0,diss_2d[:,0],axis=1)
zero_col=np.zeros(ni)

# set south boundary to zero
v1_2d[:,0]=0.
v2_2d[:,0]=0.
te_2d[:,0]=0.
diss_2d[:,0]=0.
vist_2d[:,0]=0.

x2_2d[:,0]=x2_2d[:,1]-deltaYBottom;

# duplicate last column and put it at the end (north boundary)
x1_2d=np.insert(x1_2d,-1,x1_2d[:,-1],axis=1)
x2_2d=np.insert(x2_2d,-1,x2_2d[:,-1],axis=1)
v1_2d=np.insert(v1_2d,-1,v1_2d[:,-1],axis=1)
v2_2d=np.insert(v2_2d,-1,v2_2d[:,-1],axis=1)
te_2d=np.insert(te_2d,-1,te_2d[:,-1],axis=1)
diss_2d=np.insert(diss_2d,-1,diss_2d[:,-1],axis=1)
vist_2d=np.insert(vist_2d,-1,vist_2d[:,-1],axis=1)
p_2d=np.insert(p_2d,-1,p_2d[:,-1],axis=1)
diss_2d=np.insert(diss_2d,-1,diss_2d[:,-1],axis=1)

# append a column with zeros and put it at the end (north boundary)

nj=nj+2
#--------------------------------------------------------------------------
#*************** DO NOT CHANGE ANY PART OF THE ABOVE LINES. ***************
#--------------------------------------------------------------------------
#
hmax=0.050 # Maximum hill height.
H=3.035*hmax # Cahnnel height.
L=9*hmax # Space between two hills summit. 

#**** LOADING MEASUREMENT DATA AT DIFFERENT X_1 (STREAMWISE) LOCATIONS. ****
#--------------------------------------------------------------------------
xh1=np.genfromtxt("./XYInputs/xh1.xy", comments="%")
y_1=xh1[:,0] # x_2 coordinates, wall-normal direction.
v1_Exp_1=xh1[:,1] # mean velocity in the streamwise direction (x_1) along wall-normal direction (x_2). 
v2_Exp_1=xh1[:,2] # mean velocity in the streamwise direction (x_1) along wall-normal direction (x_2). 
uu_Exp_1=xh1[:,3] # Normal Reynolds stress (Re_xx) along wall-normal direction (x_2).  
vv_Exp_1=xh1[:,4] # Normal Reynolds stress (Re_yy) along wall-normal direction (x_2).
uv_Exp_1=xh1[:,5] # Shear Reynolds stress (Re_xy) along wall-normal direction (x_2).
# The locations for the measurement data are: x/h=0.05, 0.5, 1, 2, 3, 4, 5, 6, 7 and 8.
# You should find appropriate "i" corresponds to measurement x locations.
#For example, "xh005.xy", "xh05.xy" and "xh1.xy" are the measurment data at x/h=0.05,x/h=0.5 and x/h=1, repectively.

#################################### plot v_1 vs. x_2 at x_1=hmax
fig1,ax1 = plt.subplots()
xx=hmax
i1 = (np.abs(xx-x1_2d[:,1])).argmin()  # find index which closest fits xx
plt.plot(v1_2d[i1,:],x2_2d[i1,:],'b-')
plt.plot(v1_Exp_1,y_1,'bo')
plt.xlabel("$V_1$")
plt.ylabel("$x_2$")
plt.title("Velocity")
plt.axis([-0.1,0.6,0.0225,H+0.01])
# Create inset of width 30% and height 40% of the parent axes' bounding box
# at the lower left corner (loc=3)
# upper left corner (loc=2)
# use borderpad=1, i.e.
# 22 points padding (as 22pt is the default fontsize) to the parent axes
axins1 = inset_axes(ax1, width="40%", height="30%", loc=2, borderpad=2)
plt.plot(v1_2d[i1,:],x2_2d[i1,:],'b-')
plt.plot(v1_Exp_1,y_1,'bo')
plt.axis([-0.1,0.01,0.0225,0.04])
# reduce fotnsize 
axins1.tick_params(axis = 'both', which = 'major', labelsize = 10)
# Turn ticklabels of insets off
#axins1.tick_params(labelleft=False, labelbottom=False)
# put numbers on the right y axis
axins1.yaxis.set_label_position("right")
axins1.yaxis.tick_right()

plt.savefig('Vel_python.eps')

#################################### contour p
fig2 = plt.figure("Figure 2")
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,p_2d, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("contour pressure plot")
plt.savefig('p_contour.eps')

fig3 = plt.figure("Figure 3")
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,np.sqrt(v1_2d**2+v2_2d**2), 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("contour velocity magnitude plot")
plt.savefig('vel_contour.eps')


# compute velociy gradients
dv1dx1_2d= np.zeros((ni,nj))
dv1dx2_2d= np.zeros((ni,nj))
dv2dx1_2d= np.zeros((ni,nj))
dv2dx2_2d= np.zeros((ni,nj))

# note that the routine 'dphidx_dy' wants x_1 and x_2 at faces (of size (ni-1)x(nj-1))
# Here we cheat a bit and give the coordinate at cell center (but using correct size (ni-1)x(nj-1))
dv1dx1_2d,dv1dx2_2d = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1],v1_2d)
dv2dx1_2d,dv2dx2_2d = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1],v2_2d)


#################################### vector plot
fig4 = plt.figure("Figure 4")
plt.clf() #clear the figure
k=6# plot every forth vector
ss=3.2 #vector length
plt.quiver(x1_2d[::k,::k],x2_2d[::k,::k],v1_2d[::k,::k],v2_2d[::k,::k],width=0.01)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("vector plot")
plt.savefig('vect_python.eps')
#

#Bernoulli velocity and pressure
def V_b(x_i):
    a = np.abs(x2_2d[x_i,0] - x2_2d[x_i,-1])
    b = np.trapz(v1_2d[x_i,:],x2_2d[x_i,:])
    return 1 / a * b

rho = 998.29 #kg/m^3
P_bern = np.zeros(ni)
for i in range(0,ni):
    P_bern[i] = V_b(i)**2 /2 * rho
