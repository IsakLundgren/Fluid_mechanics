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
            #print(f'Column names are {", ".join(row)}')
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
fig1,ax1 = plt.subplots(figsize=(13,6.5))
xx=hmax
i1 = (np.abs(xx-x1_2d[:,1])).argmin()  # find index which closest fits xx
plt.rcParams['font.size'] = '30'
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

plt.savefig('Vel_python.eps', bbox_inches = 'tight')

#################################### contour p
fig2 = plt.figure("Figure 2",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,p_2d, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("contour pressure plot")
plt.colorbar()
plt.savefig('p_contour.eps', bbox_inches = 'tight')

fig3 = plt.figure("Figure 3",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,np.sqrt(v1_2d**2+v2_2d**2), 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("contour velocity magnitude plot")
plt.colorbar()
plt.savefig('vel_contour.eps', bbox_inches = 'tight')

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
fig4 = plt.figure("Figure 4",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
k=6# plot every forth vector
ss=3.2 #vector length
plt.quiver(x1_2d[::k,::k],x2_2d[::k,::k],v1_2d[::k,::k],v2_2d[::k,::k],width=0.01)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("vector plot")
plt.savefig('vect_python.eps', bbox_inches = 'tight')
#

#########################################Bernoulli velocity and pressure

def V_b(x_i):
    a = np.abs(x2_2d[x_i,0] - x2_2d[x_i,-1])
    b = np.trapz(v1_2d[x_i,:],x2_2d[x_i,:])
    return 1 / a * b

fig5 = plt.figure("Figure 5",figsize=(13,6.5))
plt.rcParams['font.size'] = '20'
plt.clf() #clear the figure
plt.plot(x1_2d[:,0], V_b(range(0,ni)))
plt.legend()
plt.xlabel("$x [m]$")
plt.ylabel("$v [m/s]$")
plt.title("Bulk velocity")
plt.savefig('BulkVel.eps', bbox_inches = 'tight')

rho = 998.29 #kg/m^3
P_bern = np.zeros(ni)
for i in range(0,ni):
    P_bern[i] = (V_b(0)**2-V_b(i)**2) /2 * rho

#Starccm bulk pressure
def P_b(x_i):
    a = np.abs(x2_2d[x_i,0] - x2_2d[x_i,-1])
    b = np.trapz(p_2d[x_i,:],x2_2d[x_i,:])
    return 1/a * b

P_star = np.zeros(ni)
for i in range(0,ni):
    P_star[i] = P_b(i) - P_b(0)

fig5 = plt.figure("Figure 5",figsize=(13,6.5))
plt.rcParams['font.size'] = '20'
plt.clf() #clear the figure
plt.plot(x1_2d[:,0], P_bern, label='Bernoulli equation')
plt.plot(x1_2d[:,0], P_star, label='Bulk pressure')
plt.legend()
plt.xlabel("$x [m]$")
plt.ylabel("$P [Pa]$")
plt.title("Pressures from different approaches")
plt.savefig('Pressures.eps', bbox_inches = 'tight')

print('mean flow velocity: ', str(V_b(int(ni/2))))
print('The pressure drop according to the Bernoulli principle is', str(np.abs(P_bern[-1] - P_bern[0])), 'Pascals.')
print('The pressure drop according to a StarCCM+ simulation is', str(np.abs(P_star[-1] - P_star[0])), 'Pascals.')

###########################################Skin friction and pressure
ustar_top = np.zeros(ni)
ustar_bot = np.zeros(ni)
yplus_top = np.zeros(ni)
yplus_bot = np.zeros(ni)
shearstress_top = np.zeros(ni)
shearstress_bot = np.zeros(ni)
x_top = np.zeros(ni)
x_bot = np.zeros(ni)
y_top = np.zeros(ni)
y_bot = np.zeros(ni)
with open('./CSVOutputs/SkinFrictionTop.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = 0  
    for row in csv_reader:
        if n == 0:
            n += 1
        else:
            ustar_top[n-1]=row[0]
            yplus_top[n-1]=row[1]
            shearstress_top[n-1]=row[2]
            x_top[n-1]=row[3]
            y_top[n-1]=row[4]
            n += 1
with open('./CSVOutputs/SkinFrictionBot.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = 0  
    for row in csv_reader:
        if n == 0:
            n += 1
        else:
            ustar_bot[n-1]=row[0]
            yplus_bot[n-1]=row[1]
            shearstress_bot[n-1]=row[2]
            x_bot[n-1]=row[3]
            y_bot[n-1]=row[4]
            n += 1

C_f_top = np.zeros(ni)
C_f_bot = np.zeros(ni)
for i in range(0, ni):
    C_f_top[i] = shearstress_top[i] / (0.5 * rho * V_b(i)**2)
    C_f_bot[i] = shearstress_bot[i] / (0.5 * rho * V_b(i)**2)

fig5 = plt.figure("Figure 5",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.plot(x_top, C_f_top, label='Top wall')
plt.plot(x_bot, C_f_bot, label='Bottom wall')
plt.legend()
plt.xlabel("$x$")
plt.ylabel("$C_f$")
plt.title("Skin friction plot")
plt.savefig('skinFriction.eps', bbox_inches = 'tight')

#################################################Vorticity

w3 = dv2dx1_2d - dv1dx2_2d

fig6 = plt.figure("Figure 6",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,w3, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Contour voritcity plot")
plt.colorbar()
plt.savefig('vort_contour.eps', bbox_inches = 'tight')

##################################################Turbulent viscosity
#############Might have to be divided by rho######
mu = 0.001003 #Viscocity Pa * s
ratio = vist_2d / mu

fig7 = plt.figure("Figure 7",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,ratio, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Contour turbulent viscosity plot")
plt.colorbar()
plt.savefig('TurbViscRatio_contour.eps', bbox_inches = 'tight')

#Levels
fig8 = plt.figure("Figure 8",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf()
dni = ni / 6
for i in range(1,5): #Divide in 6 parts
    lab = 'x = ' + str(x1_2d[int(dni * i), 0])
    plt.plot(ratio[int(dni * i), :], (x2_2d[int(dni * i), :] - x2_2d[int(dni * i), 0]), label = lab)
plt.legend(fontsize=23)
plt.ylabel("$y$")
plt.xlabel("$\mu_t / \mu$")
plt.savefig('TurbVisc_levels.eps', bbox_inches = 'tight')


#y+
x_iYp = int(ni/2)
yplus_calc = ustar_bot[x_iYp] * (x2_2d[x_iYp,:] - x2_2d[x_iYp,0]) * rho / mu

fig9 = plt.figure("Figure 9",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.plot(ratio[x_iYp, :], yplus_calc)
plt.xlabel("$\mu_t / \mu$")
plt.ylabel("$y^+$")
plt.title("Turbulent viscosity on bottom wall")
plt.yscale('log')
plt.grid()
plt.axis([0, 450, 1, 1000])
plt.savefig('TurbVisc_yplus.eps', bbox_inches = 'tight')

#########################################################Turbulent Diffussion
nu_t = vist_2d / rho
symm1 = 2 * dv1dx1_2d
symm2 = 2 * dv2dx2_2d
assym = dv1dx2_2d + dv2dx1_2d
diffus1_1,temp1 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu_t * symm1)
temp2,diffus2_2 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu_t * symm2)
diffus2_1,diffus1_2 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu_t * assym)
diffus_t1 = diffus1_1 + diffus1_2
diffus_t2 = diffus2_1 + diffus2_2

nu = mu / rho
diffus1_1,temp1 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu * symm1)
temp2,diffus2_2 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu * symm2)
diffus2_1,diffus1_2 = dphidx_dy(x1_2d[0:-1,0:-1],x2_2d[0:-1,0:-1], nu * assym)
diffus_v1 = diffus1_1 + diffus1_2
diffus_v2 = diffus2_1 + diffus2_2

noOfPlots = 4
plt.rcParams['font.size'] = '10'
fig10, axs = plt.subplots(1,noOfPlots, figsize=(13,6.5))
fig10.suptitle('Viscous and Turbulent diffusions at varied x_1')
dni = ni / (noOfPlots + 1)
for i in range(0,noOfPlots):
    axs[i].plot(diffus_t1[int(dni * i + 1), :],x2_2d[int(dni * i + 1), :], label = 'Turbulent')
    axs[i].plot(diffus_v1[int(dni * i + 1), :],x2_2d[int(dni * i + 1), :], label = 'Viscous', linestyle='dotted')
plt.savefig('ViscDiffus_levels.eps', bbox_inches = 'tight')
plt.legend()

#y+
fig11 = plt.figure("Figure 11",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.plot(diffus_v1[x_iYp, :], yplus_calc, linestyle = 'dotted')
plt.plot(diffus_t1[x_iYp, :], yplus_calc)
plt.xlabel("Diffusion Magnitude m^2/s")
plt.ylabel("$y^+$")
plt.title("Turbulent viscosity on bottom wall")
plt.grid()
plt.savefig('ViscDiffus_yplus.eps', bbox_inches = 'tight')

#########################################################Production term
fig12 = plt.figure("Figure 12",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,te_2d, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Contour turbulent kinetic energy plot")
plt.colorbar()
plt.savefig('TurbKinEn_contour.eps', bbox_inches = 'tight')

P_k = vist_2d * (2*dv1dx1_2d**2 + (dv1dx2_2d + dv2dx1_2d)**2 + 2*dv2dx2_2d**2)
fig13 = plt.figure("Figure 13",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,P_k, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Contour $k_e$ production term plot")
plt.colorbar()
plt.savefig('KinEProd_contour.eps', bbox_inches = 'tight')

fig14 = plt.figure("Figure 14",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
plt.clf() #clear the figure
plt.contourf(x1_2d,x2_2d,vist_2d, 50)
plt.xlabel("$x$")
plt.ylabel("$y$")
plt.title("Contour turbulent viscosity plot")
plt.colorbar()
plt.savefig('TurbVisc_contour.eps', bbox_inches = 'tight')

###########################################################Wall boundary conditions epsilon

wDTop = np.ones(ni) * 1000
wDBot = np.ones(ni) * 1000
with open('./CSVOutputs/WallDistanceTop.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = 0  
    for row in csv_reader:
        if n == 0:
            n += 1
        else:
            wDTop[n-1]=row[0]
            n += 1
with open('./CSVOutputs/WallDistanceBot.csv') as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',')
    n = 0  
    for row in csv_reader:
        if n == 0:
            n += 1
        else:
            wDBot[n-1]=row[0]
            n += 1            

dispRateCalcTop = 2 * mu / rho * te_2d[:,-1] / (wDTop)**2
dispRateCalcBot = 2 * mu / rho * te_2d[:,1] / (wDBot)**2

print('The average value for the dissipation on both boundaries was ', str(np.average((diss_2d[:,1] + diss_2d[:,-1]) / 2)))
print('Biggest differences in dissipation rates calculated with python/starccm+')
print('Bottom: ', str(np.max(dispRateCalcBot - diss_2d[:,1])))
print('Top: ',str(np.max(dispRateCalcTop - diss_2d[:,-1])))

#############################################################Compare with experimental data for v_1

nOS = 7 #Number of stations
pts = np.zeros(nOS)
fileDir = ["./XYInputs/xh" for x in range(nOS)]

for i in range(nOS):
    pts[i] = i+1
pts[-1] = pts[-1] + 1

fig15 = plt.figure("Figure 15",figsize=(13,6.5))
plt.rcParams['font.size'] = '30'
for i in range(nOS):
    j = (np.abs(pts[i] * hmax-x1_2d[:,[i]])).argmin()
    fileDir[i] = fileDir[i] + str('{0:.0f}'.format(pts[i])) + ".xy"
    #print(fileDir[i])
    xh=np.genfromtxt(fileDir[i], comments="%")
    yExp = xh[:,0]
    v1Exp = xh[:,1]
    plt.plot(v1_2d[j,:],x2_2d[i1,:],'b-')
    plt.plot(v1Exp,yExp,'bo')

plt.xlabel("$V_1$")
plt.ylabel("$x_2$")
plt.title("Velocity")
plt.savefig('VelExpComp.eps', bbox_inches = 'tight')

plt.show(block=True)