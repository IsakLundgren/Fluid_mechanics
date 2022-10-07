import numpy as np
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size': 22})

plt.interactive(True)

# Channel flow
data=np.genfromtxt("channel_flow_data.dat", comments="%")

ni=199  # number of grid nodes in x_1 direction
nj=28  # number of grid nodes in x_2 direction
x1=data[:,0] #don't use this array
x2=data[:,1] #don't use this array
v1=data[:,2] #don't use this array
v2=data[:,3] #don't use this array
p=data[:,4]  #don't use this array

# transform the arrays from 1D fields x(n) to 2D fields x(i,j)
# the first index 'i', correponds to the x-direction
# the second index 'j', correponds to the y-direction

x1_2d=np.reshape(x1,(nj,ni)) #this is x_1 (streamwise coordinate)
x2_2d=np.reshape(x2,(nj,ni)) #this is x_2 (wall-normal coordinate)
v1_2d=np.reshape(v1,(nj,ni)) #this is v_1 (streamwise velocity component)
v2_2d=np.reshape(v2,(nj,ni)) #this is v_1 (wall-normal velocity component)
p_2d=np.reshape(p,(ni,nj))   #this is p   (pressure)

x1_2d=np.transpose(x1_2d)
x2_2d=np.transpose(x2_2d)
v1_2d=np.transpose(v1_2d)
v2_2d=np.transpose(v2_2d)
p_2d=np.transpose(p_2d)
dv1_dx2=np.zeros((ni,nj))
for i in range(0,ni-1):
   for j in range(0,nj-1):
      dx2=x2_2d[i,j+1]-x2_2d[i,j-1]
      dv1_dx2[i,j]=(v1_2d[i,j+1]-v1_2d[i,j-1])/dx2

# fix the derivative at the walls
for i in range(0,ni-1):
# lower wall
   dx2=x2_2d[i,1]-x2_2d[i,0]
   dv1_dx2[i,0]=(v1_2d[i,1]-v1_2d[i,0])/dx2
# upper wall
   dx2=x2_2d[i,nj-1]-x2_2d[i,nj-2]
   dv1_dx2[i,nj-1]=(v1_2d[i,nj-1]-v1_2d[i,nj-2])/dx2

# you can also use the built-in command
x1_1d=x1_2d[:,1] # make 1d array
x2_1d=x2_2d[1,:] # make 1d array
dv1_dx1_built_in, dv1_dx2_bulit_in=np.gradient(v1_2d,x1_1d,x2_1d)

############################## D1
print("The fully developed area:")
L = 0.6385 #m
LConv = L / ni

h = 0.01 #m
hConv = h / nj

#Distance with dv1dx1 = 0.01
dv1_dx1_center = np.abs(dv1_dx1_built_in[:,int(nj / 2)])
cond = 0.01

tempMin = dv1_dx1_center - cond
stringentCondIndex=np.argmin(tempMin)

d1result1 = np.argmin(tempMin) * LConv #Fully developed area distance from inlet
print("can be found at ", d1result1," meters from the inlet using the velocity gradient criterion.") 

#Distance with v = vmax * 0.99
cond = 0.99
tempV = np.abs(v1_2d[:,int(nj / 2)])
tempMin = tempV - np.max(tempV) * cond

d1result2 = np.argmin(np.abs(tempMin)) * LConv #Fully developed area distance from inlet
print("can be found at ", d1result2," meters from the inlet using the maximum velocity criterion.") 

#Analytical soln
nu = 1.6e-5
d1result3 = 0.016 * np.mean(tempV) * (2 * h)**2 / nu
print("can be found at ", d1result3," meters from the inlet using an analytical approximation.") 

#v2 at x2 = h/4 for fully developed region
d1result4 = v2_2d[int(d1result1 / LConv),int(h / 4 / hConv)]
print("has a vertical component of ", d1result4, " m/s at $x_2 = h/4$")

############################### D2

mu = 1.925 * 10**-5 #m**2/s

#Wall shear on the lower wall
wallShearL = mu * dv1_dx2[:, 0]

#Wall shear on the upper wall
n = np.array([0,-1])
#The transversial component of tij is given by 2 mu Sij which is just mu * dv1_dx2[:, nj - 1]

wallShearU = mu * dv1_dx2[:, nj - 1] * n[1]

fig1,ax1 = plt.subplots()
plt.plot(x1_2d[:,0], wallShearL, '-b',label='wallShearL')
plt.plot(x1_2d[:,nj -1], wallShearU,  '--r',label='wallShearU')
plt.ylabel('$t_1$ [N/m^2]') 
plt.xlabel('$x_1 [m]$') 
leg = ax1.legend();

anVal = h /2 * (p_2d[14,0] - p_2d[14,-1])/L

print("D2, numerical value: ", wallShearL[-2], ", analytical value: ", anVal)

############################### D3

fig1,ax1 = plt.subplots()
plt.plot(x1_2d[:,int(nj / 8)], v1_2d[:,int(nj / 8)],label='Near wall') #Near wall
plt.plot(x1_2d[:,int(nj / 2)], v1_2d[:,int(nj / 2)],label='Center') #Center
plt.ylabel('$V_1[m/2]$') 
plt.xlabel('$x_1[m]$')
leg = ax1.legend();

fig1, ax1 = plt.subplots()
plt.plot(x1_2d[:,int(nj / 8)], dv1_dx1_built_in[:,int(nj / 8)],label='Near wall') #Near wall
plt.plot(x1_2d[:,int(nj / 2)], dv1_dx1_built_in[:,int(nj / 2)],label='Center') #Center
plt.ylabel('$dV_1 / dX_1 [1/s]$') 
plt.xlabel('$x_1[m]$')
leg = ax1.legend();

#numerical integral
xi = np.zeros(ni)
for i in range(ni):
   sum = 0
   for j in range(nj - 1): ##using left hand side truncation
      sum += v1_2d[i,j] * (x2_2d[i,j+1] - x2_2d[i,j])
   xi[i] = sum

fig1, ax1 = plt.subplots()
plt.plot(x1_2d[:,0], xi)
plt.ylabel('ξ') 
plt.xlabel('$x_1 [m]$')

############################### D4

fig1, ax1 = plt.subplots()
plt.plot(x1_2d[:,0], v2_2d[:,int(nj / 10)])
plt.ylabel('$dV_2 [m/2]$') 
plt.xlabel('$x_1 [m]$')

############################### D5
#w1,w2 zero in entire domain
#w3 = dv2_dx1 - dv1_dx2
#cross terms need to be recalculated

x1_1d=x1_2d[:,1] # make 1d array
x2_1d=x2_2d[1,:] # make 1d array
dv1_dx1, dv1_dx2=np.gradient(v1_2d,x1_1d,x2_1d)
dv2_dx1, dv2_dx2=np.gradient(v2_2d,x1_1d,x2_1d)

w3 = dv2_dx1 - dv1_dx2

#Most stringent fully developed conditions: 
fig1, ax1 = plt.subplots()
plt.plot(x2_1d, w3[stringentCondIndex, :], '-b', label='Fully Developed') #Fully dev.
plt.plot(x2_1d, w3[0, :], '-r',label='Inlet') #Inlet.
for i in range(0,8):
   plt.plot(x2_1d, w3[int( i * np.argmin(np.abs(tempMin)) /16), :], ':g')

plt.ylabel('$Vorticity [1/s]$') 
plt.xlabel('$x_1 [m]$')
leg = ax1.legend()

############################### D6

S12 = 1/2 * (dv2_dx1 + dv1_dx2)
O12 = 1/2 * (dv2_dx1 - dv1_dx2)

fig1, ax1 = plt.subplots()

plt.contourf(x1_2d, x2_2d, S12, 100)
plt.ylabel('$dv_2 dx_1 + dv_1 dx_2$') 
plt.xlabel('$x_1$')

fig1, ax1 = plt.subplots()
plt.contourf(x1_2d, x2_2d, O12, 100)
plt.ylabel('$dv_2 dx_1 - dv_1 dx_2$') 
plt.xlabel('$x_1$')
############################### D7

S11 = dv1_dx1
S21 = S12
S22 = dv2_dx2
SKK = S11 + S22 #0 by incompressibility

tau_11 = 2 * mu * S11
tau_12 = 2 * mu * S12
tau_21 = 2 * mu * S21
tau_22 = 2 * mu * S22

dissipation=np.zeros((ni,nj))
for i in range(0,ni):
   for j in range(0,nj):
      dissipation[i,j] = tau_11[i,j] * dv1_dx1[i,j] + tau_12[i,j] * dv1_dx2[i,j] + tau_21[i,j] * dv2_dx1[i,j] + tau_22[i,j] * dv2_dx2[i,j]

#Numerical integration
dispInt = np.trapz(np.trapz(dissipation, x2_2d[0,:]), x1_2d[:,0])
vInt = np.trapz(v1_2d[0,:], x2_1d)

rho = 1.204 #kg/m3
c_p = 1006 #kJ/(kgK)
d7result1 = dispInt / (rho * c_p * vInt)

print('The change in bulk temperature is ', d7result1, ' degrees Kelvin')
############################### D8

EigVals=np.zeros((ni,nj,2))
EigVecs=np.zeros((ni,nj,2,2))
for i in range(0,ni):
   for j in range(0,nj):
      tau = np.array([[tau_11[i,j], tau_12[i,j]], [tau_21[i,j], tau_22[i,j]]])
      EigVals[i,j,:], EigVecs[i,j,:,:] = np.linalg.eig(tau)

fig1, ax1 = plt.subplots()
plt.plot(x2_1d, EigVals[stringentCondIndex, :, 0], label='1st Eigen Value')
plt.plot(x2_1d, EigVals[stringentCondIndex, :, 1], label='2nd Eigen Value')
plt.ylabel('$x_2[m]$') 
plt.xlabel('$x_1[m]$')
leg = ax1.legend();
#plt.plot(x2_1d, tau_11[stringentCondIndex, :]) removed because 0
plt.plot(x2_1d, tau_12[stringentCondIndex, :], label='\u03C4_12 & \u03C4_21')
plt.ylabel('$x_2[m]$') 
plt.xlabel('$x_1[m]$')
leg = ax1.legend();
#plt.plot(x2_1d, tau_22[stringentCondIndex, :]) removed because 0

################################## D9

fig1, ax1 = plt.subplots()
plt.quiver(x1_2d, x2_2d, EigVecs[:,:,0,0], EigVecs[:,:,1,0], scale = 5, scale_units ='inches')
plt.quiver(x1_2d, x2_2d, EigVecs[:,:,0,1], EigVecs[:,:,1,1], scale = 5, scale_units ='inches')
plt.ylabel('$x_2[m]$') 
plt.xlabel('$x_1[m]$')

#Construct t vector
t_1=np.zeros((ni,nj,2))
t_2=np.zeros((ni,nj,2))
for i in range(0,ni):
   for j in range(0,nj):
      t_1[i,j,0] = EigVals[i,j,0]*EigVecs[i,j,0,0]
      t_1[i,j,1] = EigVals[i,j,0]*EigVecs[i,j,1,0]
      t_2[i,j,0] = EigVals[i,j,1]*EigVecs[i,j,0,1]
      t_2[i,j,1] = EigVals[i,j,1]*EigVecs[i,j,1,1]
fig1, ax1 = plt.subplots()
plt.quiver(x1_2d, x2_2d, t_1[:,:,0], t_1[:,:,1], scale = None, scale_units ='inches')
plt.quiver(x1_2d, x2_2d, t_2[:,:,0], t_2[:,:,1], scale = None, scale_units ='inches')
plt.ylabel('$x_2[m]$') 
plt.xlabel('$x_1[m]$')
#plt.show(block=True)