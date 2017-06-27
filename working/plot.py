import matplotlib.pyplot as plt
import numpy as np

k=2
d=2

M=20+1
N=20+1

X=np.zeros((N*d,M*k),dtype=float)
Y=np.zeros((N*d,M*k),dtype=float)
U=np.zeros((N*d,M*k),dtype=float)
V=np.zeros((N*d,M*k),dtype=float)
P=np.zeros((N*d,M*k),dtype=float)
T=np.zeros((N*d,M*k),dtype=float)
W=np.zeros((N*d,M*k),dtype=float)
R=np.zeros((N*d,M*k),dtype=float)

for i in range(0,k*d):
	
	if(i<=9):
		x=np.loadtxt('X000%d.dat' % i)
		y=np.loadtxt('Y000%d.dat' % i)
		u=np.loadtxt('U000%d.dat' % i)
		v=np.loadtxt('V000%d.dat' % i)
        	p=np.loadtxt('P000%d.dat' % i)
                t=np.loadtxt('T000%d.dat' % i)
		w=np.loadtxt('W000%d.dat' % i)
		r=np.loadtxt('R000%d.dat' % i)
	elif(i<=99):
		x=np.loadtxt('X00%d.dat' % i)
		y=np.loadtxt('Y00%d.dat' % i)
		u=np.loadtxt('U00%d.dat' % i)
		v=np.loadtxt('V00%d.dat' % i)
		p=np.loadtxt('P00%d.dat' % i)
                t=np.loadtxt('T00%d.dat' % i)
		w=np.loadtxt('W00%d.dat' % i)
		r=np.loadtxt('R00%d.dat' % i) 
	elif(i<=999):
		x=np.loadtxt('X0%d.dat' % i)
		y=np.loadtxt('Y0%d.dat' % i)
		u=np.loadtxt('U0%d.dat' % i)
		v=np.loadtxt('V0%d.dat' % i)
		p=np.loadtxt('P0%d.dat' % i)
                t=np.loadtxt('T0%d.dat' % i)
		w=np.loadtxt('W0%d.dat' % i)
		r=np.loadtxt('R0%d.dat' % i)
	else:
		x=np.loadtxt('X%d.dat' % i)
		y=np.loadtxt('Y%d.dat' % i)
		u=np.loadtxt('U%d.dat' % i)
		v=np.loadtxt('V%d.dat' % i)
		p=np.loadtxt('P%d.dat' % i)
                t=np.loadtxt('T%d.dat' % i)
		w=np.loadtxt('W%d.dat' % i)
		r=np.loadtxt('R%d.dat' % i)

	x=np.reshape(x,[N,M])
	y=np.reshape(y,[N,M])
	u=np.reshape(u,[N,M])
	v=np.reshape(v,[N,M])
        p=np.reshape(p,[N,M])
        t=np.reshape(t,[N,M])
	w=np.reshape(w,[N,M])
	r=np.reshape(r,[N,M])

	X[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=x
	Y[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=y
	U[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=u
	V[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=v
        P[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=p
	T[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=t
	W[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=w
	R[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=r

plt.figure()
plt.title('Velocity Vector')
plt.quiver(X,Y,U,V)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal') 

plt.figure()
plt.title('V Velocity')
plt.contourf(X,Y,V,density=20)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.plot(X,Y,'g')
plt.plot(X.T,Y.T,'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.axis('equal')

plt.figure()
plt.title('Pressure')
plt.contourf(X,Y,P,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.axis('equal') 

plt.figure()
plt.title('Temperature')
plt.contourf(X,Y,T,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.axis('equal')

plt.figure()
plt.title('Viscosity')
plt.contourf(X,Y,W,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.axis('equal')

plt.figure()
plt.title('Density')
plt.contourf(X,Y,R,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.xlabel('X')
plt.ylabel('Y')
plt.colorbar()
plt.axis('equal')

plt.show()
