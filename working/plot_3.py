import matplotlib.pyplot as plt
import numpy as np

k=3
d=3

M=80+1
N=40+1

r=0.15

X=np.zeros((N*d,M*k),dtype=float)
Y=np.zeros((N*d,M*k),dtype=float)
U=np.zeros((N*d,M*k),dtype=float)
V=np.zeros((N*d,M*k),dtype=float)
P=np.zeros((N*d,M*k),dtype=float)
T=np.zeros((N*d,M*k),dtype=float)

for i in range(0,k*d):
	
	if(i<=9):
		x=np.loadtxt('X0%d.dat' % i)
		y=np.loadtxt('Y0%d.dat' % i)
		u=np.loadtxt('U0%d.dat' % i)
		v=np.loadtxt('V0%d.dat' % i)
        	p=np.loadtxt('P0%d.dat' % i)
                t=np.loadtxt('T0%d.dat' % i)
	else:
		x=np.loadtxt('X%d.dat' % i)
		y=np.loadtxt('Y%d.dat' % i)
		u=np.loadtxt('U%d.dat' % i)
		v=np.loadtxt('V%d.dat' % i)
		p=np.loadtxt('P%d.dat' % i)
                t=np.loadtxt('T%d.dat' % i)
        
	x=np.reshape(x,[N,M])
	y=np.reshape(y,[N,M])
	u=np.reshape(u,[N,M])
	v=np.reshape(v,[N,M])
        p=np.reshape(p,[N,M])
        t=np.reshape(t,[N,M])

	X[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=x
	Y[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=y
	U[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=u
	V[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=v
        P[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=p
	T[(i/k)*N:(i/k)*N+N,(i%k)*M:(i%k)*M+M]=t

x_c = np.linspace(-r,r,200)
y_c = np.sqrt(r**2-x_c**2)

x_circle = np.concatenate([x_c,np.fliplr([x_c[:-1]])[0]])
y_circle = np.concatenate([y_c,-np.fliplr([y_c[:-1]])[0]])

r = 0.10
x_c = np.linspace(-r,r,200)
y_c = np.sqrt(r**2-x_c**2)

x_circle2 = np.concatenate([x_c,np.fliplr([x_c[:-1]])[0]])
y_circle2 = np.concatenate([y_c,-np.fliplr([y_c[:-1]])[0]])

#x_circle2 = x_circle.copy()
#y_circle2 = y_circle.copy()

x_circle = x_circle-0.40
x_circle2= x_circle2-0.10

plt.figure()
plt.title('Resultant Velocity')
plt.contourf(X,Y,np.sqrt(U**2+V**2),density=20)
#plt.quiver(X,Y,U,V)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.fill(x_circle,y_circle,'w')
#plt.fill(x_circle2,y_circle2,'w')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal') 

plt.figure()
plt.title('Pressure')
plt.contourf(X,Y,P,density=5)
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.fill(x_circle,y_circle,'w')
#plt.fill(x_circle2,y_circle2,'w')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal') 

plt.figure()
plt.title('Temperature')
plt.contourf(X,Y,T,density=5)
#plt.streamplot(X,Y,U,V,density=4,color='k')
plt.plot(X[:,0],Y[:,0],'k')
plt.plot(X[:,-1],Y[:,-1],'k')
plt.plot(X[0,:],Y[0,:],'k')
plt.plot(X[-1,:],Y[-1,:],'k')
plt.plot(x_circle,y_circle,'k')
#plt.plot(x_circle2,y_circle2,'k')
#plt.plot(X,Y,'g')
#plt.plot(X.T,Y.T,'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.axis('equal')
#plt.show()

#plt.figure()
#plt.title('Velocity Profile')
#plt.plot(Y[:,(M-1)/2],U[:,(M-1)/2],'k')
plt.show()
