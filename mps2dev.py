import numpy as np
import matplotlib.pyplot as plt

output=open("./output.dat", 'w+')

print("Funkcja: e^(-x^2-y)")
print("Funkcja: e^(-x^2-y)", file=output)
def f(x,y):
    return np.exp(-x**2-y)
def fprime(x,y):
    return -2*x*f(x,y), -f(x,y)
def fbis(x,y):
    return -2*f(x,y)-2*x*fprime(x,y)[0], -fprime(x,y)[0], f(x,y)

def AnglesToCartesian(angles, radiuses=1):
    return radiuses*np.cos(angles), radiuses*np.sin(angles)
def Weight(x,y):
	return 1/(x**2+y**2)
Circle=2*np.pi
NumberOfRegions=8
AngleRange=Circle/NumberOfRegions

RegionAngles=np.linspace(0,Circle,NumberOfRegions+1)

#plot regions
x,y = AnglesToCartesian(RegionAngles)
xlines = np.vstack((np.zeros_like(x), x))
ylines = np.vstack((np.zeros_like(y), y))
plt.plot(xlines, ylines, "k--", label="Regions")

#plot points
maxradius=1
minradius=0
RandomAngles = np.arange(0,Circle,AngleRange)+AngleRange*np.random.random(NumberOfRegions)
RandomRadii=(maxradius-minradius)*np.random.random(NumberOfRegions)+minradius
xpoint, ypoint = AnglesToCartesian(RandomAngles, RandomRadii)
xpointlines = np.vstack((np.zeros_like(xpoint), xpoint))
ypointlines = np.vstack((np.zeros_like(ypoint), ypoint))
plt.plot(xpointlines, ypointlines, "bo-", label="Wylosowane punkty")

#get function values
R = np.sum(RandomRadii)/NumberOfRegions #scaling radius
alpha_R=R/RandomRadii #scaling

xpointScaled = xpoint*alpha_R
ypointScaled = ypoint*alpha_R

plt.plot(xpointScaled, ypointScaled, "ro", label="Przeskalowane punkty")

FunctionValuesScaled=f(xpointScaled, ypointScaled)
FunctionAtZero=f(0,0)


alpha1=np.sum(xpointScaled**2)
alpha2=np.sum(xpointScaled*ypointScaled)
alpha3=np.sum(ypointScaled**2)
delta1=np.sum((FunctionValuesScaled-FunctionAtZero)*xpointScaled)
delta2=np.sum((FunctionValuesScaled-FunctionAtZero)*ypointScaled)

a = -(delta1*alpha3-delta2*alpha2)/(alpha2**2-alpha1*alpha3) #d/dx
b = -(delta2*alpha1-delta1*alpha2)/(alpha2**2-alpha1*alpha3) #d/dy

print("Wyniki analityczne")
print("%.4f, %.4f" %fprime(0,0))

print("Wartosci pochodnych ze skalowaniem")
print("%.4f, %.4f" %(a,b))
print("Roznice")
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])))

print("Wyniki analityczne", file=output)
print("%.4f, %.4f" %fprime(0,0), file=output)

print("Wartosci pochodnych ze skalowaniem", file=output)
print("%.4f, %.4f" %(a,b), file=output)
print("Roznice", file=output)
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])), file=output)

FunctionValues=f(xpoint, ypoint)

alpha1=np.sum(xpoint**2)
alpha2=np.sum(xpoint*ypoint)
alpha3=np.sum(ypoint**2)
delta1=np.sum((FunctionValues-FunctionAtZero)*xpoint)
delta2=np.sum((FunctionValues-FunctionAtZero)*ypoint)

a = -(delta1*alpha3-delta2*alpha2)/(alpha2**2-alpha1*alpha3) #d/dx
b = -(delta2*alpha1-delta1*alpha2)/(alpha2**2-alpha1*alpha3) #d/dy

print("Wartosci pochodnych bez skalowania")
print("%.4f, %.4f" %(a,b))
print("Roznice")
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])))


print("Wartosci pochodnych bez skalowania", file=output)
print("%.4f, %.4f" %(a,b), file=output)
print("Roznice", file=output)
print("%.4f, %.4f" %(abs(a-fprime(0,0)[0]), abs(b-fprime(0,0)[1])), file=output)

plt.xlim(-maxradius, maxradius)
plt.ylim(-maxradius, maxradius)
plt.grid()
plt.savefig("plot.png")
plt.show()
plt.clf()



###wersja z wagami
weights=Weight(xpoint,ypoint)

Xterm = np.sum(xpoint**2*weights)
XYterm = np.sum(xpoint*ypoint*weights)
Yterm = np.sum(ypoint**2*weights)
matrix = np.array([[Xterm, XYterm], [XYterm, Yterm]])
inverse = np.linalg.inv(matrix)
functionVector = np.array([[np.sum(xpoint*weights*(FunctionValues-FunctionAtZero))], [np.sum(ypoint*weights*(FunctionValues-FunctionAtZero))]]) 
solution=inverse.dot(functionVector)
A, B = solution
print("Wyniki wersji z wagami")
print("%.4f, %.4f" %(A,B))
print("Roznice")
print("%.4f, %.4f" %(abs(solution[0]-fprime(0,0)[0]), abs(solution[1]-fprime(0,0)[1])))


######drugie pochodne

X4Y0term = np.sum(xpoint**4*weights)
X3Y1term = np.sum(xpoint**3*ypoint*weights)
X2Y2term = np.sum(ypoint**2*xpoint**2*weights)
X1Y3term = np.sum(xpoint*ypoint**3*weights)
X0Y4term = np.sum(ypoint**4)
functionVector = np.array([[np.sum(xpoint**2*weights*(FunctionValues-FunctionAtZero-A*xpoint-B*ypoint))],
                           [np.sum(xpoint*ypoint*weights*(FunctionValues-FunctionAtZero-A*xpoint-B*ypoint))],
                           [np.sum(ypoint**2*weights*(FunctionValues-FunctionAtZero-A*xpoint-B*ypoint))]])
matrix = np.array([[X4Y0term/2, X3Y1term, X2Y2term/2], [X3Y1term/2, X2Y2term, X1Y3term/2], [X2Y2term/2, X1Y3term, X0Y4term/2]])
inverse=np.linalg.inv(matrix)
solution=inverse.dot(functionVector)
C,D,E = solution
print("======Drugie pochodne======")
print("Wyniki analityczne")
print("%.4f, %.4f, %.4f" %fbis(0,0))
print("Wyniki wersji z wagami")
print("%.4f, %.4f, %.4f" %(C,D,E))
print("Roznice")
print("%.4f, %.4f, %.4f" %(abs(C-fbis(0,0)[0]), abs(D-fbis(0,0)[1]),abs(E-fbis(0,0)[2])))
