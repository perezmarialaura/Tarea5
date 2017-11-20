import numpy as np
import matplotlib.pyplot as plt

Bb = 0.2497
Bd = 5.16
Ad = 0.3105
Ah = 64.3
def V_obs(R,Mb,Md,Mh):
    return Mb**0.5*R/(R**2 + Bb**2)**(3/4) + Md**0.5*R/(R**2 + (Bd+Ad)**2)**(3/4) + Mh**0.5/(R**2 + Ah**2)**0.25

data = np.genfromtxt("RadialVelocities.dat", delimiter= " ", skip_header = 1)
R = data[:, 0]
V = data[:, 1]
mb, md, mh = np.genfromtxt("M.dat", delimiter=" ", unpack = True)
plt.figure()
plt.title("Curvas de rotación de una galaxia")
plt.scatter(R,V, label="Datos")
plt.plot(R,V_obs(R, mb,md,mh), label="Fit bayesiano: Mb = "+str(mb)+", Md = "+str(md)+", Mh = "+str(mh),linewidth=2.0)
plt.legend(prop={'size': 10})
plt.xlabel("Distancia al centro de la galaxia (kpc)")
plt.ylabel("Velocidad de rotación (km/s)")
plt.savefig("galaxias.png")
