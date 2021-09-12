# Solutions to the Van der Pol Equation: a Model of
# Aortic Blood Flow

# Aortic Blood Flow

from scipy import linspace
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def vdp(t, z):
    x, y = z
    return [y, mu*(1 - x**2)*y - 5.6*x]

a, b = 0, 10

#mus = [0, 1, 2]
mus=[1,1,1]
#styles = ["-", "--", ":"]
points=[[1,0],[4,0],[-4,4]]
t = linspace(a, b, 500)

for mu,point  in zip(mus,points):
    sol = solve_ivp(vdp, [a, b], point, t_eval=t)
    plt.plot(sol.y[0], sol.y[1], "-")
  
# make a little extra horizontal room for legend
plt.xlim([-5,5])    
plt.legend([f"$p_0={point}$" for point in points])
plt.axes().set_aspect(0.5)
