## Article: Improvemente of the Cardiac Oscillator Based Model for the Simulation of Bundle #  Branch Blocks
#
#  @author W.J. Felipe
## 

from scipy import linspace
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def SA(t,y):
    x1,y1=y[0],y[1]
    dydt = [y1,
            -a1*y1*(x1-u11)*(x1-u12) - f1*x1*(x1 + d1)*(x1 + e1)
           ]
    return dydt

def AV(t,y):
    x1,y1,x2,y2=y[0],y[1],y[2],y[3]
    dydt = [y1, 
            -a1*y1*(x1-u11)*(x1-u12) - f1*x1*(x1 + d1)*(x1 + e1),
            y2,
            -a2*y2*(x2-u21)*(x2-u22) - f2*x2*(x2+d2)*(x2+e2) + K_SA_AV*(x1-x2)
           ]
    return dydt

def RB(t,y):
    x1,y1,x2,y2,x3RB,y3RB =y[0],y[1],y[2],y[3], y[4], y[5]
    dydt = [y1, 
            -a1*y1*(x1-u11)*(x1-u12) - f1*x1*(x1 + d1)*(x1 + e1),
            y2,
            -a2*y2*(x2-u21)*(x2-u22) - f2*x2*(x2+d2)*(x2+e2) + K_SA_AV*(x1-x2),
            y3RB,
            -a3*y3RB*(x3RB-u31)*(x3RB-u32) - f3*x3RB*(x3RB + d3)*(x3RB + e3) + K_AV_RB*(x2-x3RB)
           ]
    return dydt

def LB(t,y):
    x1,y1,x2,y2,x3LB,y3LB =y[0],y[1],y[2],y[3], y[4], y[5]
    dydt = [y1, 
            -a1*y1*(x1-u11)*(x1-u12) - f1*x1*(x1 + d1)*(x1 + e1),
            y2,
            -a2*y2*(x2-u21)*(x2-u22) - f2*x2*(x2+d2)*(x2+e2) + K_SA_AV*(x1-x2),
            y3LB,
            -a3*y3LB*(x3LB-u31)*(x3LB-u32) - f3*x3LB*(x3LB + d3)*(x3LB + e3) + K_AV_LB*(x2-x3LB)
           ]
    return dydt

#PARAMETROS PARA LA SIMULACIÓN
a1, u11, u12, f1, d1, e1 = 40, 0.83, -0.83, 25,  3, 3.5
a2, u21, u22, f2, d2, e2 = 50, 0.83, -0.83, 8.4, 3, 5
a3, u31, u32, f3, d3, e3 = 50, 0.83, -0.83, 1.5, 3, 12
K_SA_AV, K_AV_RB, K_AV_LB= 100, 285, 285



#################### RETRATO DE FASE #######################
a, b = 0, 2
t= linspace(a,b,1000)
initial_points_SA=[[1,0],[2,2.5],[0,-0.01]]
initial_points_AV=[[1,0,1,0],[2,-2,1,-40],[2,2.5,0,-0.01]]
initial_points_RB=[[1,0,1,0,1,-20],[2,-2,1,-40,-1.5,-40],[2,2.5,0,-0.01,0.001,-0.001]]
initial_points_LB=[[1,0,1,0,1,-40],[2,-2,1,-40,-0.5,-40],[2,2.5,0,-0.01,0.5,-0.001]]

# SA NODE  
for initial_point_SA in initial_points_SA:
    SA_sol = solve_ivp(SA,[a, b],initial_point_SA, t_eval=t)
    plt.plot(SA_sol.y[0],SA_sol.y[1], "-")
plt.grid()
#plt.legend([f"$p_0={initial_point_SA[0],initial_point_SA[1]}$" for initial_point_SA in initial_points_SA])
plt.xlabel("$x_1$")
plt.ylabel("$y_1$")
plt.show()

# AV NODE 
for initial_point_AV in initial_points_AV:
    AV_sol = solve_ivp(AV,[a, b],initial_point_AV, t_eval=t)
    plt.plot(AV_sol.y[2],AV_sol.y[3], "-")
plt.grid()
#plt.legend([f"$p_0={initial_point_AV[2],initial_point_AV[3]}$" for initial_point_AV in initial_points_AV])
plt.xlabel("$x_2$")
#plt.xlabel("$x_2$ : potencial de acción del nodo AV. \n SA    $\overrightarrow{K_{SA-AV}}$    AV")
plt.ylabel("$y_2$")
plt.show()

# RB NODE 
for initial_point_RB in initial_points_RB:
    RB_sol = solve_ivp(RB,[a, b],initial_point_RB, t_eval=t)
    plt.plot(RB_sol.y[4],RB_sol.y[5], "-")
plt.grid()
#plt.legend([f"$p_0={initial_point_RB[4],initial_point_RB[5]}$" for initial_point_RB in initial_points_RB])
plt.xlabel("$x_3$")
#plt.xlabel("$x_2$ : potencial de acción del nodo AV. \n SA    $\overrightarrow{K_{SA-AV}}$    AV")
plt.ylabel("$y_3$")
plt.show()

# LB NODE 
for initial_point_LB in initial_points_LB:
    LB_sol = solve_ivp(LB,[a, b],initial_point_LB, t_eval=t)
    plt.plot(LB_sol.y[4],LB_sol.y[5], "-")
plt.grid()
#plt.legend([f"$p_0={initial_point_LB[4],initial_point_LB[5]}$" for initial_point_LB in initial_points_LB])
plt.xlabel("$x_4$")
#plt.xlabel("$x_2$ : potencial de acción del nodo AV. \n SA    $\overrightarrow{K_{SA-AV}}$    AV")
plt.ylabel("$y_4$")
plt.show()

##### TIME VS x_i ###
a, b=0,16
t= linspace(a,b,1000)

# SA NODE  
SA_sol = solve_ivp(SA,[a, b],[0.001,0.001], t_eval=t)
plt.plot(t,SA_sol.y[0], "-")
plt.xlim([4,16])
plt.ylim([-2,1.75])
plt.ylabel("$\mathbf{SA}$\n $x_1$")
plt.grid()
plt.axes().set_aspect(0.75)
plt.show()

# AV NODE  

AV_sol = solve_ivp(AV,[a, b],[0.001,0.001,0.001,0.001], t_eval=t)
plt.plot(t,AV_sol.y[2], "-")
plt.xlim([4,16])
plt.ylim([-2,1.75])
plt.ylabel("$\mathbf{AV}$\n $x_2$")
plt.grid()
plt.axes().set_aspect(0.75)
plt.show()

# RIGHT BUNDLE BRANCH RB
RB_sol = solve_ivp(RB,[a, b],[0.001,0.001,0.001,0.001,0.001,0.001], t_eval=t)
plt.plot(t,RB_sol.y[4], "-")
plt.xlim([4,16])
plt.ylim([-2,2])
plt.ylabel("$\mathbf{RB}$\n $x_3$")
plt.grid()
plt.axes().set_aspect(0.75)
plt.show()

# LEFT BUNDLE BRANCH LB
LB_sol = solve_ivp(LB,[a, b],[0.001,0.001,0.001,0.001,0.001,0.001], t_eval=t)
plt.plot(t,LB_sol.y[4], "-")
plt.xlim([4,16])
plt.ylim([-2,2])
plt.ylabel("$\mathbf{LB}$\n $x_4$")
plt.grid()
plt.axes().set_aspect(0.75)
plt.show()

