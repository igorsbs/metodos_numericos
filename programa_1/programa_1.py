import matplotlib.pyplot as plt
import numpy as np


def runge_kutta_4_ordem(func,x,y_init,step):

    h = step
    y = [y_init]

    for i in range(len(x)-1):
        k_1 = func(x[i],y[i])
        k_2 = func(x[i]+h/2,y[i]+k_1*h/2)
        k_3 = func(x[i]+h/2,y[i]+k_2*h/2)
        k_4 = func(x[i]+h,y[i]+k_3*h)

        y_to_append = y[i]+(k_1+2*k_2+2*k_3+k_4)*h/6

        y.append(y_to_append)

    return y

def esf_sed_edo_re0(t,u_z):

    # ni = 0.000891   # Viscosidade da água a 25°C
    ni = 0.001002   # Viscosidade da água a 25°C [Pa.s]
    a = 0.03    # Raio da esfera [m]

    g = 9.81    # Aceleração da gravidade

    rho_s = 8900    # Densidade do Cobre
    rho_l = 1000    # Densidade da Água
    delta_rho = rho_s - rho_l
    
    u_s = 2*delta_rho*g*a**2/(9*ni)

    m_p = 4/3*np.pi*a**3*rho_s

    # stk = m_p*u_s/(6*np.pi*ni*a**2)

    stk = 0.3

    u_z_e = u_z/u_s

    dv = (- u_z + 1)/stk

    # t_e = t*u_s/a

    return dv    

if __name__ == "__main__":

    u_z_e_init = 0
    plt.figure()
    for step in [0.7, 0.5, 0.3, 0.1, 0.07]:
        x = np.arange(0,5+step,step)
        u_z_e = runge_kutta_4_ordem(esf_sed_edo_re0,x,u_z_e_init,step)
        plt.plot(x,u_z_e,label='Step: '+ str(step))
    plt.legend()
    plt.show()