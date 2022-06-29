import matplotlib.pyplot as plt
import numpy as np


def Re(D_h, slope=0.014, rho=997, mu=1.002e-3, g=9.81):

    '''
    :param D_h: hydraulic diameter
    :param slope: slope of the tube
    :param rho: density of the liquid
    :param mu: liquid dynamic viscosity
    :param g: standard acceleration due to gravity
    :return: Reynolds number
    '''

    temp = 2 * g * D_h * slope
    num=rho*np.power(temp,4)*D_h
    den=np.power(0.316,4)*mu
    v=np.power(num/den,1/7)

    re=(rho*v*D_h)/mu

    return re


def mass_flow_no_pipe(f_w_perim, D=100*1e-3, s=0.014, rho=997 ):

    '''
    :param f_w_perim: fraction of wetted perimeter
    :param D: diameter of the outer tube
    :param s: slope of the tube
    :param rho: liquid density
    :return: mass flow kg/m^3 in the case of a tube with no inner pipe
    '''

    alpha = f_w_perim * np.pi
    P = alpha * D
    A = (np.power(D,2)/4)*(alpha-np.sin(2*alpha)/2)
    D_h=4*(A/P)

    f=0.316/(Re(D_h=D_h, slope=s, rho=rho) ** 0.25)
    v=np.sqrt(2*9.81*D_h*s/f)


    m=A*v*rho

    return m

def perimeter_limit_case(D,d):

    '''
    The function takes in input the outer and the inner diameters, and computes the wetted perimeter when the
    depth of the liquid is equal to the diameter of the inner pipe

    :param D: outer diameter
    :param d: inner pipe diameter
    :return: wetted perimeter

    '''

    y=d
    p=np.pi*d
    alpha=np.arccos(1-2*y/D)
    P=alpha*D

    return P+p

def mass_flow_inner_pipe_high_liquid_level(f_w_perim,D,d,s,rho,g=9.81):

    '''
    :param f_w_perim: fraction of wetted perimeter
    :param D: outer diameter
    :param d: inner diameter
    :param s: slope of the tube
    :param rho: density of the liquid
    :param g: standard acceleration due to gravity
    :return: mass flow computed for a tube where the depth of thw liquid is equal to the diameter od the inner pipe
    '''

    p=np.pi*d
    P=f_w_perim*(np.pi*D+np.pi*d)-p
    alpha=P/D
    A=(np.power(D,2)/4)*(alpha-np.sin(2*alpha)/2)
    a=(np.pi/4)*np.power(d,2)
    A_h=A-a
    p_h=P+p
    print(p_h==f_w_perim*(np.pi*D+np.pi*d))
    D_h=4*(A_h/p_h)
    f=0.316/(Re(D_h=D, slope=s, rho=rho) ** 0.25)
    v=np.sqrt(2*g*D_h*s/f)

    m=A_h*v*rho
    return m



if __name__ == '__main__':
    print('1')

