#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 26 12:20:56 2024

@author: leo
"""

from scipy.optimize import least_squares
from scipy import optimize
from scipy.special import jv, jvp, yv, iv
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed   
import os 
 
def Viktorov_matrix(omega, p , cl, ct, R):
    '''
    Equation I.28 can be written as a matrix.  Where the determinant of that matrix is zero, we have the roots
    of the equality.
    omega = frequency
    p = angular wavenumber
    cl = longitudinal wave velocity of material
    ct = shear wave velocity of material
    '''
    # re-scaling for numerical convenience
    R = R*1000
    cl = cl/1000
    ct = ct/1000
    omega = omega/1e6
    # Code starts
    kl = omega/cl
    kt = omega/ct
    x = kl*R
    y = kt*R
    # Matrix form
    a = (jv(p+2,x) + jv(p-2,x) - 2*((kt/kl)**2 -1)*jv(p,x))
    b = (jv(p+2,y) - jv(p-2,y))
    c = (jv(p+2,x) - jv(p-2,x))
    d = (jv(p+2,y) + jv(p-2,y))
    M = np.array([[a, b],[c, d]])
    return M

def Gregory_matrix(omega, p , cl, ct, R):
    '''
    Equation 6.8, p.115 can be written as a matrix.  Where the determinant of that matrix is zero, we have the roots
    of the equality.
    omega = frequency
    p = angular wavenumber
    cl = longitudinal wave velocity of material
    ct = shear wave velocity of material
    '''
    # re-scaling for numerical convenience
    R = R*1000
    cl = cl/1000
    ct = ct/1000
    omega = omega/1e6
    # Code starts
    kl = omega/cl
    kt = omega/ct
    x = kl*R
    y = kt*R
    k = p/R # wavenumber (alpha in Gregory's paper)
    # Matrix form
    a = 2*k**2 - kt**2 - 2*k/R + (2*kl/R)*jv(p+1,x)/jv(p,x)
    b = 4*k**2*(k - (1/R) - kl*jv(p+1,x)/jv(p,x))
    c = k - (1/R) - kt*jv(p+1,y)/jv(p,y)
    d = 2*k**2 - kt**2 - 2*k/R + (2*kt/R)*jv(p+1,y)/jv(p,y)
    M = np.array([[a, b],[c, d]])
    return M

def Sato_matrix62(omega, p, cl, ct, R):
    '''
    Equation 2.1 and 2.2 of Odaka and Usami (1978).
    "Some properties of spheroidal Modes of a Homogeneous Elastic sphere with Special 
    Reference to Radial Dependence of Displacement. Journal of Computational Physics, 29, 431-445, 1978
    omega = frequency
    p = angular wavenumber
    cl = longitudinal wave velocity of material
    ct = shear wave velocity of material
    '''
    # re-scaling for numerical convenience
    R = R*1000
    cl = cl/1000
    ct = ct/1000
    omega = omega/1e6
    # Code starts
    kl = omega/cl
    kt = omega/ct
    x = kl*R
    y = kt*R
    # Matrix form
    a = np.sqrt(2)*np.sqrt(np.pi)*(np.sqrt(2)*np.sqrt(np.pi)*(-p*jv(p, y) + y*jv(p - 1, y) -3*jv(p, y)/2)*np.sqrt(1/y)/R**2 + np.sqrt(2)*np.sqrt(np.pi)*(-2*y*jv(p - 1, y) + (2*p + 1)*jv(p, y) - (2*y**2 -(2*p - 1)*(2*p + 1))*jv(p, y)/4)*np.sqrt(1/y)/R**2)
    b = 2*np.pi*(p - 1/2)*(p + 1/2)*(-p*jv(p, x) + x*jv(p - 1, x) - 3*jv(p, x)/2)
    c = (-p*jv(p, y) + y*jv(p - 1, y) - 3*jv(p, y)/2)*np.sqrt(1/x)*np.sqrt(1/y)/R**4
    d = (-2*x*jv(p - 1, x) + (2*p + 1)*jv(p, x) - (2*y**2 - (2*p - 1)*(2*p + 1))*jv(p, x)/4)*np.sqrt(1/x)/R**2    
    M = np.array([[a, b], [c, d]])
    return M

def Valle_matrix(omega, kh, cL, cT, rho, nu, R):
    '''
    Determinant of coefficients for the 6x6 matrix in the thesis of Christine Valle.
    kh = dimensionless wavenumber (equivalent to kp in Cerv's code)
    b = Outer radius
    a = Inner radius
    cL1 = Longitudinal velocity of medium 1
    cT1 = Transversal velocity of medium 1
    cL2 = Longitudinal velocity of medium 2
    cT2 = Transversal velocity of medium 2
    ''' 
    # SUBROUTINES NEEDED:
    def lame_coefficients(Vp1, Vs1, den1, Vp2, Vs2, den2):
        '''
        Calculate the Lame coefficients for the two media.
        Vp1, Vs1, den1 = Velocities and densities of medium 1.
        Vp2, Vs2, den2 = Velocities and densities of medium 2.    
        Output: Lame coeficients in GPa
        lame_coeffs = [mu1, lambda1, mu2, lambda2]
        '''
        mu1 = den1*Vs1**2
        lambda1 = den1*(Vp1**2 - 2*Vs1**2)
        mu2 = den2*Vs2**2
        lambda2 = den2*(Vp2**2 - 2*Vs2**2)
        lame_coeffs = np.array([mu1, lambda1, mu2, lambda2])
        return lame_coeffs*1e-9
    
    def mu_over_lambda(lame_coeffs,medium_numerator, medium_denominator):
        '''
        Calculate mu_(medium_numerator)/lambda_(medium_denominator) for a double-layered cylinder.
        Input:
        lame_coeffs = array with the lame coefficients as calculated by function 'lame_coefficients'
        lame_coeffs = [mu1, lambda1, mu2, lambda2]
        medium =  1 (outermost layer. Integer)
        medium =  2 (outermost layer. Integer.)
        '''
        if (medium_numerator == 1) & (medium_denominator == 1):
            #print('mu1/lambda1: ')
            val = lame_coeffs[medium_numerator-1]/lame_coeffs[medium_denominator]
        elif (medium_numerator == 1) & (medium_denominator == 2):
            #print('mu1/lambda2: ')
            val = lame_coeffs[medium_numerator-1]/lame_coeffs[medium_denominator+1]
        elif (medium_numerator == 2) & (medium_denominator == 1):
            #print('mu2/lambda1: ')
            val = lame_coeffs[medium_numerator]/lame_coeffs[medium_denominator]           
        return val
    
    def lambda_over_lambda(lame_coeffs,medium_numerator = 2, medium_denominator = 1):
        '''
        Calculate mu_(medium_numerator)/lambda_(medium_denominator) for a double-layered cylinder.
        Input:
        lame_coeffs = array with the lame coefficients as calculated by function 'lame_coefficients'
        lame_coeffs = [mu1, lambda1, mu2, lambda2]
        medium =  1 (outermost layer. Integer)
        medium =  2 (outermost layer. Integer.)
        '''
        val = lame_coeffs[medium_numerator+1]/lame_coeffs[medium_denominator]
        return val
        
    # CODE BEGINS:
    cL1 = cL[0]
    cL2 = cL[1]
    cT1 = cT[0]
    cT2 = cT[1]
    rho1 = rho[0]
    rho2 = rho[1]
    # Calculate the Lame coefficients (lambda and mu)
    lame_coeffs = lame_coefficients(cL1, cT1, rho1, cL2, cT2, rho2)
    # re-scaling for numerical convenience
    R = R*1000
    cL1 = cL1/1000
    cL2 = cL2/1000
    cT1 = cT1/1000
    cT2 = cT2/1000
    omega = omega/1e6
    #omega_m1 = omega*h/cT1
    #k_m = k*h
    h = R*(1 - nu)
    #nu = a/b
    omega_h1 = omega*h/(cT1*(1-nu)) # omega*b/cT1 #  
    omega_h2 = omega*h/(cT2*(1-nu)) # omega*b/cT2 #  
    kappa1 = cL1/cT1
    kappa2 = cL2/cT2        
    # The coefficients of the matrix
    # 1st row ################################################################################################################
    d11 = (2*mu_over_lambda(lame_coeffs,1, 1)*kh*(kh-1) - (1 + 2*mu_over_lambda(lame_coeffs,1, 1))*(omega_h1/kappa1)**2)*jv(kh,omega_h1/kappa1) + 2*(omega_h1/kappa1)*mu_over_lambda(lame_coeffs,1, 1)*jv(kh+1,omega_h1/kappa1)
    d12 = (2*mu_over_lambda(lame_coeffs,1, 1)*kh*(kh-1) - (1 + 2*mu_over_lambda(lame_coeffs,1, 1))*(kh/kappa1)**2)*yv(kh,omega_h1/kappa1) + 2*(omega_h1/kappa1)*mu_over_lambda(lame_coeffs,1, 1)*yv(kh+1,omega_h1/kappa1)
    d13 = 2*mu_over_lambda(lame_coeffs,1, 1)*1j*kh*((kh-1)*jv(kh,omega_h1) - omega_h1*jv(kh+1,omega_h1))
    d14 = 2*mu_over_lambda(lame_coeffs,1, 1)*1j*kh*((kh-1)*yv(kh,omega_h1) - omega_h1*yv(kh+1,omega_h1))
    d15 = 0.0; d16 = 0.0;
    # 2nd row ################################################################################################################
    d21 = 2*1j*kh*((kh - 1)*jv(kh,omega_h1/kappa1) - (omega_h1/kappa1)*jv(kh+1,omega_h1/kappa1))
    d22 = 2*1j*kh*((kh - 1)*yv(kh,omega_h1/kappa1) - (omega_h1/kappa1)*yv(kh+1,omega_h1/kappa1))
    d23 = (omega_h1**2 - 2*kh*(kh - 1))*jv(kh,omega_h1) - 2*omega_h1*jv(kh+1,omega_h1)
    d24 = (omega_h1**2 - 2*kh*(kh - 1))*yv(kh,omega_h1) - 2*omega_h1*yv(kh+1,omega_h1)
    d25 = 0.0; d26 = 0.0;
    # 3rd row ################################################################################################################
    d31 = (2*(mu_over_lambda(lame_coeffs,1, 1)/(nu**2))*kh*(kh - 1) - (1 + 2*mu_over_lambda(lame_coeffs,1, 1))*(omega_h1/kappa1)**2)*jv(kh,omega_h1*nu/kappa1) + 2*(omega_h1/(kappa1*nu))*mu_over_lambda(lame_coeffs,1, 1)*jv(kh+1,omega_h1*nu/kappa1)
    d32 = (2*(mu_over_lambda(lame_coeffs,1, 1)/(nu**2))*kh*(kh - 1) - (1 + 2*mu_over_lambda(lame_coeffs,1, 1))*(kh/kappa1)**2)*yv(kh,omega_h1*nu/kappa1) + 2*(omega_h1/(kappa1*nu))*mu_over_lambda(lame_coeffs,1, 1)*yv(kh+1,omega_h1*nu/kappa1)
    d33 = 2*(mu_over_lambda(lame_coeffs,1, 1)/nu)*1j*kh*(((kh - 1)/nu)*jv(kh,omega_h1*nu) - omega_h1*jv(kh+1,omega_h1*nu))
    d34 = 2*(mu_over_lambda(lame_coeffs,1, 1)/nu)*1j*kh*(((kh - 1)/nu)*yv(kh,omega_h1*nu) - omega_h1*yv(kh+1,omega_h1*nu))
    d35 = -((2*(mu_over_lambda(lame_coeffs,2, 1)/nu**2)*kh*(kh - 1) - (lambda_over_lambda(lame_coeffs,2,1) + 2*mu_over_lambda(lame_coeffs,2, 1))*(omega_h2/kappa2)**2)*jv(kh,omega_h2*nu/kappa2) + 2*(omega_h2/(kappa2*nu))*mu_over_lambda(lame_coeffs,2, 1)*jv(kh+1,omega_h2*nu/kappa2))
    d36 = -2*(mu_over_lambda(lame_coeffs,2, 1)/nu)*1j*kh*(((kh - 1)/nu)*jv(kh,omega_h2*nu) - omega_h2*jv(kh+1,omega_h2*nu))
    # 4th row ################################################################################################################
    d41 = 2*1j*(kh/nu)*(((kh -1)/nu)*jv(kh,omega_h1*nu/kappa1) - (omega_h1/kappa1)*jv(kh+1,omega_h1*nu/kappa1))
    d42 = 2*1j*(kh/nu)*(((kh -1)/nu)*yv(kh,omega_h1*nu/kappa1) - (omega_h1/kappa1)*yv(kh+1,omega_h1*nu/kappa1))
    d43 = (omega_h1**2 - (kh - 1)*(2*kh/(nu**2)))*jv(kh,omega_h1*nu) - (2*omega_h1/nu)*jv(kh+1,omega_h1*nu)
    d44 = (omega_h1**2 - (kh - 1)*(2*kh/(nu**2)))*yv(kh,omega_h1*nu) - (2*omega_h1/nu)*yv(kh+1,omega_h1*nu)
    d45 = 0.0; d46 = 0.0;
    # 5th row ################################################################################################################
    d51 = 0.0; d52 = 0.0; d53 = 0.0; d54 = 0.0;
    d55 = 2*1j*(kh/nu)*(((kh - 1)/nu)*jv(kh,omega_h2*nu/kappa2) - (omega_h2/kappa2)*jv(kh+1,omega_h2*nu/kappa2))
    d56 = (omega_h2**2 - (kh - 1)*(2*kh/(nu**2)))*jv(kh,omega_h2*nu) - 2*(omega_h2/nu)*jv(kh+1,omega_h2*nu)
    # 6th row ################################################################################################################        
    d61 = kh*jv(kh,omega_h1*nu/kappa1) - (omega_h1*nu/kappa1)*jv(kh+1,omega_h1*nu/kappa1)
    d62 = kh*yv(kh,omega_h1*nu/kappa1) - (omega_h1*nu/kappa1)*yv(kh+1,omega_h1*nu/kappa1)
    d63 = 1j*kh*jv(kh,omega_h1*nu)
    d64 = 1j*kh*yv(kh,omega_h1*nu)
    d65 = (omega_h2*nu/kappa2)*jv(kh+1,omega_h2*nu/kappa2) - kh*jv(kh,omega_h2*nu/kappa2) 
    d66 = -1j*kh*jv(kh,omega_h2*nu) 
    # The dij matrix #########################################################################################################
    dij = np.array([
                    [d11 ,d12, d13, d14, d15, d16], 
                    [d21, d22, d23, d24, d25, d26],
                    [d31, d32, d33, d34, d35, d36],
                    [d41, d42, d43, d44, d45, d46],
                    [d51, d52, d53, d54, d55, d56],
                    [d61, d62, d63, d64, d65, d66],
                   ])
    return dij

def Cooke_Rand_matrix(omega, kh, cL, cT, rho, eta, R):
    '''
    Determinant of coefficients for the 6x6 matrix in the paper:
    Cooke, J. R., & Rand, R. H. (1973). A mathematical study of resonance in intact fruits and vegetables using a 3-media elastic sphere model. Journal of agricultural engineering research, 18(2), 141-157.[doi: https://doi.org/10.1016/0021-8634(73)90023-1]
    kh = dimensionless wavenumber (equivalent to kp in Cerv's code)
    R_2 = Outer radius
    R_1 = Inner radius
    eta = R_1/R_2 (ratio of R_1 to R_2)
    cL1 = Longitudinal velocity of medium 1
    cT1 = Transversal velocity of medium 1
    cL2 = Longitudinal velocity of medium 2
    cT2 = Transversal velocity of medium 2
    ''' 
    # SUROUTINES NEEDED:
    def Mn_derivatives(n, A_n, R_n, bessel_type = 'jn'):
        '''First and second derivative of coefficient M_n 
        according to Cooke and Rand (1973):
        Inputs:
        ------ 
        n = order number. 
        A_n = coeffcient multiplying R_n.
        R_n = radius n.
        bessel_type = Either jn or yn for the spherical bessel function.
        Outputs:
        -------
        M_n = coeffcient M_n.
        M_np = first derivative of coeffcient M_n.
        M_npp = second derivative of coeffcient M_n.    
        '''    
        # First derivative of coeffcient M_n
        if bessel_type == 'jn':
            M_n = np.sqrt(2*np.pi/(A_n*R_n))*jv(n + 1/2, A_n*R_n)/2 
            M_np = np.sqrt(2*np.pi/(A_n*R_n))*(-A_n*R_n*jv(n + 3/2, A_n*R_n) + n*jv(n + 1/2, A_n*R_n))/(2*R_n)   
            # M_np = (n/(A_n*R_n))*spherical_jn(n, A_n*R_n) - spherical_jn(n + 1, A_n*R_n) # Eq-89    
        elif bessel_type == 'yn':
            M_n = np.sqrt(2*np.pi/(A_n*R_n))*yv(n + 1/2, A_n*R_n)/2 
            M_np = A_n*np.sqrt(2*np.pi/(A_n*R_n))*(yv(n - 1/2, A_n*R_n)/2 - yv(n + 3/2, A_n*R_n)/2)/2 - np.sqrt(2*np.pi/(A_n*R_n))*yv(n + 1/2, A_n*R_n)/(4*R_n) 
            # M_np = (n/(A_n*R_n))*spherical_yn(n, A_n*R_n) - spherical_yn(n + 1, A_n*R_n) # Eq-89    
        else:
            print('There is a typo on the bessel_type argument! Check it out!')
        # Second derivative of coeffcient M_n
        M_npp = -(2/R_n)*M_np - (A_n**2 - (n*(n + 1)/R_n**2))*M_n # Eq-90
        return M_n, M_np, M_npp
    
    
    def poisson_ratio(c_l, c_t):
        nu = (c_l**2 - 2*c_t**2)/(2*(c_l**2 - c_t**2))
        return nu
    
    def shear_modulus(c_t, rho):
        mu = rho*c_t**2
        return mu
    # CODE BEGINS:
    cL1 = cL[0]
    cL2 = cL[1]
    cL3 = cL[2]
    cT1 = cT[0]
    cT2 = cT[1]
    cT3 = cT[2]
    rho1 = rho[0]
    rho2 = rho[1]
    rho3 = rho[2]
    # Calculate the mu coefficients 
    mu1 = shear_modulus(cT1, rho1) 
    mu2 = shear_modulus(cT2, rho2)
    mu3 = shear_modulus(cT3, rho3)
    nu1 = poisson_ratio(cL1, cT1)
    nu2 = poisson_ratio(cL2, cT2)
    nu3 = poisson_ratio(cL3, cT3)
    # re-scaling for numerical convenience
    R = R#*1000
    cL1 = cL1#/1000
    cL2 = cL2#/1000
    cL3 = cL3#/1000
    cT1 = cT1#/1000
    cT2 = cT2#/1000
    cT3 = cT3#/1000
    omega = omega#/1e6
    
    ##########################
    n = kh - 0.5 
    #theta = 30*np.pi/180
    #eta = 0.3
    R_2 = 0.99*R
    R_1 = eta*R_2
    R_3 = R #1.00005*R_2
    # For medium 1
    A_1 = omega/cL1
    A_2 = omega/cT1
    #eta = np.cos(theta)
    # For medium 2
    A_3 = omega/cL2
    A_4 = omega/cT2
    
    #spherical_jn(n, z[, derivative])
    # bessel function constants
    M_1, M_1p, M_1pp = Mn_derivatives(n, A_1, R_1, bessel_type = 'jn')
    M_2, M_2p, M_2pp = Mn_derivatives(n, A_2, R_1, bessel_type = 'jn')
    M_3, M_3p, M_3pp = Mn_derivatives(n, A_3, R_1, bessel_type = 'jn')
    M_4, M_4p, M_4pp = Mn_derivatives(n, A_4, R_1, bessel_type = 'jn')
    M_5, M_5p, M_5pp = Mn_derivatives(n, A_3, R_2, bessel_type = 'jn')
    M_6, M_6p, M_6pp = Mn_derivatives(n, A_4, R_2, bessel_type = 'jn')
    M_7, M_7p, M_7pp = Mn_derivatives(n, A_3, R_1, bessel_type = 'yn')
    M_8, M_8p, M_8pp = Mn_derivatives(n, A_4, R_1, bessel_type = 'yn')
    M_9, M_9p, M_9pp = Mn_derivatives(n, A_3, R_2, bessel_type = 'yn')
    M_10, M_10p, M_10pp = Mn_derivatives(n, A_4, R_2, bessel_type = 'yn')
    
    # Components of the 6-by-6 determinant (Eqs.42-87)
    ###############################################################################
    S11 = M_1/(A_3*R_1)
    S12 = (M_2p/A_3) + (M_2/(A_3*R_1))
    S13 = -M_3/(A_3*R_1)
    S14 = -M_7/(A_3*R_1)
    S15 = (-M_4p/A_3) - (M_4/(A_3*R_1))
    S16 = (-M_8p/A_3) - M_8/(A_3*R_1)
    ###############################################################################
    S21 = M_1p/A_3
    S22 = n*(n + 1)*M_2/(A_3*R_1)
    S23 = -M_3p/A_3
    S24 = -M_7p/A_3
    S25 = -n*(n + 1)*M_4/(A_3*R_1)
    S26 = -n*(n + 1)*M_8/(A_3*R_1)
    # New constants H_1 and H_2
    H_1 = (nu1/(1 - 2*nu1))*(A_1/A_3)**2
    H_2 = nu2/(1 - 2*nu2)
    ###############################################################################
    S31 = -H_1*M_1 + M_1pp/(A_3**2)
    S32 = -n*(n + 1)*(M_2/(A_3*R_1)**2 - M_2p/(R_1*A_3**2))
    S33 = M_3*H_2 - M_3pp/(A_3**2)
    S34 = H_2*M_7 - M_7pp/(A_3**2)
    S35 = n*(n + 1)*(M_4/(A_3*R_1)**2 - M_4p/(R_1*A_3**2))
    S36 = n*(n + 1)*(M_8/(A_3*R_1)**2 - M_8p/(R_1*A_3**2))
    # New constant H_3
    H_3 = mu1/mu2
    ###############################################################################
    S41 = -2*H_3*(M_1/(A_3*R_1)**2 - M_1p/(R_1*A_3**2))
    S42 = H_3*(M_2pp/A_3**2 - 2*M_2/(A_3*R_1)**2 + n*(n + 1)*M_2/(A_3*R_1)**2)
    S43 = 2*(M_3/(A_3*R_1)**2 - M_3p/(R_1*A_3**2))
    S44 = 2*(M_7/(A_3*R_1)**2 - M_7p/(R_1*A_3**2))
    S45 = -M_4pp/A_3**2 - (M_4/(R_1*A_3**2)**2)*(n*(n + 1) - 2)
    S46 = -M_8pp/(A_3**2) + (M_8/(A_3*R_1)**2)*(n*(n + 1) - 2)
    ###############################################################################
    # New constants
    r_2 = R_2
    r_3 = R_3
    r_23 = r_2/r_3
    mu_23 = mu2/mu3
    rho_23 = rho2/rho3
    H_4 = 2.0 - ((A_4*R_2)**2)*(mu_23*(1 - nu3)/(2*rho_23*(1 + nu3)))
    H_5 = (A_3*R_2*mu_23*r_23*(1 - nu3))/((1 - r_23)*(1 + nu3))
    ###############################################################################
    S51 = 0.0
    S52 = 0.0
    S53 = (-n*(n + 1)*M_5/(A_3*R_2)) - H_2*H_5*M_5 + H_4*M_5p/A_3 \
        + H_5*M_5pp/(A_3**2)
    S54 = (-n*(n + 1)*M_9/(A_3*R_2)) - H_2*H_5*M_9 + H_4*M_9p/A_3 \
        + H_5*M_9pp/(A_3**2)
    S55 = (n*(n + 1)*M_6p/A_3)*(-1 + H_5/(A_3*R_2)) \
        + (n*(n + 1)*M_6/(A_3*R_2))*(-1 + H_4 - H_5/(A_3*R_2))
    S56 = (n*(n + 1)*M_10p/A_3)*(-1 + H_5/(A_3*R_2)) \
        + (n*(n + 1)*M_10/(A_3*R_2))*(-1 + H_4 - H_5/(A_3*R_2))
    ###############################################################################
    # New Constants
    H_6 = (A_3*R_2*mu_23*r_23*(1 - nu3))/(2*(1 - r_23))
    H_7 = (1 - nu3)*(1 + (mu_23*(A_4*R_2)**2)/(2*rho_23))
    H_8 = (1 + nu3) - H_6/(A_3*R_2)
    ###############################################################################
    S61 = 0.0
    S62 = 0.0
    S63 = -n*(n + 1)*M_5/(A_3*R_2) + H_7*M_5/(A_3*R_2) \
        + 2*H_6*M_5/(A_3*R_2)**2 + (H_8 - H_6/(A_3*R_2))*(M_5p/A_3)
    S64 = -n*(n + 1)*M_9/(A_3*R_2) + H_7*M_9/(A_3*R_2) \
        + 2*H_6*M_9/(A_3*R_2)**2 + (H_8 - H_6/(A_3*R_2))*(M_9p/A_3)
    S65 = (-n*(n + 1)*M_6/(A_3*R_2))*(1 - H_8) + H_7*M_6/(A_3*R_2) \
        + 2*H_6*M_6/(A_3*R_2)**2 + (-n*(n + 1) + H_7)*(M_6p/A_3) \
            - H_6*M_6pp/A_3**2 
    S66 = (-n*(n + 1)*M_10/(A_3*R_2))*(1 - H_8) + H_7*M_10/(A_3*R_2) \
        + (-n*(n + 1) + H_7)*(M_10p/A_3) - H_6*M_10pp/A_3**2

    # Dispersion Matrix
    dij = np.array([
        [S11, S12, S13, S14, S15, S16 ],
        [S21, S22, S23, S24, S25, S26 ],
        [S31, S32, S33, S34, S35, S36 ],
        [S41, S42, S43, S44, S45, S46 ],
        [S51, S52, S53, S54, S55, S56 ],
        [S61, S62, S63, S64, S65, S66 ]
        ])   
    
    return dij