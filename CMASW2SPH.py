from scipy.optimize import least_squares
from scipy import optimize
from scipy.special import jv, jvp, yv, iv
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed   
import os 

from scipy.special import spherical_jn, spherical_yn

def Cooke_Rand_matrix(omega, kh, cL, cT, rho, eta, R):
    '''
    Determinant of coefficients for the 6x6 matrix in the thesis of Christine Valle.
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
    n = kh - 0.5 #1
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

class cmasw2():
    def __init__(self, filename, xinital,  lb, ub, R, Print_r=True, Savefile= True):
        self.filename = filename
        self.R = R
        self.lb = lb
        self.ub = ub
        self.init = xinital 
        self.Print_r = Print_r
        self.Savefile = Savefile
   
        self.measure = np.loadtxt(filename) 
        self.omega = self.measure[:,0]*2*np.pi
        self.y_data =self.measure[:,1] 
        self.loss = []
        self.para = []
        self.jac = [] # save Jacovian matrix (added on 21/08/2024
        self.num = 3 #(len(lb) + 1) // 4 

     
    def aimfuc(self, x, omega, y_data):
        '''
        Calculate the loss function. Perform the forward model and save 
        outputs of a fitting iteration.
        '''
        vp = np.array(x[:self.num])
        vs = np.array(x[self.num:2*self.num])  
        rho = np.array(x[2 * self.num:3 * self.num]) 
        nu = np.array(x[3 * self.num])           
        yt =  self.velocity(omega, vp, vs, rho, nu) 
        loss = np.abs(yt - y_data)    
        
        if self.Print_r: 
            plt.clf()  
            plt.scatter(omega/(2*np.pi), y_data, label= 'measurement')    
            plt.plot(omega/(2*np.pi), yt, ':', label='Prediction')    
            self.loss.append(np.sqrt(np.sum(loss**2))) 
            plt.title('TRF Algorithm, Iteration Number: %s; L2 Norm: %s' %(len(self.loss), int(self.loss[-1])))
            plt.legend() 
            plt.xlabel('frequency (Hz)') 
            plt.ylabel('Phase Velocity (m/s)')
            plt.savefig('dispersion_update.png')  
            print(np.sqrt(np.sum(loss ** 2)))
        
        if self.Savefile: 
            np.savetxt('prediction.txt', yt)
            self.para.append(x)
        return loss  


    def inverse(self):
        # I added an ftol for the least_squares algorithm below!
        popt = least_squares(self.aimfuc, self.init, bounds=(self.lb, self.ub), method='trf', args=(self.omega, self.y_data), ftol = 0.01) # 0.5 for example
        print("inversion is done") 
        r = popt.x   
        jac = popt.jac # Added 21/08/2024
        self.jac.append(jac) # Added 21/08/2024
        Vp = np.array(r[:self.num])
        Vs = np.array(r[self.num:self.num*2])
        Rho =  np.array(r[self.num*2:self.num*3]) 
        Nu = np.array(r[self.num * 3])
        #E = 2*rho*Vs**2*2*(1+mu)         
        
        print('Shear Velocity', Vs) 
        #print('Poisson Ratio', mu)  
        print('P wave velocity:', Vp)
        #print('Youngs Modulus (MPa):', E/10**6) 
        print('Density',Rho)
        print('ra/rb', Nu)  

        if self.Savefile:
            np.savetxt('loss.csv', self.loss)
            np.savetxt('parameters.csv', self.para)
            
    def velocity(self, omega, vp, vs, rho, nu):
        '''
        Function to perform the forward model.
        '''
        def opt(i):
            def fun(k): 
                matrix1 = Cooke_Rand_matrix(i, k, vp, vs, rho, nu, self.R)  
                matrix2 = (matrix1)#/10**11                
                sign, logdet = np.linalg.slogdet(matrix2)
                return np.real(sign * np.exp(logdet))
            c_test = 0.9*np.min(self.y_data) # 90% of the min. vel. of the dispersion curve is a sensible starting choice
            incre = self.R*i / c_test
            root = 0.00001
            for j in range(10**6):
                past = incre
                val1 = fun(incre)
                incre =  incre - 0.01
                val2 = fun(incre)
                if (np.real(val1) * np.real(val2) <= 0):
                    root =  optimize.brentq(fun,incre, past)      
                    break 
            return (self.R*i/root)     #give one value at a frequency
    
        def final(n):
            y = opt(omega[n])
            return y   
        
        z2 = Parallel(n_jobs=-1)(delayed(final)(i) for i in range(len(omega))) 
        return z2

        
 
 

 
