from scipy.optimize import least_squares
from scipy import optimize
from scipy.special import jv, jvp, yv, iv
import numpy as np
import matplotlib.pyplot as plt
from joblib import Parallel, delayed   
import os 

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
        self.num = (len(lb) + 1) // 4 

     
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
                matrix1 = Valle_matrix(i, k, vp, vs, rho, nu, self.R)  
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

        
 
 

 
