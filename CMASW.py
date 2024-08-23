from scipy.optimize import least_squares
from scipy import optimize
from scipy.special import jv, jvp, yv, iv
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve
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


class cmasw():
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
        vp = np.array(x[0])
        vs = np.array(x[1])
        yt =  self.velocity(omega, vp, vs) 
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
        popt = least_squares(self.aimfuc, self.init, bounds=(self.lb, self.ub), method='trf', args=(self.omega, self.y_data))
        print("inversion is done") 
        r = popt.x
        Vp = np.array(r[0])
        Vs = np.array(r[1])
                
        print('Shear Velocity', Vs)         
        print('P wave velocity:', Vp)        

        if self.Savefile:
            np.savetxt('loss.csv', self.loss)
            np.savetxt('parameters.csv', self.para)
            
    def velocity(self, omega, vp, vs):
        '''
        Function to perform the forward model.
        '''
        def opt(i):
            def fun(k): 
                matrix1 = Gregory_matrix(i, k , vp, vs, self.R) #Viktorov_matrix(i, k , vp, vs, self.R) 
                matrix2 = matrix1 #np.real(matrix1)#/10**11
                sign, logdet = np.linalg.slogdet(matrix2)
                return sign * np.exp(logdet)
            c_test = 0.5*np.min(self.y_data) # 90% of the min. vel. of the dispersion curve is a sensible starting choice
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

        
 
 

 
