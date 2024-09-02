# Spheroidal vibrational eigenfrequencies for simple shapes

The spheroidal eigenfrequencies for order $n$ and mode $l$ follow this notation ${}_nS_l$.

- Homogeneous Sphere <a target="_blank" href="https://colab.research.google.com/github/edur409/Circumferential-Surface-Waves/blob/main/Fundamental_Modes/Sphere_Fundamental_Modes.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- Homogeneous Sphere (Correct $n = 0$ mode) <a target="_blank" href="https://colab.research.google.com/github/edur409/Circumferential-Surface-Waves/blob/main/Fundamental_Modes/Sphere_Fundamental_Modes_Corrected.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- Homogeneous Cylinder <a target="_blank" href="https://colab.research.google.com/github/edur409/Circumferential-Surface-Waves/blob/main/Fundamental_Modes/Cylinder_Fundamental_Modes.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

- Cooke and Rand(1973) <a target="_blank" href="https://colab.research.google.com/github/edur409/Circumferential-Surface-Waves/blob/main/Fundamental_Modes/Cooke_Rand_1973_Fundamental_Modes.ipynb">
  <img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open In Colab"/>
</a>

# Comments on the fundamental eigenmodes of a sphere

You can compare the values to the ones of Lucien Saviot [here](https://saviot.cnrs.fr/lamb/index.en.html). In these notebooks, frequency intervals for root searching are tailored for ultrasonic frecuencies on spherical objects of radii on the centimeter scale, you can modify the codes for more general purposes. 

The fundamental modes for order $n = 0$ are a special case covered by Lamb (1881).   Saviot et al. (2004) is a comment to remind us that in a homogeneous sphere:

> For materials with positive Poisson ratio, it is impossible to have the energy for the fundamental $n = 0$ mode smaller than the energy for the fundamental $n = 2$ one.

Lamb (1881) implies the same on the last paragraph of his article in page 211 (italics here for emphasis):

> 13. As an application of the preceeding results we may calculate the frequency of vibration of a steel ball one centimetre in radius, *for the slowest of those fundamental modes in which the surface oscillates in the form of a harmonic spheroid of the second order.*  In $\S$ 12 we obtained for this case $ka/\pi = .842$.  Now, $ka/\pi = T_0/r$, where $r$ is the period, and $T_0 = 2a/\sqrt{\mu \rho^{-1}}$. Making then $a = 1$, and adopting from Everett the values $\mu = 8.19\times10^{11}$, $\rho = 7.85$, in C.G.S. measure, I find that the frequency $r^{-1} = 136000$, about. 

If we use $v_p = 6009$ m/s and $v_s = 3212$ m/s for steel, the value of $ka/\pi$ is the same as that of Lamb (1881), which corresponds to a frequency of $135276$ Hz or about $136000$ Hz as he estimates. In other words, the gravest frequency of vibration for a homogenous elastic sphere with positive Poisson ratio corresponds to the $_2S_0$ mode.  This mode ($_2S_0$) is the change of shape from an oblate to a prolate ellipsoid not the so called breathing mode $_0S_0$ (volumetric change), for which the frequency is remarkably larger in comparison.

## Fundamental eigenfrequencies: the roots of these equations
* For $n = 0$ and $l > 0$:
  
  - From Lamb (1881), page 201, Eq.59:
    
$$tan(\theta) = \frac{4\theta}{4 - \frac{\theta^2}{\lambda^2}}$$ 

where

> $\lambda^2 = \left(\frac{v_s}{v_p}\right)^2$  

> $\theta = \frac{2 \pi f R}{v_p} = \frac{\omega}{v_p} = k_l R$$

  - or equivalently from Aki and Richards:
    
$$cot(x) = \frac{1}{x} - \frac{1}{4}\lambda^2 x$$

* For $n > 0$ and $l > 0$:

$$\begin{aligned} & {\left[\left(\frac{1}{2}-\frac{n(n-1)}{\eta^2}\right) J_{n+1 / 2}(\xi)-\frac{2 \xi}{\eta^2} J_{n+3 / 2}(\xi)\right]\left[\frac{2}{\eta} J_{n+3 / 2}(\eta)+\left(\frac{2\left(n^2-1\right)}{\eta^2}-1\right) J_{n+1 / 2}(\eta)\right]} \\ 
& \quad+2 n(n+1)\left[\frac{n-1}{\eta^2} J_{n+1 / 2}(\xi)-\frac{\xi}{\eta^2} J_{n+3 / 2}(\xi)\right]\left[\frac{n-1}{\eta^2} J_{n+1 / 2}(\eta)-\frac{1}{\eta} J_{n+3 / 2}(\eta)\right]=0\end{aligned}$$

where

> $R$: radius of the sphere

> $f$: frequency

> $\eta = \frac{2 \pi f R}{v_p} = \frac{\omega}{v_p} = k_l R$

> $\xi = \frac{2 \pi f R}{v_s} = \frac{\omega}{v_p} = k_t R$

> $J_{\nu}(z)$: Bessel function of the first kind 

## References:

- Lamb, H. (1881). On the vibrations of an elastic sphere. Proceedings of the London Mathematical Society, 1(1), 189-212. [doi: https://doi.org/10.1112/plms/s1-13.1.189](https://doi.org/10.1112/plms/s1-13.1.189) 

- Sato, Y. and Usami, T. (1962). Basic study on the oscillation of a homogeneous elastic sphere-Part I: Frequency of the Free Oscillations, Geophysical Magazine, 31(15-24).

- Saviot, L., Murray, D. B., Mermet, A., & Duval, E. (2004). Comment on “Estimate of the vibrational frequencies of spherical virus particles”. Physical Review E, 69(2), 023901. [doi: https://doi.org/10.1103/PhysRevE.69.023901](https://doi.org/10.1103/PhysRevE.69.023901)
