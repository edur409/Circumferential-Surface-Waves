# Discrepancy in the circumferential phase velocity of the Rayleigh mode for a sphere of steel

![Clorennec_Gallezot](https://github.com/user-attachments/assets/54fca755-f51e-4450-9547-29021644ef1e)

The image above illustrates the difference between: 1) assuming that an integer number of wavelengths fit around the circumference (Clorennec & Royer 2004), and 2) assuming the number of wavelengths around the circumference is an integer and a half (Gallezot et al. (2020)), as per seismological observations.  In seismology, the fact that the number of wavelengths around the circumference is an integer and a half is known as *Jeans' rule* or *Jeans' equation*.  

## Jeans' formula

The standing waves formed by the interference of seismic surface waves travelling several times around the circumference of the Earth should satisfy:

$$ka = l + \frac{1}{2}$$ 

where $l$ is the mode number, $k = w/c$ , $w = 2 \pi/ T$ , $T$ is the period of free oscillations, $a$ is the radius of the Earth, and $c$ is the phase velocity of the corresponding traveling waves.  This is the asymptotic theory of the Earth's normal modes that relates the phase velocity at a particular eigenfrequency to the colatitudinal mode number associated with that frequency. 

## The polar phase shift in seismology

From Brune et al. (1961) we know that every polar passage of a surface wave around the earth produces a phase shift of $\pi/2$ ($\lambda/4$).  They demonstrated this also over the surface of a sphere of steel. Below is a sketch of the waveforms collected on a steel sphere of $4.46~in$ in diameter (Figure 3 of their paper). 

![Brune_et_al_1961_Fig3_sketch](https://github.com/user-attachments/assets/399f314d-fd67-4c6b-accc-a4f417d9c63a)

> "R1 has not made a polar passage, R2 has made one polar passage and is shifted by $\pi/2$, and R3
has made two polar passages and is shifted by $\pi$ relative to R1. This figure clearly
shows that the phase shift on each polar passage is $\pi/2$."

Brune (1964) reminds us that for free oscillations of the earth, the surface wavelength, $\lambda$, is given by:

$$\lambda = 2\pi a/(n + \frac{1}{2})$$

where $n$ is an integer and $a$ is again the Earth's radius.  He claims (page 2107) that Sato an Matumoto (1954) showed that Jeans' formula is a requirement for standing waves.

Sato and Matumoto (1954) explain in page 252 the meaning of $m$ and $n$ in the Legendre polynomials:

> Here it should be added that $m$ must be an integer, while $n$ can take any arbitrary positive real value. When $n$ is an integer, (3.1) and (3.2) represent a stationary vibration, while if $n$ does not take an integral value, (3.1) and (3.2) represent a vibration in which the nodal lines move with the lapse of time, or a wave.

it is worth noting that (3.1) and (3.2) are the equations for the torsional modes (vibrations of the first class) of a sphere made up of two layers (a mantle and a core).

## The laser ultrasound experimental discrepancy on the polar phase shift

Royer et al. (1988) experimentally showed that instead of a polar shift of $\pi/2$, the polar shift of a spherical surface acoustic wave (SSAW) on its polar passage is $\pi$, as shown in the sketch of their Figure 4 below.  They collected experimental evidence using a sphere of steel of $25~mm$ in diameter.  

> "Figure 4 shows the signal detected at a point about $8.4~mm$ away from the pole. The first part (1) corresponds to the passage of the SSAW ring converging toward the pole, the second one (2) to the passage of the SSAW ring diverging from the pole."

![Royer_et_al_1988_Fig4_sketch](https://github.com/user-attachments/assets/00d1be51-e8e9-43bc-b49e-aefa2cd1baf1)

> As it appears clearly, pulse 2 undergoes a phase shift of 180 degrees on passing through the pole. Such a $\pi$-phase shift is well known in optics for focused spherical waves. Similar results have also been reported for bulk acoustic waves. For surface acoustic waves, an experiment has been carried out with cylindrical waves propagating in a plane. A phase jump of $\pi/2$ was demonstrated.

Hsieh (1993) argues in his doctoral thesis (pages 42 and 43) that the polar shift is also $\pi$ and therefore, to form a surface wave resonance, the circumference of the sphere must be an integral number ($M$) of surface wavelengths given by:

$$ka = M - 1$$

Clearly, the academic literature on this topic seems to be in disagreement. 

The dispersion curve as a function of frequency is shown below.

![Steel_Sphere_Discrepancy_Vph_vs_Frequency](https://github.com/user-attachments/assets/d7061f59-e4c4-4403-9ea2-f03c2ad81ffd)

# References

- Jeans, J. H. (1923). The Propagation of Earthquake Waves. Proceedings of the Royal Society of London. Series A, Containing Papers of a Mathematical and Physical Character, 102(718), 554–574. [https://doi.org/10.1098/rspa.1923.0015](https://doi.org/10.1098/rspa.1923.0015)

- Matumoto, T., & Sato, Y. (1954). On the vibration of an elastic globe with one layer. The vibration of the first class. Bull. Earthq. Res. Inst, 32, 247-258. https://doi.org/10.15083/0000034050 [https://repository.dl.itc.u-tokyo.ac.jp/records/34050](https://repository.dl.itc.u-tokyo.ac.jp/records/34050)
  
- Brune, J. N., Nafe, J. E., & Alsop, L. E. (1961). The polar phase shift of surface waves on a sphere. Bulletin of the Seismological Society of America, 51(2), 247-257. [doi: https://doi.org/10.1785/BSSA0510020247](https://doi.org/10.1785/BSSA0510020247)

- Brune, J. N. (1964). Travel times, body waves, and normal modes of the earth. Bulletin of the Seismological Society of America, 54(6A), 2099-2128. [https://doi.org/10.1785/BSSA05406A2099](https://doi.org/10.1785/BSSA05406A2099)
  
- Sato, Y. and Usami, T. (1962). Basic study on the oscillation of a homogeneous elastic sphere-Part I: Frequency of the Free Oscillations, Geophysical Magazine, 31(15-24).

- D. Royer, E. Dieulesaint, X. Jia, Y. Shui; Optical generation and detection of surface acoustic waves on a sphere. Appl. Phys. Lett. 29 February 1988; 52 (9): 706–708. [doi: https://doi.org/10.1063/1.99353](https://doi.org/10.1063/1.99353)

- Hsieh, C. K. P. (1993). Laser-ultrasound characterization of spherical objects. Stanford University. 

- Clorennec, D., & Royer, D. (2004). Investigation of surface acoustic wave propagation on a sphere using laser ultrasonics. Applied physics letters, 85(12), 2435-2437. [doi: https://doi.org/10.1063/1.1791331](https://doi.org/10.1063/1.1791331)

- Gallezot, M., Treyssede, F., & Abraham, O. (2020). Forced vibrations and wave propagation in multilayered solid spheres using a one-dimensional semi-analytical finite element method. Wave Motion, 96, 102555. [doi: https://doi.org/10.1016/j.wavemoti.2020.102555](https://doi.org/10.1016/j.wavemoti.2020.102555)
