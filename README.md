# rhoLDFSSFoam
A fast and accurate density-based compressible flow solver in OpenFOAM.

There are three advantages of using this solver over rhoCentralFoam available in OpenFOAM.
1. It transforms discretized governing equations in conservative variable to premitive variables (p, $\vec{u}$, T) exactly which eliminates spurious error from decoupling of equations and improves solution stability. All simulations can be performed at CFL 1 which is a challenege in rhoCentralFoam.
2. It uses $3^{rd}$ order low-storage TVD Runge-Kutta time integration scheme which makes it more accurate and stable than rhoCentralFoam.
   
   Kumar, G. & De, A. (2022). An Improved Density-Based Compressible Flow Solver in OpenFOAM for Unsteady Flow Calculations. In Advances in Fluid Mechanics: Modelling and Simulations (pp. 43-66). Singapore: Springer Nature Singapore. (https://link.springer.com/chapter/10.1007/978-981-19-1438-6_2)
   
3. Low-Diffusion Flux splitting scheme of "Edwards, J.R. & Liou, M.-S. “Low-Diffusion Flux-Splitting Methods for Flows at All Speeds,” AIAA Journal, Vol. 36, No. 9, 1998, pp. 1610-1617 (https://doi.org/10.2514/2.587)" is implemented. As an option one can choose fluxScheme: 1) LDFSS, 2) Kurganov, or 3) Tadmor in the fvScheme file to use one of the three schemes.

The comparison of numerical diffusion between Kurganov and LDFSS schemes is made by spanwise vorticity comparison from numerical solution as shown bellow. LDFSS has significantly less diffusion than Kurganov scheme. 

![forwardStepComp1](https://github.com/gauravkumar463/rhoLDFSSFoam/assets/4538589/14f9ae43-8a2b-49da-b9a4-b023ce84f255)







### Please cite this code using the following reference:
Kumar, G. & De, A. (2022). An Improved Density-Based Compressible Flow Solver in OpenFOAM for Unsteady Flow Calculations. In Advances in Fluid Mechanics: Modelling and Simulations (pp. 43-66). Singapore: Springer Nature Singapore.
