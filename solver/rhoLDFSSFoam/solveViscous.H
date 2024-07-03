solve
(
	fvm::ddt(rho, U) - fvc::ddt(rho, U)
  - fvm::laplacian(muEff, U)
  - fvc::div(tauMC)
);
U.correctBoundaryConditions();

e.ref() -= 0.5*magSqr(U.ref()) + k_.ref() - q.ref();
//e.ref() -= 0.5*magSqr(U.ref()) - q.ref();
e.correctBoundaryConditions();

solve
(
	fvm::ddt(rho,e) - fvc::ddt(rho, e)
  - fvm::laplacian(thermo.gamma()*(thermo.alpha()+turbulence->alphat()), e)
  - fvc::div
    (
        (fvc::interpolate(muEff)*mesh.magSf()*fvc::snGrad(U)
          + fvc::dotInterpolate(mesh.Sf(), tauMC)
        ) & linearInterpolate(U)
    )
);

e.correctBoundaryConditions();

thermo.correct(); 

p.ref() = rho.ref()/psi();
p.correctBoundaryConditions();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();
