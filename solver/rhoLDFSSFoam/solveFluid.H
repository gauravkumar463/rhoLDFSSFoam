if(!LTS){
	p += runTime.deltaT()*pFlux;

	U += runTime.deltaT()*UFlux;

	e += runTime.deltaT()*eFlux;
}
else{
	p += pFlux/trDeltaT.ref();

	U += UFlux/trDeltaT.ref();

	e += eFlux/trDeltaT.ref();
}
if(stage>0){
	p = ((1-beta[stage])*tp)+(beta[stage]*p);
	U = ((1-beta[stage])*tU)+(beta[stage]*U);
	e = ((1-beta[stage])*te)+(beta[stage]*e);
}

for (direction dir=0; dir<directions.nComponents; dir++)
{
    if (directions[dir] == -1){
       U.replace(dir, 0.0*U.component(dir));
    }
}

p.correctBoundaryConditions();
U.correctBoundaryConditions();
e.correctBoundaryConditions();

thermo.correct();

rho.ref() = psi()*p.ref();
rho.boundaryFieldRef() == psi.boundaryField()*p.boundaryField();

q.ref() = 0.5*magSqr(U.ref()) + k_.ref();
q.boundaryFieldRef() == 0.5*magSqr(U.boundaryField()) + k_.boundaryField();
//q.ref() = 0.5*magSqr(U.ref());
//q.boundaryFieldRef() == 0.5*magSqr(U.boundaryField());
