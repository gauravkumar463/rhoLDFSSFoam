//scalar fbodySwitch=Foam::neg(mesh.time()-T0).value();

rhoFlux = - fvc::div(phi);

rhoUFlux = - fvc::div(phiUp);// + rho*fbody*fbodySwitch;

rhoEFlux = - fvc::div(phiEp);// + rho*(U & fbody)*fbodySwitch;

pFlux = ((q-(e-e0)+(thermo.Cv()*T))*rhoFlux - (U & rhoUFlux) + rhoEFlux)*(thermo.gamma()-1.0);

UFlux = (-(rhoFlux*U) + rhoUFlux)/(rho);

eFlux = ((q-(e-e0))*rhoFlux - (U & rhoUFlux) + rhoEFlux)/rho;
