volScalarField V
(
    IOobject
    (
        "V",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("V", dimVolume, 1.0)
);
V.ref() = mesh.V();
V.storeOldTime();

volScalarField sensor
(
    IOobject
    (
       "sensor",
       runTime.timeName(),
       mesh,
       IOobject::READ_IF_PRESENT,
       IOobject::AUTO_WRITE
   ),
   mag(fvc::grad(rho))
);
