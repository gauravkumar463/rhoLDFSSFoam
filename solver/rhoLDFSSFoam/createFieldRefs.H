volScalarField& T = thermo.T();
const volScalarField& psi = thermo.psi();

volScalarField rpsi
(
    IOobject
    (
        "rpsi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    1/psi
);

inviscidFlux InviscidFlux;

const Vector<label>& directions = mesh.solutionD();

/*volVectorField fbody
(
    IOobject
    (
        "fbody",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE
    ),
    mesh,
    acceleration
);*/
