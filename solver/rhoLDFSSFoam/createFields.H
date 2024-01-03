#include "createRDeltaT.H"

Info<< "Reading thermophysical properties\n" << endl;

autoPtr<psiThermo> pThermo
(
    psiThermo::New(mesh)
);
psiThermo& thermo = pThermo();


volScalarField& p = thermo.p(); 

volScalarField& e = thermo.he();

Info<< "Reading field U\n" << endl;
volVectorField U
(
    IOobject
    (
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE
    ),
    mesh
);

volScalarField rho
(
    IOobject
    (
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    thermo.rho()
);

volScalarField q
(
    IOobject
    (
        "q",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    0.5*magSqr(U)
);

volScalarField Ma
(
    IOobject
    (
        "Ma",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mag(U)/sqrt(thermo.gamma()*p/rho)
);

volVectorField vorticity
(
    IOobject
    (
        "vorticity",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::curl(U)
);

volVectorField gradRho
(
    IOobject
    (
        "gradRho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    fvc::grad(rho)
);


surfaceScalarField pos
(
    IOobject
    (
        "pos",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("pos", dimless, 1.0)
);

surfaceScalarField neg
(
    IOobject
    (
        "neg",
        runTime.timeName(),
        mesh
    ),
    mesh,
    dimensionedScalar("neg", dimless, -1.0)
);

surfaceScalarField phi
(
    IOobject
    (
        "phi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(rho*U) & mesh.Sf()
);

surfaceVectorField phiUp
(
    IOobject
    (
        "phiUp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    (phi*linearInterpolate(U))
);

surfaceScalarField phiEp
(
    IOobject
    (
        "phiEp",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    phi*linearInterpolate(e+0.5*magSqr(U))
);

surfaceScalarField amaxSf
(
    IOobject
    (
        "amaxSf",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    linearInterpolate(U) & mesh.Sf()
);


Info<< "Creating turbulence model\n" << endl;
autoPtr<compressible::turbulenceModel> turbulence
(
    compressible::turbulenceModel::New
    (
        rho,
        U,
        phi,
        thermo
    )
);

volScalarField muEff("muEff", turbulence->muEff());
volTensorField tauMC("tauMC", muEff*dev2(Foam::T(fvc::grad(U))));


volScalarField k_
(
    IOobject
    (
        "k",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE
    ),
    0.0*turbulence->k()
);

volScalarField rhoFlux
(
    IOobject
    (
        "rhoFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoFlux", dimensionSet(1,-3,-1,0,0,0,0), 1.0)
);

volVectorField rhoUFlux
(
    IOobject
    (
        "rhoUFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("rhoUFlux", dimensionSet(1,-2,-2,0,0,0,0), vector::one)
);

volScalarField rhoEFlux
(
    IOobject
    (
        "rhoEFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("rhoEFlux", dimensionSet(1,-1,-3,0,0,0,0), 1.0)
);

volScalarField pFlux
(
    IOobject
    (
        "pFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("pFlux", dimensionSet(1,-1,-3,0,0,0,0), 1.0)
);

volVectorField UFlux
(
    IOobject
    (
        "UFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedVector("UFlux", dimensionSet(0,1,-2,0,0,0,0), vector::one)
);

volScalarField eFlux
(
    IOobject
    (
        "eFlux",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("eFlux", dimensionSet(0,2,-3,0,0,0,0), 1.0)
);

//scalarList beta = {0.11, 0.2766, 0.5, 1.0};

label numberRungeKuttaStages = 1;
scalarList beta(numberRungeKuttaStages, 1.0);

if(mesh.solutionDict().found("multiStage"))
{
    // Read multi-stage coeffs dict
    dictionary multiStage = mesh.solutionDict().subDict("multiStage");

    // Read number of multi-stages
    if (multiStage.found("numberRungeKuttaStages"))
    {
         multiStage.lookup("numberRungeKuttaStages") >> numberRungeKuttaStages;
         beta.setSize(numberRungeKuttaStages, 1.0);
    }

    // Read multi-stage coeffs
    if(multiStage.found("RKCoeff") && numberRungeKuttaStages > 1)
    {
        ITstream rkCoeffInput = multiStage.lookup("RKCoeff");

        forAll(beta,betai)
        {
            rkCoeffInput >> beta[betai];
        }
    }
}

surfaceScalarField meshPhi
(
    IOobject
    (
        "meshPhi",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar("meshPhi", dimensionSet(0,3,-1,0,0,0,0), 0.0)
);
#include "createFvOptions.H"

volScalarField tp = p;
volVectorField tU = U;
volScalarField te = e;
