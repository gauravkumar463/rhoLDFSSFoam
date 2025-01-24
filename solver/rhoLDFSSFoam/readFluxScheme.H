word fluxScheme("LDFSS");
if (mesh.schemesDict().readIfPresent("fluxScheme", fluxScheme))
{}

if ((fluxScheme == "Tadmor") || (fluxScheme == "Kurganov") || (fluxScheme == "LDFSS") || (fluxScheme == "AUSMPWplus"))
{
    Info<< "fluxScheme: " << fluxScheme << endl;
}
else
{
    FatalErrorInFunction
        << "fluxScheme: " << fluxScheme
        << " is not a valid choice. "
        << "Options are: AUSMPWplus, LDFSS, Tadmor, Kurganov"
        << abort(FatalError);
}


bool inviscid = false;
if (mesh.schemesDict().readIfPresent("inviscid", inviscid))
{
    if (inviscid == true)
    {
        Info<< "Running solver in inviscid mode" << endl;
    }
}

dimensionedScalar e0("e0",dimVelocity*dimVelocity,0.0);
mesh.schemesDict().readIfPresent("e0", e0);

/*dimensionedVector acceleration("accelration",dimVelocity/dimTime,vector::zero);
mesh.schemesDict().readIfPresent("acceleration", acceleration);
Info << acceleration << endl;

dimensionedScalar T0("T0",dimTime,0.0);
mesh.schemesDict().readIfPresent("T0", T0);
Info << T0 << endl;*/
