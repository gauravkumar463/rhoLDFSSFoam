// calculate square of speed of sound
volScalarField c2 = (thermo.gamma()-1)*thermo.Cp()*T;

rpsi = 1/psi;

surfaceScalarField rpsi_pos(interpolate(rpsi, pos, T.name()));
surfaceScalarField rpsi_neg(interpolate(rpsi, neg, T.name()));

surfaceScalarField e_pos(interpolate(e, pos, T.name()));
surfaceScalarField e_neg(interpolate(e, neg, T.name()));

surfaceScalarField T_pos(interpolate(T, pos, T.name()));
surfaceScalarField T_neg(interpolate(T, neg, T.name()));

surfaceScalarField p_pos(interpolate(p, pos, p.name()));
surfaceScalarField p_neg(interpolate(p, neg, p.name()));

surfaceVectorField U_pos(interpolate(U, pos, U.name()));
surfaceVectorField U_neg(interpolate(U, neg, U.name()));

surfaceScalarField gammaSf(linearInterpolate(thermo.gamma()));

forAll(phi,faceI){
    InviscidFlux.evaluateFlux
    (
        phi[faceI],
        phiUp[faceI],
        phiEp[faceI],
        amaxSf[faceI],
        p_pos[faceI],
        p_neg[faceI],
        T_pos[faceI],
        T_neg[faceI],
        U_pos[faceI],
        U_neg[faceI],
        rpsi_pos[faceI],
        rpsi_neg[faceI],
        e_pos[faceI]-e0.value(),
        e_neg[faceI]-e0.value(),
        k_pos[faceI],
        k_neg[faceI],
        gammaSf[faceI],
        mesh.Sf()[faceI],
        meshPhi[faceI],
	fluxScheme
    );
}

surfaceScalarField::Boundary& bfield = phi.boundaryFieldRef();

forAll(bfield, patchi)
{
    const scalarField& currPatch = phi.boundaryField()[patchi];

    fvsPatchScalarField& phip = phi.boundaryFieldRef()[patchi];
    fvsPatchVectorField& phiUpp = phiUp.boundaryFieldRef()[patchi];
    fvsPatchScalarField& phiEpp = phiEp.boundaryFieldRef()[patchi];
    fvsPatchScalarField& amaxSfp = amaxSf.boundaryFieldRef()[patchi];
    const fvPatchScalarField& pp = p.boundaryField()[patchi];
    const fvPatchScalarField& Tp = T.boundaryField()[patchi];
    const fvPatchVectorField& Up = U.boundaryField()[patchi];
    const fvPatchScalarField& rpsip = rpsi.boundaryField()[patchi];
    const fvPatchScalarField& ep = e.boundaryField()[patchi];
    const fvPatchScalarField& kp = k_.boundaryField()[patchi];
    const fvsPatchScalarField& gammap = gammaSf.boundaryField()[patchi];
    const fvsPatchVectorField& Sfp = mesh.Sf().boundaryField()[patchi];
    const fvsPatchScalarField& meshPhip = meshPhi.boundaryField()[patchi];
    
    if (bfield[patchi].coupled())
    {
        const scalarField pcp = p.boundaryField()[patchi].patchInternalField();
        const vectorField Ucp = U.boundaryField()[patchi].patchInternalField();
        const scalarField rpsicp = rpsi.boundaryField()[patchi].patchInternalField();
        const scalarField ecp = e.boundaryField()[patchi].patchInternalField();
        const scalarField kcp = k_.boundaryField()[patchi].patchInternalField();
        const scalarField pcnb = p.boundaryField()[patchi].patchNeighbourField();
        const scalarField Tcnb = T.boundaryField()[patchi].patchNeighbourField();
        const vectorField Ucnb = U.boundaryField()[patchi].patchNeighbourField();
        const scalarField rpsicnb = rpsi.boundaryField()[patchi].patchNeighbourField();
        const scalarField ecnb = e.boundaryField()[patchi].patchNeighbourField();
        const scalarField kcnb = k_.boundaryField()[patchi].patchNeighbourField();

        forAll(currPatch,faceI)
        {
            InviscidFlux.evaluateFlux
            (
                phip[faceI],
                phiUpp[faceI],
                phiEpp[faceI],
                amaxSfp[faceI],
                pcp[faceI],
                pcnb[faceI],
                Tcp[faceI],
                Tcnb[faceI],
                Ucp[faceI],
                Ucnb[faceI],
                rpsicp[faceI],
                rpsicnb[faceI],
                ecp[faceI]-e0.value(),
                ecnb[faceI]-e0.value(),
                kcp[faceI],
                kcnb[faceI],
                gammap[faceI],
                Sfp[faceI],
                meshPhip[faceI],
		fluxScheme
            );
        }
    }
    else
    {
        forAll(currPatch,faceI)
        {
            InviscidFlux.evaluateFlux
            (
                phip[faceI],
                phiUpp[faceI],
                phiEpp[faceI],
                amaxSfp[faceI],
                pp[faceI],
                pp[faceI],
                Tp[faceI],
                Tp[faceI],
                Up[faceI],
                Up[faceI],
                rpsip[faceI],
                rpsip[faceI],
                ep[faceI]-e0.value(),
                ep[faceI]-e0.value(),
                kp[faceI],
                kp[faceI],
                gammap[faceI],
                Sfp[faceI],
                meshPhip[faceI],
		fluxScheme
            );
        }

    }
}
