/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    rhoCentralFoam

Description
    Density-based compressible flow solver based on central-upwind schemes of
    Kurganov and Tadmor.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "dynamicFvMesh.H"
#include "psiThermo.H"
#include "turbulentFluidThermoModel.H"
#include "fixedRhoFvPatchScalarField.H"
#include "directionInterpolate.H"
#include "localEulerDdtScheme.H"
#include "fvcSmooth.H"
#include "motionSolver.H"
#include "inviscidFlux/inviscidFlux.H"
#include "fvOptions.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #define NO_CONTROL
    #include "postProcess.H"

    #include "setRootCase.H"
    #include "createTime.H"
    #include "createDynamicFvMesh.H"
    #include "readFluxScheme.H"
    #include "createFields.H"
    #include "createFieldRefs.H"
    #include "createDynFields.H"
    #include "createTimeControls.H"

    turbulence->validate();
    thermo.correct();

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    dimensionedScalar v_zero("v_zero", dimVolume/dimTime, 0.0);

    // Courant numbers used to adjust the time-step
    scalarField sumAmaxSf;
    scalar CoNum = 0.0;
    scalar meanCoNum = 0.0;

    Info<< "\nStarting time loop\n" << endl;
    
    while (runTime.run())
    {
        surfaceScalarField k_pos(interpolate(k_, pos, T.name()));
        surfaceScalarField k_neg(interpolate(k_, neg, T.name()));
        
        #include "inviscidFlux.H"
        #include "transformFluxes.H"     
        #include "readTimeControls.H"
        if(!LTS)
        {
            #include "centralCourantNo.H"
            #include "setDeltaT.H"
        }
        else
        {
            #include "setRDeltaT.H"
        }  
        runTime++;
        mesh.update();

        Info<< "Time = " << runTime.timeName() << nl << endl;
        tp = p;
        tU = U;
        te = e; 
        if (mesh.moving())
        {
            meshPhi = mesh.phi();
        }
        volScalarField rv = V.oldTime()/V;
        muEff = turbulence->muEff();
        tauMC = muEff*dev2(Foam::T(fvc::grad(U)));       
        for(label stage=0;stage<beta.size();stage++)
        {
            if(stage>0){
                #include "inviscidFlux.H"
                #include "transformFluxes.H"
            }
            #include "solveDynFluid.H"
        }

        if(!inviscid == true)
        {
            #include "solveViscous.H"
            turbulence->correct();
            k_ = turbulence->k();   
        }
        sensor = mag(fvc::grad(rho));
        Ma = mag(U)/sqrt(thermo.gamma()*p/rho);
        Info << "min/max(T) = " << gMin(T) << "/" << gMax(T) << endl;
        Info << "min/max(p) = " << gMin(p) << "/" << gMax(p) << endl;
        Info << "min/max(rho) = " << gMin(rho) << "/" << gMax(rho) << endl;
        Info << "min/max(Ma) = " << gMin(Ma) << "/" << gMax(Ma) << endl;
        Info << "min/max(k) = " << gMin(k_) << "/" << gMax(k_) << endl;        
        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}

// ************************************************************************* //
