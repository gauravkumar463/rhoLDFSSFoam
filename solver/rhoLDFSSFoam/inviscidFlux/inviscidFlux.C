/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright held by original author
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    inviscidFlux

Author
    Oliver Borm  All rights reserved.

\*---------------------------------------------------------------------------*/

#include "inviscidFlux.H"
#include <iostream>

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //+


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::inviscidFlux::evaluateFlux
(
    scalar& rhoFlux,
    vector& rhoUFlux,
    scalar& rhoEFlux,
    scalar& amaxSf,
    const scalar pLeft,
    const scalar pRight,
    const scalar TLeft,
    const scalar TRight,
    const scalar CvLeft,
    const scalar CvRight,
    const vector ULeft,
    const vector URight,
    const scalar rpsiLeft,
    const scalar rpsiRight,
    const scalar eLeft,  // internal energy
    const scalar eRight,
    const scalar kLeft,  // turbulent kinetic energy
    const scalar kRight,
    const scalar gamma,
    const vector Sf,
    const scalar phi,
    const word scheme
) const
{
    scalar magSf = mag(Sf);
    vector n = Sf/magSf;

    scalar rhoLeft = pLeft/rpsiLeft;
    scalar rhoRight = pRight/rpsiRight;

    scalar aLeft = sqrt(gamma*rpsiLeft);
    scalar aRight = sqrt(gamma*rpsiRight);

    scalar uLeft =  (ULeft & n)  - (phi/magSf); // phi is the "flux vector" == u * Sf
    scalar uRight = (URight & n) - (phi/magSf);

    // M+ == MLeft
    // M- == MRight
    if(scheme == "LDFSS"){
        // 
    	scalar rhoHalf = 0.5*(rhoLeft+rhoRight);
    	scalar aHalf = sqrt(0.5*(sqr(aLeft)+sqr(aRight)));

    	scalar MLeft = uLeft/aHalf;
    	scalar MRight = uRight/aHalf;
    	scalar Mavg = sqrt(0.5*(sqr(MLeft)+sqr(MRight)));

    	//scalar Ur2Sf = sqr(aHalf);
    	scalar Ur2Sf = min(sqr(aHalf),max(0.5*(magSqr(ULeft)+magSqr(URight)),0.09*sqr(aHalf)));
    	//scalar Ur2Sf = sqr(aHalf*min(max(1.0,0.5*(qLeft+qRight)),1.0));
    	scalar Mr2Half = Ur2Sf/sqr(aHalf);

    	scalar fHalf =  sqrt((sqr(1-Mr2Half)*sqr(Mavg))+(4*Mr2Half))/(1+Mr2Half);
    	scalar aTildeHalf = fHalf*aHalf;

    	MLeft = MLeft/fHalf;
    	MRight = MRight/fHalf;
    	Mavg = sqrt(0.5*(sqr(MLeft)+sqr(MRight)));

    	amaxSf = 0.5*aTildeHalf*(1+Mr2Half)*(Mavg+1)*magSf;
    	//amaxSf = aHalf*max(max(Mavg,(Mavg+1.0)),Mavg-1.0)*magSf;
    	//amaxSf = aHalf*(Mavg+1.0)*magSf;

    	scalar alphaLeft = 0.5*(1.0 + sign(MLeft));
    	scalar alphaRight = 0.5*(1.0 - sign(MRight));

    	scalar betaLeft = -max(0.0,1-int(mag(MLeft)));
    	scalar betaRight = -max(0.0,1-int(mag(MRight)));

    	scalar M4Left = 0.25*sqr(MLeft+1) + (1.0/8.0)*sqr(sqr(MLeft)-1);
    	scalar M4Right= -0.25*sqr(MRight-1) - (1.0/8.0)*sqr(sqr(MRight)-1);

    	scalar Mplus = (alphaLeft*(1+betaLeft)*MLeft) - (betaLeft*M4Left);
    	scalar Mminus = (alphaRight*(1+betaRight)*MRight) - (betaRight*M4Right);

    	scalar MHalf = 0.25*betaLeft*betaRight*sqr(sqrt(0.5*(sqr(MLeft)+sqr(MRight)))-1.0);


    	scalar UPlus = aTildeHalf*(Mplus - (MHalf*(1.0-((pLeft-pRight)/(2*rhoLeft*Ur2Sf)))));
    	scalar Uminus = aTildeHalf*(Mminus + (MHalf*(1.0+((pLeft-pRight)/(2*rhoRight*Ur2Sf)))));

    	scalar M5Left = 0.25*sqr(MLeft+1)*(2-MLeft) + (3.0/16.0)*MLeft*sqr(sqr(MLeft)-1);
    	scalar M5Right = 0.25*sqr(MRight-1)*(2+MRight) - (3.0/16.0)*MRight*sqr(sqr(MRight)-1);

    	scalar p1Left = alphaLeft*(1+betaLeft) - betaLeft*M5Left;
    	scalar p1Right = alphaRight*(1+betaRight) - betaRight*M5Right;

    	scalar pHalf = 0.5*(pLeft+pRight) + 0.5*(p1Left-p1Right)*(pLeft-pRight) 
                                + rhoHalf*Ur2Sf*(p1Left+p1Right-1.0);

    	scalar HLeft = eLeft + 0.5*magSqr(ULeft) + kLeft + pHalf/rhoLeft;
    	scalar HRight = eRight + 0.5*magSqr(URight) + kRight + pHalf/rhoRight;

    	rhoFlux  = magSf*(rhoLeft*UPlus + rhoRight*Uminus);

    	rhoUFlux = magSf*(rhoLeft*ULeft*UPlus + rhoRight*URight*Uminus) + pHalf*Sf;

    	rhoEFlux = magSf*(rhoLeft*HLeft*UPlus + rhoRight*HRight*Uminus);
    }
    else if(scheme == "AUSMPWplus"){
        // Hong Kim, Chongam Kim, Oh-Hyun Rho,
        // Methods for the Accurate Computations of Hypersonic Flows: I. AUSMPW+Scheme,
        // Journal of Computational Physics, Volume 174, Issue 1, 2001, Pages 38-80,
        // ISSN 0021-9991, https://doi.org/10.1006/jcph.2001.6873.

        // std::cout << "\n#################################################";
        // std::cout << "\nThanks for viewing my code!";
        // std::cout << "\naLeft: "  << aLeft;
        // std::cout << "\naRight: " << aRight;
        // std::cout << "\npLeft: "  << pLeft;
        // std::cout << "\npRight: " << pRight;
        // std::cout << "\nTLeft: "  << TLeft;
        // std::cout << "\nTRight: " << TRight;
        // std::cout << "\nuLeft: "  << uLeft;
        // std::cout << "\nuRight: " << uRight;
        // std::cout << "\n----------------------------";

        const scalar alpha = 3.0/16.0;
        const scalar beta  = 1.0/8.0;

        scalar vLeftSqr  = magSqr(ULeft)  - sqr(uLeft);
        scalar vRightSqr = magSqr(URight) - sqr(uRight);

        // std::cout << "\nvLeftSqr: " << vLeftSqr;
        // std::cout << "\nvRightSqr: " << vRightSqr;
        // std::cout << "\n----------------------------";


    	scalar HLeft   = pLeft/(rhoLeft*(gamma - 1.)) + 0.5*magSqr(ULeft) + kLeft + pLeft/rhoLeft;
    	scalar HRight  = pRight/(rhoRight*(gamma - 1.)) + 0.5*magSqr(URight) + kRight + pRight/rhoRight;
//      	scalar HnLeft  = CvLeft*TLeft   + 0.5*magSqr(uLeft)  + kLeft  + pLeft/rhoLeft;
        // std::cout << "\nHnLeft: " << HnLeft;
//        scalar HnRight = CvRight*TRight + 0.5*magSqr(uRight) + kRight + pRight/rhoRight;
        // std::cout << "\nHnRight: " << HnRight;
        scalar Hn      = 0.5 * (HLeft+HRight - 0.5*(vLeftSqr + vRightSqr));
        // std::cout << "\nHn: "      << Hn;

        scalar aStar = sqrt(max(0.0, (2.0 *(gamma - 1.0) / (gamma + 1.0)) * Hn));

    	// scalar aStar     = sqrt(0.5*(sqr(aLeft)+sqr(aRight)));
        // std::cout << "\naStar: "    << aStar;

        scalar a12 = (0.5 * (uLeft + uRight)) > 0. ?
                     sqr(aStar) / max(uLeft, aStar) :
                     sqr(aStar) / max(uRight, aStar);

        scalar MLeft   = uLeft/a12;   
        scalar MRight  = uRight/a12;   

        // std::cout << "\na12: "    << a12;
        // std::cout << "\nMLeft: "  << MLeft;
        // std::cout << "\nMRight: " << MRight;
        // std::cout << "\n----------------------------";

        scalar psiLeftp =  mag(MLeft) > 1. ? 
                           0.5*(1 + sign(MLeft)) : 
                           0.25*sqr(MLeft + 1.) * (2. - MLeft) + alpha*MLeft*sqr(sqr(MLeft) - 1.);
        scalar psiRightm = mag(MRight) > 1. ? 
                           0.5*(1 - sign(MRight)) : 
                           0.25*sqr(MRight - 1.) * (2. + MRight) - alpha*MRight*sqr(sqr(MRight) - 1.);

        // std::cout << "\npsiLeftp: " << psiLeftp;
        // std::cout << "\npsiRightp: " << psiRightm;

        scalar Mhat = min(1., 1./a12 * sqrt((magSqr(uLeft) + magSqr(uRight))/2.));
        scalar h = sqr(1. - Mhat);
        scalar fp = sqr(1. - h);

        // std::cout << "\nMhat: " << Mhat;
        // std::cout << "\nh: "  <<  h;
        // std::cout << "\nfp: " << fp;
        
        scalar p12 = (pLeft + pRight)/2. + (psiLeftp - psiRightm)/2. * (pLeft - pRight) + fp * (psiLeftp + psiRightm - 1.) * (pLeft + pRight)/2.;

        // std::cout << "\np12: "    << p12;

        scalar Mplus =  mag(MLeft) > 1. ?
                 0.5*(MLeft + mag(MLeft)) :
                 0.25*sqr(MLeft + 1.) + beta * sqr(sqr(MLeft) - 1.);
        scalar Mminus = mag(MRight) > 1. ?
                 0.5*(MRight - mag(MRight)) :
                 -0.25*sqr(MRight - 1.) - beta * sqr(sqr(MRight) - 1.);

        // std::cout << "\nMplus: "  << Mplus;
        // std::cout << "\nMminus: " << Mminus;

        scalar ps = psiLeftp*pLeft + psiRightm*pRight;
        scalar fL = ps != 0. ? pLeft/ps - 1. : 0.;
        scalar fR = ps != 0. ? pRight/ps - 1. : 0.;
        scalar w  = 1 - pow(min(pLeft/pRight, pRight/pLeft), 3.);

        scalar m12 = Mplus + Mminus;
        // std::cout << "\nps: "  << ps;
        // std::cout << "\nM12: " << m12;
        // std::cout << "\n----------------------------";

        scalar MLeftBar  = m12 >= 0. ? 
                           Mplus + Mminus - Mminus * w * (1. + fR) + (fL*Mplus + fR*Mminus) :
                           Mplus * w * (1. + fL);
        scalar MRightBar = m12 >= 0. ? 
                           Mminus * w * (1. + fR) :
                           Mplus + Mminus - Mplus * w * (1. + fL) + (fL*Mplus + fR*Mminus);

        // std::cout << "\nMLeftBar: "  << MLeftBar;
        // std::cout << "\nMRightBar: " << MRightBar;

    	// scalar HLeft  = eLeft + 0.5*magSqr(ULeft) + kLeft + p12/rhoLeft;
    	// scalar HRight = eRight + 0.5*magSqr(URight) + kRight + p12/rhoRight;
    	// scalar HLeft  = CvLeft*TLeft + 0.5*magSqr(ULeft) + kLeft + p12/rhoLeft;
    	// scalar HRight = CvRight*TRight + 0.5*magSqr(URight) + kRight + p12/rhoRight;

    	rhoFlux  = (MLeftBar * rhoLeft        + MRightBar * rhoRight       ) * magSf * a12;
    	rhoUFlux = (MLeftBar * rhoLeft*ULeft  + MRightBar * rhoRight*URight) * magSf * a12 + p12*Sf;
    	rhoEFlux = (MLeftBar * rhoLeft*HLeft  + MRightBar * rhoRight*HRight) * magSf * a12;
        // std::cout << "\n############ NEXT-ITERATION #############\n";
    }
    else if((scheme == "Tadmor") || (scheme == "Kurganov")){
        scalar ap = max(max(uLeft+aLeft,uRight+aRight),0.0);
        scalar am = min(min(uLeft-aLeft,uRight-aRight),0.0);

        scalar aPos = ap/(ap-am);

        amaxSf = magSf*max(mag(am),mag(ap));

        scalar aSf = am*aPos;

        if(scheme == "Tadmor"){
            aSf = -0.5*amaxSf/magSf;
            aPos = 0.5;
        }

        scalar aNeg = 1.0-aPos;

        uLeft *= aPos;
        uRight *= aNeg;

        scalar auLeft = uLeft - aSf;
        scalar auRight = uRight + aSf;

        amaxSf = magSf*max(mag(auLeft),mag(auRight));

	//scalar HLeft = (rpsiLeft/(gamma-1)) + 0.5*magSqr(ULeft) + pLeft/rhoLeft;
        //scalar HRight = (rpsiRight/(gamma-1)) + 0.5*magSqr(URight) + pRight/rhoRight;

        scalar HLeft = eLeft + 0.5*magSqr(ULeft) + kLeft + pLeft/rhoLeft;
        scalar HRight = eRight + 0.5*magSqr(URight) + kRight + pRight/rhoRight;

        rhoFlux = magSf*(rhoLeft*auLeft+ rhoRight*auRight);

        rhoUFlux = magSf*(rhoLeft*ULeft*auLeft + rhoRight*URight*auRight) + (aPos*pLeft+aNeg*pRight)*Sf;

        rhoEFlux = magSf*(rhoLeft*HLeft*auLeft + rhoRight*HRight*auRight + aSf*(pLeft-pRight)) + phi*(aPos*pLeft + aNeg*pRight);

    //    Uf = aPos*ULeft + aNeg*URight;
    }
}

// ************************************************************************* //
