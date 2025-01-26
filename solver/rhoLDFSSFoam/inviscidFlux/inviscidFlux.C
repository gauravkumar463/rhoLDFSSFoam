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

        // Flux scheme constants
        const scalar alpha = 3.0/16.0;
        const scalar beta  = 1.0/8.0;

        // Tangential velocity components, squared
        const scalar vLeftSqr  = magSqr(ULeft)  - sqr(uLeft);
        const scalar vRightSqr = magSqr(URight) - sqr(uRight);

        // "Normal" enthalpy Eq (31 ff.)
    	const scalar HnLeft  = eLeft  + 0.5*magSqr(ULeft)  + kLeft  + pLeft/rhoLeft;
        const scalar HnRight = eRight + 0.5*magSqr(URight) + kRight + pRight/rhoRight;
        const scalar Hn      = 0.5 * (HnLeft+HnRight - 0.5*(vLeftSqr + vRightSqr));
        
        const scalar aStar   = sqrt(2.0*(gamma - 1)/(gamma + 1)*Hn);

        // Eq (32)
        scalar a12;

        if ( (0.5*(uLeft + uRight)) > 0 )
        {
            a12 = sqr(aStar)/max(fabs(uLeft), aStar);
        }
        else
        {
            a12 = sqr(aStar)/max(fabs(uRight), aStar);
        }

        // Eq (29)
        const scalar MLeft   = uLeft/a12;   
        const scalar MRight  = uRight/a12;   

        // Mach number and pressure splitting functions Eq (27) and (28) ("Mach_Left_plus, ...")
        // Eq (27)
        scalar Mlp, Plp;

        if (MLeft >= 1.0)
        {
            // scalar Mlp = 0.5*(MLeft  + mag(MLeft));
            Mlp = MLeft;
            Plp = 1.0;
        }
        else if (MLeft > -1.0)
        {
            Mlp = sqr(1 + MLeft)/4 + beta*sqr(1 - sqr(MLeft));
            Plp = sqr(1 + MLeft)*(2 - MLeft)/4 + alpha * MLeft*sqr(1 - sqr(MLeft));
        }
        else
        {
            Mlp = 0.0;
            Plp = 0.0;
        }

        scalar Mrm, Prm;

        if (MRight <= -1.0)
        {
            Mrm = MRight;
            Prm = 1.0;
        }
        else if (MRight < 1.0)
        {
            Mrm = - sqr(1 - MRight)/4 - beta * sqr(1 - sqr(MRight));
            Prm = sqr(1 - MRight)*(2 + MRight)/4 - alpha * MRight * sqr(1 - sqr(MRight));
        }
        else
        {
            Mrm = 0.0;
            Prm = 0.0;
        }
          
        // scalar Mlp = mag(MLeft)  > 1 ? 0.5 * (MLeft  + mag(MLeft))  :  0.25 * sqr(MLeft  + 1) ;
        // scalar Mrm = mag(MRight) > 1 ? 0.5 * (MRight - mag(MRight)) : -0.25 * sqr(MRight - 1) ;

        // // Eq (28)
        // scalar Plp = mag(MLeft)  > 1 ? 0.5 * (1 + sign(MLeft))  :  0.25 * sqr(MLeft  + 1) * (2 - MLeft)  + alpha * MLeft  * sqr(sqr(MLeft)  - 1);
        // scalar Prm = mag(MRight) > 1 ? 0.5 * (1 - sign(MRight)) : -0.25 * sqr(MRight - 1) * (2 + MRight) - alpha * MRight * sqr(sqr(MRight) - 1);

        // Eq (26)
        const scalar ps = Plp*pLeft + Prm*pRight;

        scalar fL, fR;
        if (ps > 0)
        {
            fL = pLeft/ps  - 1;
            fR = pRight/ps - 1;
        }
        else
        {
            fL = 0.0;
            fR = 0.0;
        }
        // const scalar fL = ps > 0 ? pLeft/ps  - 1 : 0.0;
        // const scalar fR = ps > 0 ? pRight/ps - 1 : 0.0;
        // Eq (10) in:  Weigang Yao, Min Xu, "Modified AUSMPW+ Scheme and its Application" 
        // scalar fL  = mag(MLeft)  < 1 ? pLeft/ps  - 1 : 0;
        // scalar fR  = mag(MRight) < 1 ? pRight/ps - 1 : 0;

        // Eq (25)
        const scalar w   = 1 - pow(min(pLeft/pRight, pRight/pLeft), 3);

        // m12 below Eq (13)
        const scalar m12 = Mlp + Mrm;

        // Below Eq (24)
        scalar barMlp, barMrm;

        if (m12 > 0)
        {
            barMlp = Mlp + Mrm*(1 - w*(1 + fL) + (fL - fR));
            barMrm = Mrm * w*(1+fR);
        }
        else
        {
            barMlp = Mlp * w*(1+fL);
            barMrm = Mrm + Mlp*(1 - w*(1 + fL) + (fL - fR));
        }
        // scalar barMlp = m12 > 0 ? Mlp + Mrm * ((1 - w)*(1 + fR) - fL) : Mlp * w * (1 + fL);
        // scalar barMrm = m12 > 0 ? Mrm * w * (1 + fR) : Mrm + Mlp * ((1 - w)*(1 + fL) - fR);


    	rhoFlux  = (barMlp * a12 * rhoLeft        + barMrm * a12 * rhoRight        ) * magSf;
    	rhoUFlux = (barMlp * a12 * rhoLeft*ULeft  + barMrm * a12 * rhoRight*URight ) * magSf + ps*Sf;
    	rhoEFlux = (barMlp * a12 * rhoLeft*HnLeft + barMrm * a12 * rhoRight*HnRight) * magSf;
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
