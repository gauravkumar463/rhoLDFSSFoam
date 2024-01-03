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
    const vector ULeft,
    const vector URight,
    const scalar rpsiLeft,
    const scalar rpsiRight,
    const scalar eLeft,
    const scalar eRight,
    const scalar kLeft,
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

    scalar qLeft = (ULeft & n)-(phi/magSf);
    scalar qRight = (URight & n)-(phi/magSf);

    if(scheme == "LDFSS"){
    	scalar rhoHalf = 0.5*(rhoLeft+rhoRight);
    	scalar aHalf = sqrt(0.5*(sqr(aLeft)+sqr(aRight)));

    	scalar MLeft = qLeft/aHalf;
    	scalar MRight = qRight/aHalf;
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

    	scalar betaLeft = -max(0.0,1-(int)mag(MLeft));
    	scalar betaRight = -max(0.0,1-(int)mag(MRight));

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

    	rhoFlux = magSf*(rhoLeft*UPlus + rhoRight*Uminus);

    	rhoUFlux = magSf*(rhoLeft*ULeft*UPlus + rhoRight*URight*Uminus) + pHalf*Sf;

    	rhoEFlux = magSf*(rhoLeft*HLeft*UPlus + rhoRight*HRight*Uminus);
    }
    else if((scheme == "Tadmor") || (scheme == "Kurganov")){
        scalar ap = max(max(qLeft+aLeft,qRight+aRight),0.0);
        scalar am = min(min(qLeft-aLeft,qRight-aRight),0.0);

        scalar aPos = ap/(ap-am);

        amaxSf = magSf*max(mag(am),mag(ap));

        scalar aSf = am*aPos;

        if(scheme == "Tadmor"){
            aSf = -0.5*amaxSf/magSf;
            aPos = 0.5;
        }

        scalar aNeg = 1.0-aPos;

        qLeft *= aPos;
        qRight *= aNeg;

        scalar aqLeft = qLeft - aSf;
        scalar aqRight = qRight + aSf;

        amaxSf = magSf*max(mag(aqLeft),mag(aqRight));

	//scalar HLeft = (rpsiLeft/(gamma-1)) + 0.5*magSqr(ULeft) + pLeft/rhoLeft;
        //scalar HRight = (rpsiRight/(gamma-1)) + 0.5*magSqr(URight) + pRight/rhoRight;

        scalar HLeft = eLeft + 0.5*magSqr(ULeft) + kLeft + pLeft/rhoLeft;
        scalar HRight = eRight + 0.5*magSqr(URight) + kRight + pRight/rhoRight;

        rhoFlux = magSf*(rhoLeft*aqLeft+ rhoRight*aqRight);

        rhoUFlux = magSf*(rhoLeft*ULeft*aqLeft + rhoRight*URight*aqRight) + (aPos*pLeft+aNeg*pRight)*Sf;

        rhoEFlux = magSf*(rhoLeft*HLeft*aqLeft + rhoRight*HRight*aqRight + aSf*(pLeft-pRight)) + phi*(aPos*pLeft + aNeg*pRight);

    //    Uf = aPos*ULeft + aNeg*URight;
    }
}

// ************************************************************************* //
