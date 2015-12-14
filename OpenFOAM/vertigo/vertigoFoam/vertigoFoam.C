/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011 OpenFOAM Foundation
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
    vertigoFoam

Description
    Solves a simple convection equation with const. convective velocity

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "simpleControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"

    #include "createTime.H"
    #include "createMesh.H"
    #include "createFields.H"

    simpleControl simple(mesh);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nCalculating vertigo\n" << endl;

    volScalarField r( mag( (mesh.C()-Center) ^ Omega) / mag(Omega)); 

    volVectorField U( pos(radius-r) * sqr(1-r/radius) * ( (mesh.C()-Center) ^ Omega ) );

    surfaceScalarField phi( 
			   IOobject
			   (
			    "phi",
			    runTime.timeName(),
			    mesh,
			    IOobject::READ_IF_PRESENT,
			    IOobject::AUTO_WRITE
			    ),
			   
			   linearInterpolate(U) & mesh.Sf() );
    
    while (simple.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;


	solve
	  (
	   fvm::ddt(T) + fvm::div(phi, T)
	   );

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
