/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | December 2013, Zagoskin Andrey
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

\*---------------------------------------------------------------------------*/

#include "IOstream.H"
#include "fvCFD.H"
#include "simpleControl.H"

//***************************************************************************//
using namespace Foam;

int main(int argc, char *argv[])
{
#include "setRootCase.H"
#include "createTime.H"
#include "fvMesh.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "createFields.H"
#include "createPhi.H"

simpleControl simple(mesh);

const volScalarField radius (mesh.C().component(vector::X));

Info << "\nCalculating circulation distribution\n" << endl;

while (simple.loop())
{
	Info<< "Time = " << runTime.timeName() << nl << endl;

        while (simple.correctNonOrthogonal())
        {
		fvScalarMatrix circulation_
		(
			fvm::ddt(circulation)
			+ fvm::div(phi,circulation)
			==
			fvm::laplacian(nut,circulation)
			- (nut * pow(radius,-1) * mag(fvc::grad(circulation)))
		);
		circulation_.relax();
		circulation_.solve();
	}
	nut = pow(alpha * radius,2)*mag(fvc::grad(circulation * pow(radius,-1)));
	runTime.write();
}

Info<< "End\n" << endl;
return 0;
}
