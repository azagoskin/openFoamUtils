/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | January, 2014 / Zagoskin Andrey
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
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
	heatCalc

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "turbulenceModel.H"
#include "solidThermo.H"
#include "wallFvPatch.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
	timeSelector::addOptions();
	#include "setRootCase.H"
	#include "createTime.H"
	instantList timeDirs = timeSelector::select0(runTime, args);
	#include "createNamedMesh.H"

	forAll(timeDirs, timeI)
	{
		runTime.setTime(timeDirs[timeI], timeI);
		Info<< "Time = " << runTime.timeName() << endl;
		mesh.readUpdate();

		#include "createFields.H"

		surfaceScalarField heatFlux
		(
			fvc::interpolate(kappat+nu/Pr)*rhoCp*fvc::snGrad(T)
		);

		const surfaceScalarField::GeometricBoundaryField& patchHeatFlux =
			heatFlux.boundaryField();

		Info<< "\nWall heat fluxes [W]" << endl;
		forAll(patchHeatFlux, patchi)
		{
			Info<< mesh.boundary()[patchi].name()
			<< " "
			<< gSum
			(
				mesh.magSf().boundaryField()[patchi]
				*patchHeatFlux[patchi]
			)
			<< endl;
		}
		Info<< endl;


	        volScalarField wallHeatFlux
        	(
	            IOobject
	            (
	                "wallHeatFlux",
	                runTime.timeName(),
	                mesh
	            ),
	            mesh,
	            dimensionedScalar("wallHeatFlux", heatFlux.dimensions(), 0.0)
		    );

	}


	Info<< "End" << endl;

	return 0;
}

// ************************************************************************* //
