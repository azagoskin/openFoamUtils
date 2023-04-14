/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | February, 2021 / Zagoskin Andrey 
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
	correctAspectRatio

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "timeSelector.H"
#include "cellSet.H"
#include "argList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// initialise a vector field to store deltas


int main(int argc, char *argv[])
{
	#include "setRootCase.H"
	#include "createTime.H"
	Foam::instantList timeDirs = Foam::timeSelector::selectIfPresent(runTime, args);
	#include "createNamedMesh.H"

	volScalarField aspectRatio
	(
		IOobject
		(
			"aspectRatio",
			runTime.timeName(),
			mesh,
			IOobject::NO_READ,
			IOobject::AUTO_WRITE
		),
		mesh,
		dimensionedScalar("aspectRatio", dimensionSet(0, 0, 0, 0, 0), 0)
	);

	Foam::timeSelector::addOptions();

        // read the dimensions of the cells

	forAll(timeDirs, timeI)
	{
	        runTime.setTime(timeDirs[timeI], timeI);
        	Foam::Info<< "Time = " << runTime.timeName() << Foam::endl;
		surfaceVectorField deltas = mesh.delta();
		faceList faces = mesh.faces();
		pointField points = mesh.points();
		surfaceVectorField faceCenters = mesh.Cf();
		
		forAll(mesh.cells(), cellI)
		{
			cell currentCell = mesh.cells()[cellI];
			edgeList currentEdges = currentCell.edges(faces);
			pointField currentPoints = currentCell.points(faces,points);
			vectorField edges(currentEdges.size(), vector::zero);
			forAll(currentEdges, edgeI)
			{
				edge currentEdge = currentEdges[edgeI];
				label lstart = currentEdge.start();
				label lend = currentEdge.end();
				point start = points[lstart];
				point end = points[lend];
				edges[edgeI] = start - end;
			}
			scalar maxDim = max(mag(edges));
			scalar minDim = min(mag(edges));
			aspectRatio[cellI] = maxDim / minDim;
		}
		aspectRatio.write();
	}

	Info<< "End" << endl;

	return 0;
}

// ************************************************************************* //
