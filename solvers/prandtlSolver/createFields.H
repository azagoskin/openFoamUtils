#include "createMesh.H"

Info<< "Reading transportProperties\n" << endl;

IOdictionary transportProperties
        (
        IOobject
                (
                "transportProperties",
                runTime.constant(),
                mesh,
                IOobject::MUST_READ_IF_MODIFIED,
                IOobject::NO_WRITE
                )
        );

Info<< "Reading circulation tangential velocity field\n" << endl;

volScalarField circulation
        (
        IOobject
                (
                "circulation",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
        mesh
        );

Info<< "Reading radial velocity field\n" << endl;

volVectorField U
        (
        IOobject
                (
                "U",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
        mesh
        );


Info<< "Reading alpha\n" << endl;

	dimensionedScalar alpha
	(
	transportProperties.lookup("alpha")
	);

volScalarField nut
        (
        IOobject
                (
                "nut",
                runTime.timeName(),
                mesh,
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
                ),
        mesh
        );

