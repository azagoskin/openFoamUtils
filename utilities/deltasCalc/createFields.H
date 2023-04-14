// initialise a vector field to store cell dimensions
volVectorField cellDim
(
    IOobject
    (
        "cellDim",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
	IOobject::AUTO_WRITE
    ),
    mesh,
    Foam::vector(0, 0, 0)
);

