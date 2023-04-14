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

Info<< "Reading rhoCp nu Pr\n" << endl;

dimensionedScalar rhoCp
(
	transportProperties.lookup("rhoCp")
);

dimensionedScalar nu
(
	transportProperties.lookup("nu")
);

dimensionedScalar Pr
(
	transportProperties.lookup("Pr")
);

volScalarField kappat
(
	IOobject
	(
		"kappat",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);

volScalarField T
(
	IOobject
	(
		"T",
		runTime.timeName(),
		mesh,
		IOobject::MUST_READ,
		IOobject::AUTO_WRITE
	),
	mesh
);
