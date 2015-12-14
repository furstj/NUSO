    Info<< "Reading field T\n" << endl;

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


    Info<< "Reading diffusivity DT\n" << endl;

    dimensionedScalar DT
    (
        transportProperties.lookup("DT")
    );

    dimensionedScalar K
    (
        transportProperties.lookup("K")
    );

    dimensionedScalar p
    (
        transportProperties.lookup("p")
    );

    surfaceScalarField c
    (
        IOobject
        (
            "c",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
	 ),
        mesh,
	DT
     );
