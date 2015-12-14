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


    Info<< "Reading max. angular velocity omega\n" << endl;

    dimensionedVector Omega
    (
        transportProperties.lookup("Omega")
    );

    Info<< "Reading center of rotation\n" << endl;

    dimensionedVector Center
    (
        transportProperties.lookup("Center")
    );

    Info<< "Reading vertigo radius\n" << endl;

    dimensionedScalar radius
    (
        transportProperties.lookup("radius")
    );

