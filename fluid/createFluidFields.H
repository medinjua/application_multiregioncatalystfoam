// Initialise fluid field pointer lists
PtrList<fluidReactionThermo> thermoFluid(fluidRegions.size());
PtrList<volScalarField> rhoFluid(fluidRegions.size());
PtrList<volVectorField> UFluid(fluidRegions.size());
PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> hRefFluid(fluidRegions.size());
PtrList<uniformDimensionedScalarField> pRefFluid(fluidRegions.size());
PtrList<volScalarField> ghFluid(fluidRegions.size());
PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
PtrList<compressible::momentumTransportModel>
    turbulenceFluid(fluidRegions.size());
PtrList<fluidReactionThermophysicalTransportModel>
    thermophysicalTransportFluid(fluidRegions.size());
PtrList<combustionModel> reactionFluid(fluidRegions.size());
PtrList<volScalarField> p_rghFluid(fluidRegions.size());
PtrList<volScalarField> KFluid(fluidRegions.size());
PtrList<volScalarField> dpdtFluid(fluidRegions.size());
PtrList<multivariateSurfaceInterpolationScheme<scalar>::fieldTable>
    fieldsFluid(fluidRegions.size());

List<scalar> initialMassFluid(fluidRegions.size());

PtrList<IOMRFZoneList> MRFfluid(fluidRegions.size());
PtrList<Foam::fvModels> fvModelsFluid(fluidRegions.size());
PtrList<fvConstraints> fvConstraintsFluid(fluidRegions.size());

// *************************************************************************
// Additions for porous catalyst simulations
// chemistryModel required
PtrList<basicChemistryModel> chemistryFluid(fluidRegions.size());
// Heat release rate
PtrList<volScalarField> QdotFluid(fluidRegions.size());
// Porosity field
PtrList<volScalarField> epsporeFluid(fluidRegions.size());
// Auxiliary field for species i RR
PtrList<volScalarField> omegarr_i_Fluid(fluidRegions.size());
// Isothermal flow boolean list
List<bool> isothermalFluidRegions(fluidRegions.size());
// *************************************************************************

// Populate fluid field pointer lists
forAll(fluidRegions, i)
{
    Info<< "*** Reading fluid mesh thermophysical properties for region "
        << fluidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoFluid\n" << endl;
    thermoFluid.set(i, fluidReactionThermo::New(fluidRegions[i]).ptr());

    Info<< "    Adding to rhoFluid\n" << endl;
    rhoFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "rho",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            thermoFluid[i].rho()
        )
    );

    Info<< "    Adding to UFluid\n" << endl;
    UFluid.set
    (
        i,
        new volVectorField
        (
            IOobject
            (
                "U",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    Info<< "    Adding to phiFluid\n" << endl;
    phiFluid.set
    (
        i,
        new surfaceScalarField
        (
            IOobject
            (
                "phi",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            linearInterpolate(rhoFluid[i]*UFluid[i])
                & fluidRegions[i].Sf()
        )
    );

    Info<< "    Adding to gFluid\n" << endl;
    gFluid.set
    (
        i,
        new uniformDimensionedVectorField
        (
            IOobject
            (
                "g",
                runTime.constant(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::NO_WRITE
            )
        )
    );

    Info<< "    Adding to hRefFluid\n" << endl;
    hRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "hRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimLength, 0)
        )
    );

    Info<< "    Adding to pRefFluid\n" << endl;
    pRefFluid.set
    (
        i,
        new uniformDimensionedScalarField
        (
            IOobject
            (
                "pRef",
                runTime.constant(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::NO_WRITE
            ),
            dimensionedScalar(dimPressure, 0)
        )
    );

    dimensionedScalar ghRef(- mag(gFluid[i])*hRefFluid[i]);

    Info<< "    Adding to ghFluid\n" << endl;
    ghFluid.set
    (
        i,
        new volScalarField
        (
            "gh",
            (gFluid[i] & fluidRegions[i].C()) - ghRef
        )
    );

    Info<< "    Adding to ghfFluid\n" << endl;
    ghfFluid.set
    (
        i,
        new surfaceScalarField
        (
            "ghf",
            (gFluid[i] & fluidRegions[i].Cf()) - ghRef
        )
    );

    Info<< "    Adding to turbulenceFluid\n" << endl;
    turbulenceFluid.set
    (
        i,
        compressible::momentumTransportModel::New
        (
            rhoFluid[i],
            UFluid[i],
            phiFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to thermophysicalTransport\n" << endl;
    thermophysicalTransportFluid.set
    (
        i,
        fluidReactionThermophysicalTransportModel::New
        (
            turbulenceFluid[i],
            thermoFluid[i]
        ).ptr()
    );

    Info<< "    Adding to reactionFluid\n" << endl;
    reactionFluid.set
    (
        i,
        combustionModel::New(thermoFluid[i], turbulenceFluid[i])
    );

    p_rghFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "p_rgh",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    fluidRegions[i].setFluxRequired(p_rghFluid[i].name());

    initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

    Info<< "    Adding to KFluid\n" << endl;
    KFluid.set
    (
        i,
        new volScalarField
        (
            "K",
            0.5*magSqr(UFluid[i])
        )
    );

    Info<< "    Adding to dpdtFluid\n" << endl;
    dpdtFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "dpdt",
                runTime.timeName(),
                fluidRegions[i]
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                thermoFluid[i].p().dimensions()/dimTime,
                0
            )
        )
    );

    Info<< "    Adding to fieldsFluid\n" << endl;
    fieldsFluid.set
    (
        i,
        new multivariateSurfaceInterpolationScheme<scalar>::fieldTable
    );
    forAll(thermoFluid[i].composition().Y(), j)
    {
        fieldsFluid[i].add(thermoFluid[i].composition().Y()[j]);
    }
    fieldsFluid[i].add(thermoFluid[i].he());

    Info<< "    Adding MRF\n" << endl;
    MRFfluid.set
    (
        i,
        new IOMRFZoneList(fluidRegions[i])
    );

    Info<< "    Adding fvModelsFluid\n" << endl;
    fvModelsFluid.set
    (
        i,
        new Foam::fvModels(fluidRegions[i])
    );

    Info<< "    Adding fvConstraintsFluid\n" << endl;
    fvConstraintsFluid.set
    (
        i,
        new fvConstraints(fluidRegions[i])
    );

    turbulenceFluid[i].validate();


    // *********************************************************************
    // Additions for porous catalyst simulations
    // *********************************************************************

    // chemistryModel required
    Info<< "    Adding to chemistryFluid\n" << endl;
    chemistryFluid.set(i, basicChemistryModel::New(thermoFluid[i]));

    // Heat release rate
    Info<< "    Adding to QdotFluid\n" << endl;
    QdotFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "Qdot",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar
            (
                dimensionSet(1,-1,-3,0,0,0, 0),
                0
            )
        )
    );

    // porosity field
    Info<< "    Adding to epsporeFluid\n" << endl;
    epsporeFluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epspore",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            fluidRegions[i]
        )
    );

    // Field for species i reaction rate [kg/m3/s]
    Info<< "    Adding to omegarr_i_Fluid\n" << endl;
    omegarr_i_Fluid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "omegarr_i",
                runTime.timeName(),
                fluidRegions[i],
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            fluidRegions[i],
            dimensionedScalar("omegarr_i",dimensionSet(1,-3,-1,0,0,0, 0),0)
        )
    );

    // Check the existence of the cellZone "catalyst"
    const word cellSetCatalyst = "catalyst";
    label zoneID = fluidRegions[i].cellZones().findZoneID(cellSetCatalyst);
    if (zoneID == -1)
    {
        FatalErrorIn("createFluidFields")
            << "Cannot find cellZone " << cellSetCatalyst << endl
            << "Valid cellZones are " << fluidRegions[i].cellZones().names()
            << exit(FatalError);
    }

    // Check for option to run isothermal simulation (deactivates energy eqn)
    isothermalFluidRegions[i] = readBool(thermoFluid[i].properties().lookup("isothermalFlow"));

}


