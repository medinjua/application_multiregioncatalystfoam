    const fvMesh& mesh = fluidRegions[i];

    combustionModel& reaction = reactionFluid[i];

    fluidReactionThermo& thermo = thermoFluid[i];
    thermo.validate(args.executable(), "h", "e");

    basicSpecieMixture& composition = thermo.composition();
    PtrList<volScalarField>& Y = composition.Y();

    volScalarField& rho = rhoFluid[i];
    volVectorField& U = UFluid[i];
    surfaceScalarField& phi = phiFluid[i];

    compressible::momentumTransportModel& turbulence = turbulenceFluid[i];
    fluidReactionThermophysicalTransportModel& thermophysicalTransport =
        thermophysicalTransportFluid[i];

    volScalarField& K = KFluid[i];
    volScalarField& dpdt = dpdtFluid[i];

    volScalarField& p = thermo.p();
    const volScalarField& psi = thermo.psi();

    volScalarField& p_rgh = p_rghFluid[i];

    const dimensionedVector& g = gFluid[i];
    const dimensionedScalar& pRef = pRefFluid[i];
    const volScalarField& gh = ghFluid[i];
    const surfaceScalarField& ghf = ghfFluid[i];

    multivariateSurfaceInterpolationScheme<scalar>::fieldTable& fields =
        fieldsFluid[i];

    IOMRFZoneList& MRF = MRFfluid[i];
    Foam::fvModels& fvModels = fvModelsFluid[i];
    Foam::fvConstraints& fvConstraints = fvConstraintsFluid[i];

    #include "checkRadiationModel.H"

    const dimensionedScalar initialMass
    (
        "initialMass",
        dimMass,
        initialMassFluid[i]
    );

    pimpleNoLoopControl& pimple = pimples.pimple(i);

    pressureReference& pressureReference = pressureReferenceFluid[i];

    scalar cumulativeContErr = cumulativeContErrs[i];

    // This solver does not support moving mesh but it uses the pressure
    // equation of one which does, so we need a dummy face-momentum field
    autoPtr<surfaceVectorField> rhoUf(nullptr);


    // *********************************************************************
    // Additions for porous catalyst simulations
    // *********************************************************************

    // chemistry model required
    basicChemistryModel& chemistry = chemistryFluid[i];

    // Define cell set list with catalyst zone
    const word cellSetCatalyst = "catalyst";
    label zoneID = fluidRegions[i].cellZones().findZoneID(cellSetCatalyst);
    const labelList& cellsCatalyst = fluidRegions[i].cellZones()[zoneID];

    // species reaction rate auxiliary field
    volScalarField& omegarr_i = omegarr_i_Fluid[i]; 

    // heat release rate field
    volScalarField& Qdot = QdotFluid[i];