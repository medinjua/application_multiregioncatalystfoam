// Initialise solid field pointer lists
PtrList<solidThermo> thermoSolid(solidRegions.size());
PtrList<Foam::fvModels> fvModelsSolid(solidRegions.size());
PtrList<fvConstraints> fvConstraintsSolid(solidRegions.size());

// *************************************************************************
// Additions for porous catalyst simulations
// Porosity field
PtrList<volScalarField> epsporeSolid(solidRegions.size());
// Isothermal flow boolean list
List<bool> isothermalSolidRegions(solidRegions.size());
// *************************************************************************

// Populate solid field pointer lists
forAll(solidRegions, i)
{
    Info<< "*** Reading solid mesh thermophysical properties for region "
        << solidRegions[i].name() << nl << endl;

    Info<< "    Adding to thermoSolid\n" << endl;
    thermoSolid.set(i, solidThermo::New(solidRegions[i]));

    Info<< "    Adding to fvModelsSolid\n" << endl;
    fvModelsSolid.set
    (
        i,
        new Foam::fvModels(solidRegions[i])
    );

    Info<< "    Adding fvConstraintsSolid\n" << endl;
    fvConstraintsSolid.set
    (
        i,
        new fvConstraints(solidRegions[i])
    );

    // *********************************************************************
    // Additions for porous catalyst simulations
    // *********************************************************************

    // porosity field
    Info<< "    Adding to epsporeSolid\n" << endl;
    epsporeSolid.set
    (
        i,
        new volScalarField
        (
            IOobject
            (
                "epspore",
                runTime.timeName(),
                solidRegions[i],
                IOobject::MUST_READ,
                IOobject::AUTO_WRITE
            ),
            solidRegions[i]
        )
    );

    // Check for option to run isothermal simulation (deactivates energy eqn)
    isothermalSolidRegions[i] = readBool(thermoSolid[i].properties().lookup("isothermalFlow"));
}
