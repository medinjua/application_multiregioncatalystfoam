tmp<fv::convectionScheme<scalar>> mvConvection(nullptr);

if (Y.size())
{
    mvConvection = tmp<fv::convectionScheme<scalar>>
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.divScheme("div(phi,Yi_h)")
        )
    );
}

reaction.correct();

// Construct correction to the diffusion velocity
surfaceScalarField JF
(
    IOobject
    (
        "JF",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::NO_WRITE
    ),
    mesh,
    dimensionedScalar
    (
        "zero",
        dimMass/dimArea/dimTime,
        0
    )
);

forAll(Y, i)
{
    JF += thermophysicalTransport.j(Y[i]);
}

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport.divj(Yi)
          - fvc::div(fvc::interpolate(Yi)*JF*mesh.magSf())  // Correction to the diffusion velocity
         ==
            // Following calculates reaction rate only on "catalyst" section
            epspore*reaction.R(Yi)        // Added multiplication by epspore
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
