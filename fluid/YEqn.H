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

// Following is commented out for explicit calculation of RR in catalysts
//reaction.correct();
//**************************************************************************
// Modification catalyst simulations
// Instead, calculate regular reaction rate with the chemistry pointer
chemistry.calculate();
// End modification
//**************************************************************************

forAll(Y, i)
{
    if (composition.solve(i))
    {
        volScalarField& Yi = Y[i];

        //******************************************************************
        // regular reaction matrix
        // reaction rate obtained by chemistry.RR(i), 
        // previously set by chemistry.calculate()
        // NOTE: Reaction rates for species are in units kg/m3/s
        volScalarField::Internal RRi(chemistry.RR(i));

        // First, reinitialize omegarr_i to zero
        omegarr_i *= 0.0;
        forAll(cellsCatalyst, cc)
        {
            label celli       = cellsCatalyst[cc];
            omegarr_i[celli]  = RRi[celli];
        }
        //******************************************************************

        fvScalarMatrix YiEqn
        (
            fvm::ddt(rho, Yi) + mvConvection->fvmDiv(phi, Yi)
          + thermophysicalTransport.divj(Yi)
         ==
            //reaction.R(Yi)           // commented out for replacement ****
            omegarr_i                  // replaced RR matrix ***************
          + fvModels.source(rho, Yi)
        );

        YiEqn.relax();

        fvConstraints.constrain(YiEqn);

        YiEqn.solve("Yi");

        fvConstraints.constrain(Yi);
    }
}

composition.normalise();
