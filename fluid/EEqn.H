{
    if(!isothermalFluidRegions[i]){
    
        //**********************************************************************
        // Modification catalyst simulations
        // Reinitialize Qdot to zero
        Qdot *= 0.0;

        // First, calculate the heat release rate by the chemistry pointer
        // Note the HRR units (J/m3/s)
        // Note: a call to chemistry.calculate() must have been done before
        volScalarField::Internal HRR(chemistry.Qdot().ref());
        // Now set the HRR in the catalyst cells
        forAll(cellsCatalyst, cc)
        {
            label celli = cellsCatalyst[cc];
            Qdot[celli] = HRR[celli];
        }
        //******************************************************************

        volScalarField& he = thermo.he();

        fvScalarMatrix EEqn
        (
            fvm::ddt(rho, he)
          + (
                mvConvection.valid()
                ? mvConvection->fvmDiv(phi, he)
                : fvm::div(phi, he)
            )
          + fvc::ddt(rho, K) + fvc::div(phi, K)
          + (
                he.name() == "e"
                ? mvConvection.valid()
                ? mvConvection->fvcDiv(fvc::absolute(phi, rho, U), p/rho)
                : fvc::div(fvc::absolute(phi, rho, U), p/rho)
                : -dpdt
            )
          + thermophysicalTransport.divq(he)
          ==
            rho*(U&g)
            //+ reaction.Qdot()        // commented out for replacement **********
            + Qdot                     // replaced HRR             ***************
            + fvModels.source(rho, he)
        );

        EEqn.relax();

        fvConstraints.constrain(EEqn);

        EEqn.solve();

        fvConstraints.constrain(he);

        thermo.correct();

    }

    Info<< "Min/max T:" << min(thermo.T()).value() << ' '
        << max(thermo.T()).value() << endl;
}
