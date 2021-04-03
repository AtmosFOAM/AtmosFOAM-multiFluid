// Transfers to keep sigma[0] = 0 where w[0]>0 and vice verca
// This should be combined with transferring the mean of other variables

if (transferType == wTransfer && nParts == 2)
{
    massTransfer(0,1) = rateZero;
    massTransfer(1,0) = rateZero;

    volScalarField w = u[0].component(2);
    forAll(w, cellI)
    {
        if (w[cellI] > SMALL)
        {
            massTransfer(0,1)[cellI] += max(sigma[0][cellI]-minSigma[0], scalar(0))
                                        *transferScale/dt.value();
        }
    }
    
    w = u[1].component(2);
    forAll(w, cellI)
    {
        if (w[cellI] < -SMALL)
        {
            massTransfer(1,0)[cellI] += max(sigma[1][cellI]-minSigma[1], scalar(0))
                                        *transferScale/dt.value();
        }
    }

    // Apply mass transfers
    sigma.transferMass(massTransfer, dt);
    interpolate(sigmaf, sigma);

    // Transfer mean bouyancy
    b.transferField(massTransfer, dt);

    // Transfer mean momentum
    volFlux.transferField(interpolate(massTransferf, massTransfer), dt);
}
