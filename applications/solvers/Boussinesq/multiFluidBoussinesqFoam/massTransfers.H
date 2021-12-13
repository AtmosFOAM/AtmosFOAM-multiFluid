// Set all mass transfers to zero
if (transferType != noTransfer) for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer(ip,jp) == rateZero;
        massTransfer(jp,ip) == rateZero;
    }
}

if (transferType == dwdzTransfer && nParts == 2)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        massTransfer(ip,1-ip) = transferScale
                                *max(sigma[ip] - minSigma[ip], scalar(0))
                                *max(-dwdz, rateZero);
    }
}
else if (transferType == dwdzTransfer && nParts == 3)
{
    // Transfer from fluid 0 to 1 or 2 if dbdz<0
    volScalarField dbdz = fvc::div(bf[0]*magk);
    volScalarField DbDt = Q + fvc::laplacian(alphaf,b[0]);
    volScalarField transferTmp = max(sigma[0]-minSigma[0], scalar(0))
                       *transferScale*sqrt(max(-dbdz, sqr(rateZero)));
    
    // Transfer out of 0 if dbdz<0. Transfer to 1 if DbDt>0, otherwise to 2
    forAll(dbdz, cellI)
    {
        if (dbdz[cellI] < 0)
        {
            if (DbDt[cellI] > SMALL) // (b[0][cellI] > SMALL) //
            {
                massTransfer(0,1)[cellI] += transferTmp[cellI];
            }
            else if (DbDt[cellI] < -SMALL) // (b[0][cellI] < -SMALL) //
            {
                massTransfer(0,2)[cellI] += transferTmp[cellI];
            }
        }
    }
    
    // Transfer out of 1 or 2 into 0 if w is decelerating (detrainment)
    // and in if w accelerating (entrainment)
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        // Detrainment
        massTransfer(ip,0) += max(sigma[ip]-minSigma[ip], scalar(0))
            *max(-dwdz, rateZero);
        // Entrainment
        massTransfer(0,ip) += entRate*min
        (
            max(sigma[0]-minSigma[0], scalar(0)),
            sigma[ip]
        )*max(dwdz, rateZero);
    }
}

else if (transferType == dbdzTransfer && nParts == 3)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        volScalarField dbdz = fvc::div(bf[ip]*magk);
        volScalarField transferTmp = max
        (
            (sigma[ip]-minSigma[ip])/(1-minSigma[ip]), scalar(0)
        )*transferScale*sqrt(max(-dbdz, sqr(rateZero)));
        //volScalarField DbDt = Q + fvc::laplacian(alphaf,b[ip]);
        volScalarField w0 = u[0].component(2);
    
        // Transfer out of 0 if dbdz<0. Transfer to 1 if dbdt>0, otherwise to 2
        if (ip == 0)
        {
            // Calculate which fluid to transfer to. j=1 for rising, j=2 falling
            volScalarField wdwdz = u[0].component(2)*fvc::div(volFlux[0]*magk);
            
            volScalarField w0big = positive
            (
                mag(w0*dbdz) - mag(fvc::laplacian(alphaf,b[0]))
            );
            //w0big.rename("w0big"); w0big.write();
            //volScalarField dbdt("dbdt", b.ddt()[0]);
            
            volScalarField toOne = w0big*positive(wdwdz)
                                 + (1-w0big)*positive(b.ddt()[0]);
            //toOne.rename("toOne"); toOne.write();
        
            massTransfer(0,1) += toOne*transferTmp;
            massTransfer(0,2) += (1-toOne)*transferTmp;
        }
        // Transfers between 1 and 2 if unstable and wrong order
        else if (ip == 1)
        {
            volScalarField w1 = u[1].component(2);
            massTransfer(1,2) += transferTmp*0.5*(1+sign(b[0] - b[1])
                                                + 1+sign(-w1));
        }
        else if (ip == 2)
        {
            volScalarField w2 = u[2].component(2);
            massTransfer(2,1) += transferTmp*0.5*(1+sign(b[2] - b[0])
                                                + 1+sign(w2));
        }
    }

    // Detrainment and Entrainment
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        // Detrainment
        massTransfer(ip,0) += detrainRate
            *max
            (
                (sigma[ip]-minSigma[ip])/(1-minSigma[ip]),
                scalar(0)
            )
            *max(-dwdz, rateZero);
    
        // Entrainment
        if (ip == 1)
        {
            massTransfer(0,1) += entRate*0.5*(1+sign(u[1].component(2)))* //min
            (
                max((sigma[0]-minSigma[0])/(1-minSigma[0]), scalar(0))/*,
                max((sigma[1] - minSigma[1])/(1-minSigma[1]), scalar(0))*/
            )*max(dwdz, rateZero);
        }
        else if (ip == 2)
        {
            massTransfer(0,2) += entRate*0.5*(1-sign(u[2].component(2)))*//min
            (
                max((sigma[0]-minSigma[0])/(1-minSigma[0]), scalar(0))/*,
                max((sigma[2] - minSigma[2])/(1-minSigma[0]), scalar(0))*/
            )*max(dwdz, rateZero);
        }
    }
}

// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Smooth
            if (massTransferSmooth > 0)
            {
                fvc::smooth(massTransfer(ip,jp), massTransferSmooth);
            }
            // Do not transfer out more than exists
            massTransfer(ip,jp) = min
            (
                massTransfer(ip,jp),
                max(sigma[ip]-minSigma[ip], dimlessZero)/(dt*(nParts-1))
            );
        }
    }
}

