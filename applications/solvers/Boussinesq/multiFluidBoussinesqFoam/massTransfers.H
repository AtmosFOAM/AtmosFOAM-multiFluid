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
        volScalarField transferTmp = max(sigma[ip]-minSigma[ip], scalar(0))
                       *transferScale*sqrt(max(-dbdz, sqr(rateZero)));
        //volScalarField DbDt = Q + fvc::laplacian(alphaf,b[ip]);
    
        // Transfer out of 0 if dbdz<0. Transfer to 1 if dbdt>0, otherwise to 2
        if (ip == 0)
        {
            volScalarField d2bdz2 = volVectorField(fvc::grad(dbdz)).component(2);
//            volScalarField wdwdz = //u[0].component(2)*fvc::div(volFlux[0]*magk);
//                                  fvc::div(sqr(volFlux[0])*kdir);
//            massTransfer(0,1) += transferTmp*0.5*(1+sign(b.ddt()[0]-dbdtSMALL));
//            massTransfer(0,2) += transferTmp*0.5*(1-sign(b.ddt()[0]+dbdtSMALL));
//            massTransfer(0,1) += transferTmp*0.5*(1+sign(b[0]));
//            massTransfer(0,2) += transferTmp*0.5*(1-sign(b[0]));
//            massTransfer(0,1) += transferTmp*0.5*(1+sign(DbDt));
//            massTransfer(0,2) += transferTmp*0.5*(1-sign(DbDt));
            massTransfer(0,1) += transferTmp*0.5*(1+sign(d2bdz2));
            massTransfer(0,2) += transferTmp*(1-0.5*(1+sign(d2bdz2)));
//            massTransfer(0,1) += transferTmp*0.5*(1+sign(wdwdz));
//            massTransfer(0,2) += transferTmp*(1-0.5*(1+sign(wdwdz)));
        }
        // Transfers between 1 and 2 if unstable and wrong order
        else if (ip == 1)
        {
            massTransfer(1,2) += transferTmp*0.5*(1+sign(b[0] - b[1]));
        }
        else if (ip == 2)
        {
            massTransfer(2,1) += transferTmp*0.5*(1+sign(b[2] - b[0]));
        }
    }

    // Detrainment and Entrainment
    for(label ip = 1; ip < nParts; ip++)
    {
        volScalarField dwdz = fvc::div(volFlux[ip]*magk);
        // Detrainment
        massTransfer(ip,0) += detrainRate
            *max(sigma[ip]-minSigma[ip], scalar(0))
            *max(-dwdz, rateZero);
    
        // Entrainment
        if (ip == 1)
        {
            massTransfer(0,1) += entRate*0.5*(1+sign(u[1].component(2)))*min
            (
                max(sigma[0]-minSigma[0], scalar(0)),
                max(sigma[1] - minSigma[1], scalar(0))
            )*max(dwdz, rateZero);
        }
        else if (ip == 2)
        {
            massTransfer(0,2) += entRate*0.5*(1-sign(u[2].component(2)))*min
            (
                max(sigma[0]-minSigma[0], scalar(0)),
                max(sigma[2] - minSigma[2], scalar(0))
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

