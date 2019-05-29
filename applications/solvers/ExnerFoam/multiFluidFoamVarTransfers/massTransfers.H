// The transfer terms are only formulated for nParts == 2
if (nParts != 2)
{
    FatalErrorIn("massTransfer.H")
        << " mass transfer terms only formulated for nParts = 2, not nParts = "
        << nParts << exit(FatalError);
}

// Set all mass transfers to zero
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=ip+1; jp < nParts; jp++)
    {
        massTransfer[ip][jp] == massTransferZero;
        massTransfer[jp][ip] == massTransferZero;
    }
}

// Mass tranfer terms as diffusion between partitions
if (Ksigma.value() > SMALL)
{
    for(label ip = 0; ip < nParts; ip++)
    {
        for(label jp=ip+1; jp < nParts; jp++)
        {
            massTransferTmp = 0.5*Ksigma*fvc::laplacian
            (
                sigmaRho[jp] -sigmaRho[ip]
            );
            massTransfer[ip][jp] += max(massTransferTmp, massTransferZero);
            massTransfer[jp][ip] -= min(massTransferTmp, massTransferZero);
        }
    }
}

// Mass transfer to avoid changes in sigmaRho due to partition divergence
if (divTransfer)
{
    partitionedSurfaceScalarField verticalVolFlux = volFlux;
    verticalVolFlux[0] = (fvc::interpolate( w[0]*Foam::vector(0,0,1), "Uf") & mesh.Sf());
    verticalVolFlux[1] = (fvc::interpolate( w[1]*Foam::vector(0,0,1), "Uf") & mesh.Sf());
    massTransfer[0][1] += divCoeff*sigmaRho[0]*max(-fvc::div(verticalVolFlux[0]), rateZero);
    massTransfer[1][0] += divCoeff*sigmaRho[1]*max(-fvc::div(verticalVolFlux[1]), rateZero);
    
// REverser it
//    massTransfer[0][1] += sigmaRho[0]*max(-fvc::div(volFlux[1]), rateZero);
//    massTransfer[1][0] += sigmaRho[1]*max(-fvc::div(volFlux[0]), rateZero);
}

// Mass transfer based on vertical velocity and horizontal divergence
if (wTransfer)
{
    massTransfer[0][1] -= max(sign(u[0].component(vector::Z)), dimlessZero)*
    min
    (
        fvc::div(sigmaFlux[0] - sigmaRhof[0]*((Uf[0] & ghat)*ghatSf)),
        massTransferZero
    );
    
    massTransfer[1][0] -= min(sign(u[1].component(vector::Z)), dimlessZero)*
    max
    (
        fvc::div(sigmaFlux[1] - sigmaRhof[1]*((Uf[1] & ghat)*ghatSf)),
        massTransferZero
    );
}

// Mass transfer based on buoyancy
if (thetaTransfer)
{
    massTransfer[0][1] -= sigmaRho[0]*thetaTransferDiffusivity*min
    (
        fvc::laplacian(theta[0])/theta[0],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );

    massTransfer[1][0] += sigmaRho[1]*thetaTransferDiffusivity*max
    (
        fvc::laplacian(theta[1])/theta[1],
        dimensionedScalar("", dimensionSet(0,-2,0,0,0), scalar(0))
    );
}

// Mass transfer based mixing due to drag
if (dragTransfer)
{
    massTransfer[0][1] += dragTransferCoeff * 
        sigmaRho[0]/max(cloudRadiusMin, cloudRadiusMax*sigma[0]*sigma[1]) * 
        sigmaRho.sum()/rho[1] * mag(u[0] - u[1]);

    massTransfer[1][0] += dragTransferCoeff * 
        sigmaRho[1]/max(cloudRadiusMin, cloudRadiusMax*sigma[0]*sigma[1]) * 
        sigmaRho.sum()/rho[0] * mag(u[0] - u[1]);
}

if (0)
{
    wVar[0] = max(wVar[0], 0*wVar[0]);
    wVar[1] = max(wVar[1], 0*wVar[1]);
    massTransfer[0][1] += 0.05 * sigmaRho[1] * 0.5 / wVarTimescale *
    (
        1 + Foam::erf
        (
            (w[0] - w[1]) / 
            max
            (
                dimensionedScalar("",w[0].dimensions(),scalar(1e-16)), 
                Foam::sqrt( 2 * (wVar[0] + wVar[1]) ) 
            )
        )
    );
    massTransfer[1][0] += 0.05 * sigmaRho[0] * 0.5 / wVarTimescale *
    (
        1 + Foam::erf
        (
            (w[0] - w[1]) / 
            max
            (
                dimensionedScalar("",w[0].dimensions(),scalar(1e-16)), 
                Foam::sqrt( 2 * (wVar[0] + wVar[1]) ) 
            )
        )
    );
}

sigmaRhoPreTransfer = sigmaRho;

#include "applyMassTransfer.H"

/*transfer[0] *= 0.;
transfer[1] *= 0.;

transferRate[0] *= 0.;
transferRate[1] *= 0.;*/

// Update the sums
sigmaRho.updateSum();
transfer.updateSum();

