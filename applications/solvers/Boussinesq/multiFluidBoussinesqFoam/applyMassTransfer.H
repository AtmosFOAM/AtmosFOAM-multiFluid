{
// Constraint to keep sigma => 0 and to ensure not more than all transferred
// in one time step (dt sigma S <= sigma)
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp = 0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            // Do not transfer out more than there exists
            massTransfer[ip][jp] = min
            (
                massTransfer[ip][jp],
                max(sigma[ip]-minSigma, dimlessZero)/dt
            );

            // Increase values that are <minSigma (shouldn't be needed)
            massTransfer[ip][jp] = max
            (
                massTransfer[ip][jp],
                max(minSigma-sigma[jp], dimlessZero)/dt
            );
        }
    }
}

// Constrain massTransferD and remove from massTransfer
if (Ksigma.value() > SMALL)
{
    massTransferD[0] = min(massTransferD[0], massTransfer[0][1]);
    massTransferD[1] = min(massTransferD[1], massTransfer[1][0]);
    massTransfer[0][1] -= massTransferD[0];
    massTransfer[1][0] -= massTransferD[1];
}

// Transfer terms for other variables
transfer[0] = dt*massTransfer[0][1]/max(sigma[1], minSigmaDiv);
transfer[1] = dt*massTransfer[1][0]/max(sigma[0], minSigmaDiv);
S[0] = dt*massTransfer[0][1]/max(sigma[0], minSigmaDiv);
S[1] = dt*massTransfer[1][0]/max(sigma[1], minSigmaDiv);
transfer.updateSum();
S.updateSum();

transferD[0] = massTransferD[0]/max(sigma[1], minSigmaDiv);
transferD[1] = massTransferD[1]/max(sigma[0], minSigmaDiv);

// Transfer terms on faces
massTransferf[0] = fvc::interpolate(massTransfer[0][1], "transfer");
massTransferf[1] = fvc::interpolate(massTransfer[1][0], "transfer");

Sf[0] = dt*massTransferf[0]/max(sigmaf[0], minSigmaDiv);
Sf[1] = dt*massTransferf[1]/max(sigmaf[1], minSigmaDiv);

transferf[0] = dt*massTransferf[0]/max(sigmaf[1], minSigmaDiv);
transferf[1] = dt*massTransferf[1]/max(sigmaf[0],minSigmaDiv);
transferf.updateSum();

// Apply mass transfer terms (operator split) to sigma equation
for(label ip = 0; ip < nParts; ip++)
{
    for(label jp=0; jp < nParts; jp++)
    {
        if (ip != jp)
        {
            sigma[ip] += dt*
            (
                massTransfer[jp][ip] - massTransfer[ip][jp]
              + massTransferD[jp] - massTransferD[ip]
            );
        }
    }
    sigma[ip].correctBoundaryConditions();
}
}
