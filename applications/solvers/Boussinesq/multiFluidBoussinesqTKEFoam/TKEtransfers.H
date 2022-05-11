// TKE and kinetic energy transfers

// Set TKE transferred to the mean and transfer resolved velocity differences
for(label ip = 0; ip < nParts; ip++) 
    for(label jp = 0; jp < nParts; jp++) if (ip != jp)
{
    tkeT(ip,jp) == TKE[ip];
    tkeAdd(ip,jp) = 0.5*magSqr(fvc::reconstruct(volFluxT(ip,jp) - volFlux[jp]));
    tkeRemove(ip,jp) = 0.5*magSqr(fvc::reconstruct(volFluxT(ip,jp) - volFlux[ip]));
}

// Transfer the TKE
TKE.asymetricTransferField(massTransfer, tkeT, tkeAdd, tkeRemove, dt);

/*for(label ip = 0; ip < nParts; ip++) 
{
    Info << "In TKEtransfers.H TKE[" << ip << "] goes from " << min(TKE[ip]).value()
         << " to " << max(TKE[ip]).value() << endl;
}*/
