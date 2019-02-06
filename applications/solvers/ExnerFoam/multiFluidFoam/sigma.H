// Updating sigma and rho in each partition and sigma on the face

rhoTheta = pRef/R*pow(Exner, (1-kappa)/kappa);
for(label ip = 0; ip < nParts; ip++)
{
    // Update the density so that pressure matches in both partitions
    rho[ip] = rhoTheta/theta[ip];

    // update sigma from rho and sigmaRho
    sigma[ip] = sigmaRho[ip]/rho[ip];
}
sigma.updateSum();
for(label ip = 0; ip < nParts; ip++)
{
    sigma[ip] /= sigma.sum();
    sigmaRho[ip] = sigma[ip]*rho[ip];
}
sigma.updateSum();
rho.updateSum();
sigmaRho.updateSum();
theta.updateSum();

// Update sigma on the face and limit to between zero and one
for(label ip = 0; ip < nParts; ip++)
{
    sigmaf[ip] = fvc::interpolate(sigma[ip], "sigma");
    sigmaf[ip] = min(max(sigmaf[ip], scalar(0)), scalar(1));
    sigmaRhof[ip] = fvc::interpolate(sigmaRho[ip]);
}
sigmaf.updateSum();
sigmaRhof.updateSum();

