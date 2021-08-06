{
    Info << "rho minmax: [" << min(sigmaRho.sum()).value() << ", " 
                               << max(sigmaRho.sum()).value() << "]" << endl;
    Info << "rho Tot: " << (fvc::domainIntegrate(sigmaRho.sum())/Vtot).value() << endl;
    
    Info << "sigma0 minmax: [" << min(sigmaRho[0]/sigmaRho.sum()).value() << ", " 
                               << max(sigmaRho[0]/sigmaRho.sum()).value() << "]" << endl;
    Info << "sigma0 minmax: [" << min(sigma[0]).value() << ", " 
                               << max(sigma[0]).value() << "]" << endl;
    Info << "sigma1 minmax: [" << min(sigma[1]).value() << ", " 
                               << max(sigma[1]).value() << "]" << endl;
    
    
    Info << "theta0 minmax: [" << min(theta[0]).value() << ", " 
                               << max(theta[0]).value() << "]" << endl;
    Info << "theta1 minmax: [" << min(theta[1]).value() << ", " 
                               << max(theta[1]).value() << "]" << endl;
                               
    Info << "volFlux0 minmax: [" << min(volFlux[0]).value() << ", " 
                                 << max(volFlux[0]).value() << "]" << endl;
    Info << "volFlux1 minmax: [" << min(volFlux[1]).value() << ", " 
                                 << max(volFlux[1]).value() << "]" << endl;
    
    volScalarField PE0 = -sigmaRho[0]*(g & mesh.C());
    volScalarField PE1 = -sigmaRho[1]*(g & mesh.C());
    Info << "PE0 minmax: [" << min(PE0).value() << ", " 
                               << max(PE0).value() << "]" << endl;
    Info << "PE1 minmax: [" << min(PE1).value() << ", " 
                               << max(PE1).value() << "]" << endl;
    Info << "PE Tot: " << (fvc::domainIntegrate(PE0+PE1)/Vtot).value() << endl;
    
    volScalarField IE0 = sigmaRho[0]*theta[0]*Exner*Cv;
    volScalarField IE1 = sigmaRho[1]*theta[1]*Exner*Cv;
    Info << "IE0 minmax: [" << min(IE0).value() << ", " 
                               << max(IE0).value() << "]" << endl;
    Info << "IE1 minmax: [" << min(IE1).value() << ", " 
                               << max(IE1).value() << "]" << endl;
    Info << "IE Tot: " << (fvc::domainIntegrate(IE0+IE1)/Vtot).value() << endl;
    
    
    surfaceScalarField KE0 = 0.5*sigmaRhof[0]*sqr(volFlux[0]/mesh.magSf());
    surfaceScalarField KE1 = 0.5*sigmaRhof[1]*sqr(volFlux[1]/mesh.magSf());
    Info << "KE0 minmax: [" << min(KE0).value() << ", " 
                               << max(KE0).value() << "]" << endl;
    Info << "KE1 minmax: [" << min(KE1).value() << ", " 
                               << max(KE1).value() << "]" << endl;
    Info << "KE Tot: " << sum(KE0+KE1).value()/KE0.size() << endl;
    
    
    Info << endl;

}
