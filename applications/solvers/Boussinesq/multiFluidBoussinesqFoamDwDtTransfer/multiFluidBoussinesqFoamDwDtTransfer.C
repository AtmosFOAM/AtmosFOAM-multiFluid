/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 1991-2009 OpenCFD Ltd.
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM; if not, write to the Free Software Foundation,
    Inc., 51 Franklin St, Fifth Floor, Boston, MA 02110-1301 USA

Application
    multiFluidBoussinesqFoam

Description
    Transient Solver for dry, multi-fluid Boussinesq equations

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"
#include "fvcCurlf.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "zeros.H"
    #include "readTransferCoeffs.H"
    #include "readEnvironmentalProps.H"
    #define dt runTime.deltaT()
    #include "createFields.H"
    //#include "print_buoyancy.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);
    const scalar offCentre = readScalar(mesh.schemesDict().lookup("offCentre"));

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    /*u[0] *= 0;
    u[1] *= 0;
    volFlux[0] *= 0;
    volFlux[1] *= 0;
    forAll (b[0], celli)
    {
        if (celli >= 20000)
        {
            b[0][celli] *= 0;
            b[1][celli] *= 0;
        }
    }*/
    // forAll (volFlux[0], celli)
    // {
    //     if (1)
    //     {
    //         Info << endl;
    //         Info << celli << endl;
    //         Info << "P: " << P[celli] << endl;
    //         Info << "vF0: " << volFlux[0][celli]/2e7 << endl;
    //         Info << "vF1: " << volFlux[1][celli]/2e7 << endl;
    //         Info << "vFddt0: " << volFlux.ddt()[0][celli]/2e7 << endl;
    //         Info << "vFddt1: " << volFlux.ddt()[1][celli]/2e7 << endl;
            
    //     }
    // }
    
    
    Info<< "\nStarting time loop\n" << endl;

    

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;
        int counter = 0;
        #include "partitionedCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "print_buoyancy.H"
            #include "sigmaEqn.H"
            if (!noTransfers)
            {
                #include "massTransfers.H"
                #include "applyMassTransfer.H"
            }
            #include "print_buoyancy.H"
            #include "calculateDrag.H"
            #include "bEqn.H"
            #include "print_buoyancy.H"
            // Pressure and velocity updates
            for (int corr=0; corr<nCorr; corr++)
            {
                #include "PEqn.H"
                #include "momentumTransfers.H"
                //#include "PiEqn.H"
                // Update velocities based on the flux
                for(label ip = 0; ip < nParts; ip++)
                {
                    u[ip] = fvc::reconstruct(volFlux[ip]);
                }
                volScalarField u0 = u[0].component(2);
                volScalarField u1 = u[1].component(2);
                forAll (u0, celli)
                {
                    if (1)
                    {
                        Info << endl;
                        Info << celli << endl;
                        Info << "vF0: " << u0[celli] << endl;
                        Info << "vF1: " << u1[celli] << endl;
                    }
                }
            }
        }

        // Apply mass transfer terms (operator split) to sigmaf
        for(label ip = 0; ip < nParts; ip++)
        {
            sigmaf[ip] += dt*(massTransferf[1-ip] - massTransferf[ip]);
        }
        sigmaf.updateSum();
        volFlux.updateSum();

        // Update diagnositcs
        for(label ip = 0; ip < nParts; ip++)
        {
            Uf[ip] = linearInterpolate(u[ip]);
            Uf[ip] += (volFlux[ip] - (Uf[ip] & mesh.Sf()))
                      *mesh.Sf()/sqr(mesh.magSf());
            divu[ip] = fvc::div(sigmaf[ip]*volFlux[ip]);

            // Pressure gradient in each fluid including drag
            dPdz[ip] = mesh.Sf().component(2)/mesh.magSf()*
            (
                fvc::snGrad(P+Pi[ip])
              + (1-sigmaf[ip])*dragCommon*(volFlux[ip] - volFlux[1-ip])
                /mesh.magSf()
            );
        }
        Uf.updateSum();
        divu.updateSum();
        u.updateSum();

        runTime.write();

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
