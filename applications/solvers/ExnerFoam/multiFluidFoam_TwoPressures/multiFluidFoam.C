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
    multiFluidFoam

Description
    Transient Solver for dry, buoyant, Navier-Stokes equations with constant
    dynamic viscosity partitioned into multiple fluids. 
    Solution using two Exner pressures and multiple theta and velocity, u
    Advective form momentum and theta equations. 

\*---------------------------------------------------------------------------*/

#include "HodgeOps.H"
#include "fvCFD.H"
#include "ExnerTheta.H"
#include "PartitionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "setRootCase.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "readEnvironmentalProperties.H"
    #include "readThermoProperties.H"
    #include "readTransferCoeffs.H"
    HodgeOps H(mesh);
    #define dt runTime.deltaT()
    #include "createFields.H"
    
    const dictionary& itsDict = mesh.solutionDict().subDict("iterations");
    const int nOuterCorr = itsDict.lookupOrDefault<int>("nOuterCorrectors", 2);
    const int nCorr = itsDict.lookupOrDefault<int>("nCorrectors", 1);
    const int nNonOrthCorr =
        itsDict.lookupOrDefault<int>("nNonOrthogonalCorrectors", 0);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    Info<< "\nStarting time loop\n" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "partitionedCourantNo.H"

        for (int ucorr=0; ucorr < nOuterCorr; ucorr++)
        {
            #include "sigmaEqn.H"
            #include "thetaEqn.H"
            #include "calculateDrag.H"
            #include "exnerEqn.H"
        }
        runTime.write();
        Info << "sigma[1] goes from " << min(sigma[1].internalField()).value()
             << " to " << max(sigma[1].internalField()).value()
             << " sumSigma goes from " << min(sigma.sum().internalField()).value()
             << " to "  << max(sigma.sum().internalField()).value() << endl;

        Info<< "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
            << "  ClockTime = " << runTime.elapsedClockTime() << " s"
            << nl << endl;
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
