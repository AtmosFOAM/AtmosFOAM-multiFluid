/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2020 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    conditionalAveraging

Description
    Separates a fluid into rising and falling fields

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    Foam::argList::validArgs.append("conditionField");
    Foam::argList::validArgs.append("conditionValue");
    Foam::argList::validArgs.append("conditionNameBelow");
    Foam::argList::validArgs.append("conditionNameAbove");

    argList::addOption
    (
        "fields",
        "list",
        "specify a list of fields to be interpolated. Eg, '(U T p)' - "
        "regular expressions not currently supported"
    );
    
    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    const word conditionFieldName = args.args()[1].c_str();
    const scalar conditionValue = readScalar(IStringStream(args.args()[2])());
    wordList conditionName(2);
    conditionName[0] = args.args()[3].c_str();
    conditionName[1] = args.args()[4].c_str();
    
    Info << "Conditionally averaging based on " << conditionFieldName
         << ".\nBelow " << conditionValue << " is " << conditionName[0]
         << " above " << conditionValue << " is " << conditionName[1] << endl;

    HashSet<word> selectedFields;
    if (args.optionFound("fields"))
    {
        args.optionLookup("fields")() >> selectedFields;
    }
    if (selectedFields.size())
    {
        Info<< "Conditioanlly averaging fields " << selectedFields << endl;
    }
    else
    {
        Info<< "Calculating volume ratios" << endl;
    }

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField conditionField
        (
            IOobject(conditionFieldName, runTime.timeName(), mesh,
                     IOobject::MUST_READ),
            mesh
        );
        volVectorField gradCondition = fvc::grad(conditionField);
        
        partitionedVolScalarField sigma
        (
            "sigma",
            conditionName, 
            volScalarField
            (
                IOobject("sigma", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar("", dimless, scalar(0)),
                "zeroGradient"
            )
        );
        
        forAll(conditionField, cellI)
        {
            // Estimate the volume ratio rising and falling for each cell
            
            // Find the projection of the cell centre onto the w=0 surface
            point xw0 = mesh.C()[cellI]
                      - conditionField[cellI]*gradCondition[cellI]
                      /magSqr(gradCondition[cellI]);

            // Calulate the distance to the w=0 line for each vertex
            scalarList d(mesh.cellPoints()[cellI].size());
            labelList nNegPos(2, label(0));
            scalarList dTot(2, scalar(0));
            forAll(d, i)
            {
                const point& v = mesh.points()[mesh.cellPoints()[cellI][i]];
                d[i] = (v - xw0) & gradCondition[cellI]
                                   /mag(gradCondition[cellI]);
                if (d[i] < 0)
                {
                    nNegPos[0]++;
                    dTot[0] -= d[i];
                }
                else
                {
                    nNegPos[1]++;
                    dTot[1] += d[i];
                }
            }
            
            // The approximate volume fraction rising and falling
            sigma[0][cellI] = dTot[0]/(dTot[0] + dTot[1]);
            sigma[1][cellI] = dTot[1]/(dTot[0] + dTot[1]);
        }
        sigma.write();
    }
    
    Info<< nl << "ExecutionTime = " << runTime.elapsedCpuTime() << " s"
        << "  ClockTime = " << runTime.elapsedClockTime() << " s"
        << nl << endl;

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
