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
    conditionalAveraging3

Description
    Calculates the volume fraction (sigma) in 3 partitions
    part 0 named name0 if conditionField1 < conditionValue1
    part 1 named name1 if conditionField1 > conditionValue1
                       && conditionField2 > conditionValue2
    part 2 named name2 if conditionField1 > conditionValue1
                       && conditionField2 < conditionValue2

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "PartitionedFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


int main(int argc, char *argv[])
{
    timeSelector::addOptions();
    #include "addRegionOption.H"
    Foam::argList::validArgs.append("name0");
    Foam::argList::validArgs.append("name1");
    Foam::argList::validArgs.append("name2");
    Foam::argList::validArgs.append("conditionField1");
    Foam::argList::validArgs.append("conditionValue1");
    Foam::argList::validArgs.append("conditionField2");
    Foam::argList::validArgs.append("conditionValue2");

    #include "setRootCase.H"
    #include "createTime.H"
    instantList timeDirs = timeSelector::select0(runTime, args);

    wordList partName(3);
    partName[0] = args.args()[1].c_str();
    partName[1] = args.args()[2].c_str();
    partName[2] = args.args()[3].c_str();
    const word conditionFieldName1 = args.args()[4].c_str();
    const scalar conditionValue1 = readScalar(IStringStream(args.args()[5])());
    const word conditionFieldName2 = args.args()[6].c_str();
    const scalar conditionValue2 = readScalar(IStringStream(args.args()[7])());
    
    Info << "Conditionally averaging based on " << conditionFieldName1 << " and "
         << conditionFieldName2 << endl;

    #include "createNamedMesh.H"

    forAll(timeDirs, timeI)
    {
        runTime.setTime(timeDirs[timeI], timeI);
        Info<< "Time = " << runTime.timeName() << endl;
        mesh.readUpdate();

        volScalarField conditionField1
        (
            IOobject(conditionFieldName1, runTime.timeName(), mesh,
                     IOobject::MUST_READ),
            mesh
        );
        volScalarField conditionField2
        (
            IOobject(conditionFieldName2, runTime.timeName(), mesh,
                     IOobject::MUST_READ),
            mesh
        );
        
        partitionedVolScalarField sigma
        (
            "sigma",
            partName, 
            volScalarField
            (
                IOobject("sigma", runTime.timeName(), mesh),
                mesh,
                dimensionedScalar("", dimless, scalar(1)),
                "zeroGradient"
            )
        );
        
        forAll(conditionField1, cellI)
        {
            // Find which partition this cell is in
            if (conditionField1[cellI] < conditionValue1)
            {
                sigma[0][cellI] = 1;
                sigma[1][cellI] = 0;
                sigma[2][cellI] = 0;
            }
            else if (conditionField2[cellI] > conditionValue2)
            {
                sigma[0][cellI] = 0;
                sigma[1][cellI] = 1;
                sigma[2][cellI] = 0;
            }
            else
            {
                sigma[0][cellI] = 0;
                sigma[1][cellI] = 0;
                sigma[2][cellI] = 1;
            }
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
