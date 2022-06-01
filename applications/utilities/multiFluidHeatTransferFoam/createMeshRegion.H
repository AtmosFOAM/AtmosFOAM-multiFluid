    // Check for non-default mesh region
    const string meshRegion = args.optionFound("region") ?
                              args.optionRead<string>("region") :
                              fvMesh::defaultRegion;

    Info << "Create mesh for time = " << runTime.timeName() <<  " region "
         << meshRegion << endl;

    Foam::fvMesh mesh
    (
        Foam::IOobject
        (
            meshRegion,
            runTime.timeName(),
            runTime,
            IOobject::MUST_READ
        )
    );

    // List of cells to sum
    labelList sumCells;
    if (args.optionFound("cellSet"))
    {
        const word cellSetName(args.optionRead<string>("cellSet"));
        cellSet cells(mesh, cellSetName);
        sumCells = cells.toc();
    }
    else
    {
        // Select all cells
        sumCells.setSize(mesh.nCells());

        forAll(mesh.cells(), cellI)
        {
            sumCells[cellI] = cellI;
        }
    }

