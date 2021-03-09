/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016 OpenFOAM Foundation
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

Class
    Foam::TransferField

Description
    A matrix of GeometricFields to transfer between PartitionedFields
    without diagonal elements.

SourceFiles
    TransferField.C

\*---------------------------------------------------------------------------*/

#ifndef TransferField_H
#define TransferField_H

#include "GeometricField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class TransferField Declaration
\*---------------------------------------------------------------------------*/

template<class Type, template<class> class PatchField, class GeoMesh>
class TransferField
:
    public PtrList<GeometricField<Type, PatchField, GeoMesh> >
{
    // Private data

        //- The base name for all of the partitions
        const word baseName_;

        //- The names of the parts
        const wordList partNames_;

    // Dissallow copy construct
        TransferField(const TransferField&);


public:

    // Public typedefs

        //- Type of mesh on which this TransferField is instantiated
        //- Optionally write sum
        typedef typename GeoMesh::Mesh Mesh;
        
    // Constructors

        //- Read in based on names
        TransferField
        (
            const word& baseName__,
            const wordList& partNames__,
            const Mesh& mesh,
            const word& timeName,
            IOobject::writeOption writeOpt = IOobject::AUTO_WRITE
        );

        //- Initialise all parts to "field"
        TransferField
        (
            const word& baseName__,
            const wordList& partNames__,
            const GeometricField<Type, PatchField, GeoMesh>& field,
            IOobject::writeOption writeOpt = IOobject::AUTO_WRITE
        );


    //- Destructor
    ~TransferField();

    // Access functions
    
        //- The base name for all of the partitions
        const word& baseName() const { return baseName_; }

        //- The names of the parts
        const wordList& partNames() const { return partNames_; }

        //- The number of parts
        label nParts() const { return partNames().size(); }

    // Member Functions

        // Output
        void write();


    // Functions from PtrList

        const GeometricField<Type, PatchField, GeoMesh>& operator()
        (
            const label ip, label jp
        ) const;

        GeometricField<Type, PatchField, GeoMesh>& operator()
        (
            const label ip, label jp
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "TransferField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
