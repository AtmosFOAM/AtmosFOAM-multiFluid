/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2014-2018 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "multiFluidThetaBuoyantKEpsilon.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "addToRunTimeSelectionTable.H"
#include "Partitioned.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace RASModels
{

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
multiFluidThetaBuoyantKEpsilon<BasicTurbulenceModel>::multiFluidThetaBuoyantKEpsilon
(
    const alphaField& alpha,
    const rhoField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    kEpsilon<BasicTurbulenceModel>
    (
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName,
        type
    ),

    Cg_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cg",
            this->coeffDict_,
            1.0
        )
    )
{
    if (type == typeName)
    {
        this->printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class BasicTurbulenceModel>
bool multiFluidThetaBuoyantKEpsilon<BasicTurbulenceModel>::read()
{
    if (kEpsilon<BasicTurbulenceModel>::read())
    {
        Cg_.readIfPresent(this->coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


template<class BasicTurbulenceModel>
tmp<volScalarField>
multiFluidThetaBuoyantKEpsilon<BasicTurbulenceModel>::Gcoef() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    const volScalarField& T = this->transport_.T();
    const volScalarField& rho = this->rho_;

    return 
        (Cg_*this->Cmu_)*this->alpha_*rho*this->k_*
        (
            g &
            (
                fvc::grad(rho)/rho
              - 1/(this->transport_.gamma()-1)*fvc::grad(T)/T
            )
        )
       /(this->epsilon_ + this->epsilonMin_);
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
multiFluidThetaBuoyantKEpsilon<BasicTurbulenceModel>::kSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    const 
     Partitioned<ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>>&
       turbulence =
       Partitioned<ThermalDiffusivity<PhaseCompressibleTurbulenceModel<fluidThermo>>>::New
       (this->mesh_);

    Info << "From turbulence " << this->k_.name()
         << " Turbulence partNames = " << turbulence.partNames() << endl;

    if (mag(g.value()) > small)
    {
        return -fvm::SuSp(Gcoef(), this->k_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::kSource();
    }
}


template<class BasicTurbulenceModel>
tmp<fvScalarMatrix>
multiFluidThetaBuoyantKEpsilon<BasicTurbulenceModel>::epsilonSource() const
{
    const uniformDimensionedVectorField& g =
        this->mesh_.objectRegistry::template
        lookupObject<uniformDimensionedVectorField>("g");

    if (mag(g.value()) > small)
    {
        vector gHat(g.value()/mag(g.value()));

/*        volScalarField w(gHat & this->U_);
        volScalarField uh
        (
            mag(this->U_ - gHat*w)
          + dimensionedScalar("small", dimVelocity, small)
        );

        return -fvm::SuSp(this->C1_*tanh(mag(w)/uh)*Gcoef(), this->epsilon_);
*/        return -fvm::SuSp(this->C1_*Gcoef(), this->epsilon_);
    }
    else
    {
        return kEpsilon<BasicTurbulenceModel>::epsilonSource();
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace Foam

// ************************************************************************* //
