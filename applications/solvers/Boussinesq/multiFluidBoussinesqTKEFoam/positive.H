// Return 1 if positive, zero otherwise
namespace Foam
{
tmp<volScalarField> positive(const volScalarField& f)
{
    tmp<volScalarField> tpos
    (
        new volScalarField
        (
            IOobject("positive("+f.name()+')', f.instance(), f.mesh()),
            f.mesh(),
            dimensionedScalar("0", dimless, scalar(0)),
            calculatedFvPatchField<scalar>::typeName
        )
    );
    volScalarField& pos = tpos.ref();
    
    forAll(pos, cellI)
    {
        if (f[cellI] > 0) pos[cellI] = 1;
    }
    forAll(f.boundaryField(), patchi)
    {
        forAll(f.boundaryField()[patchi], faceI)
        {
            if (f.boundaryField()[patchi][faceI] > 0)
            {
                pos.boundaryFieldRef()[patchi][faceI] = 1;
            }
        }
    }

    return tpos;
}

}
