// Function which approximates inverse of error function
// to accuracy of ~6e-3.
// Input field must be dimensionless

volScalarField erfinv(volScalarField inputField)
{
    volScalarField inputFieldSign = sign(inputField);
    
    volScalarField tempField0 = -sqr(inputField) + scalar(1);
    volScalarField tempField1 = Foam::log( min(max(tempField0, 1e-16),1-1e-16) );
    volScalarField tempField2 = 0.5*tempField1 + 2/(0.147*M_PI);
    volScalarField tempField3 = tempField1/0.147;
    
    volScalarField outputField = inputFieldSign * sqrt
    ( 
        sqrt
        (
            sqr(tempField2) - tempField3
        )
      - tempField2 
    );
    
    return outputField;
}

/*template <typename T>
T erfinv(T inputField)
{
    T inputFieldSign = sign(inputField);
    
    //T tempField0 = -inputField*inputField;
    //tempField0 = tempField0 +  + scalar(1);
    
    T tempField1 = Foam::log( -sqr(inputField) + 1 );
    T tempField2 = 0.5*tempField1 + 2/(0.147*M_PI);
    T tempField3 = tempField1/0.147;
    
    T outputField = inputFieldSign * sqrt
    ( 
        sqrt
        (
            sqr(tempField2) - tempField3
        )
      - tempField2 
    );
    
    return outputField;
}*/
