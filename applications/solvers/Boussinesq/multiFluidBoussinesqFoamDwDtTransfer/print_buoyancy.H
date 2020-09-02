{
    forAll (b[0], celli)
    {
        if (0)
        // if ((celli > 20000) && (celli < 20033))
        {
            if (mag(b[0][celli]-b[1][celli]) > 1e-15)
            // if (1)
            {
                Info << endl;
                Info << "t=" << runTime.timeName() << ", ucorr=" << ucorr << ", counter=" << counter << ": sigma0[" << celli << "] = " << sigma[0][celli] << ", sigma1[" << celli << "] = " << sigma[1][celli] << endl;
                Info << "t=" << runTime.timeName() << ", ucorr=" << ucorr << ", counter=" << counter << ": b0[" << celli << "] = " << b[0][celli] << ", b1[" << celli << "] = " << b[1][celli] << endl;
                Info << endl;
            }
        }
    }
}