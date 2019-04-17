/*--------------------------------*- C++ -*----------------------------------*\
| =========                 |                                                 |
| \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox           |
|  \\    /   O peration     | Version:  4.x                                   |
|   \\  /    A nd           | Web:      www.OpenFOAM.org                      |
|    \\/     M anipulation  |                                                 |
\*---------------------------------------------------------------------------*/

singleGraph
{
    start   (0.5 0.0005 0);
    end     (0.5   0.0005 1);
    fields  (U T);

    //start   (0 0.0005 0.5);
    //end     (0 0.0005 0.5);
    //fields  (T);

    #includeEtc "caseDicts/postProcessing/graphs/sampleDict.cfg"

    setConfig
    {
        axis    z;
    }

    // Must be last entry
    #includeEtc "caseDicts/postProcessing/graphs/graph.cfg"
}

// ************************************************************************* //
