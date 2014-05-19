    // Initialise fluid field pointer lists
    PtrList<rhoThermo> thermoFluid(fluidRegions.size());
    PtrList<volScalarField> rhoFluid(fluidRegions.size());
    PtrList<volVectorField> UFluid(fluidRegions.size());
    PtrList<surfaceScalarField> phiFluid(fluidRegions.size());
    PtrList<uniformDimensionedVectorField> gFluid(fluidRegions.size());
    PtrList<compressible::turbulenceModel> turbulence(fluidRegions.size());
    PtrList<volScalarField> p_rghFluid(fluidRegions.size());
    PtrList<volScalarField> ghFluid(fluidRegions.size());
    PtrList<surfaceScalarField> ghfFluid(fluidRegions.size());
    PtrList<radiation::radiationModel> radiation(fluidRegions.size());
    PtrList<volScalarField> KFluid(fluidRegions.size());
    PtrList<volScalarField> dpdtFluid(fluidRegions.size());

    List<scalar> initialMassFluid(fluidRegions.size());

    PtrList<fv::IOoptionList> fluidFvOptions(fluidRegions.size());

    // Populate fluid field pointer lists
    forAll(fluidRegions, i)
    {
        Info<< "*** Reading fluid mesh thermophysical properties for region "
            << fluidRegions[i].name() << nl << endl;

        Info<< "    Adding to thermoFluid\n" << endl;

        thermoFluid.set
        (
            i,
            rhoThermo::New(fluidRegions[i]).ptr()
        );



        Info<< "    Adding to UFluid\n" << endl;
        UFluid.set
        (
            i,
            new volVectorField
            (
                IOobject
                (
                    "U",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        Info<< "    Adding to phiFluid\n" << endl;
        phiFluid.set
        (
            i,
            new surfaceScalarField
            (
                IOobject
                (
                    "phi",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                linearInterpolate(rhoFluid[i]*UFluid[i])
                    & fluidRegions[i].Sf()
            )
        );
//-----------------------------------------------------------------------------
// here we can try to add the additional lines for 2 phases fluid and interface
    Info<< "Constructing twoPhaseMixtureThermo\n" << endl;
// instead if the original "mesh" from compressibleInterFoam we use fluidRegions[i]     
    twoPhaseMixtureThermo twoPhaseProperties(fluidRegions[i]);
    volScalarField& alpha1(twoPhaseProperties.alpha1());
    volScalarField& alpha2(twoPhaseProperties.alpha2());
    
    Info<< "Reading thermophysical properties\n" << endl;

    volScalarField& p = twoPhaseProperties.p();
    volScalarField& T = twoPhaseProperties.T();
    volScalarField& rho1 = twoPhaseProperties.thermo1().rho();
    const volScalarField& psi1 = twoPhaseProperties.thermo1().psi();
    volScalarField& rho2 = twoPhaseProperties.thermo2().rho();
    const volScalarField& psi2 = twoPhaseProperties.thermo2().psi();
          
//-----------------------------------------------------------------------------

        Info<< "    Adding to rhoFluid\n" << endl;
        rhoFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "rho",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::READ_IF_PRESENT,
                    IOobject::AUTO_WRITE
                ),
                alpha1*rho1 + alpha2*rho2 // thermoFluid[i].rho()
                // does not matter the i counter since in this preliminary trial we will just have one fluid domain with two phases inside.
            )
        );


// commented for debugging
/*
        Info<< "    Adding to gFluid\n" << endl;
        gFluid.set
        (
            i,
            new uniformDimensionedVectorField
            (
                IOobject
                (
                    "g",
                    runTime.constant(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE
                )
            )
        );
*/

//----------------------------------------from compressibleInterFoam
    dimensionedScalar pMin(twoPhaseProperties.lookup("pMin"));
/*
by doing a grep -R "fluidRegions" .
we get fluidRegions[i].C()
and fluidRegions[i].Cf()
therefore is reasonable to subtitute mesh con fluidRegions nelle seguenti righe
*/    


//  commented for debugging
//    Info<< "Calculating field g.h\n" << endl;
//    volScalarField gh("gh", g & fluidRegions[i].C());
//    surfaceScalarField ghf("ghf", g & fluidRegions[i].Cf());
    
    
    
    
//    volScalarField gh("gh", g & mesh.C());
//    surfaceScalarField ghf("ghf", g & mesh.Cf());
// ci da pero' errore 
// ./fluid/createFluidFields.H:137: error: ‘g’ was not declared in this scope
// proviamo ad effettuare un grep -R ' g ' . in compressibleInterFoam per vedere in quale file venga dichiarato.
/*
OK risolto aggiungendo in the main body of interChtMultiRegionFoam.C the line 
    #include "readGravitationalAcceleration.H"
*/  
 
/*
altri errori da correggere 
./fluid/createFluidFields.H:158: error: ‘rho’ was not declared in this scope
./fluid/createFluidFields.H:158: error: ‘phi’ was not declared in this scope
./fluid/createFluidFields.H:167: error: ‘U’ was not declared in this scope
./fluid/createFluidFields.H:194: error: no matching function for call to ‘
si in effetti in createFluidFields.H
rho is called  rhoFluid[i]
phi is called phiFluid[i]
U   is called  UFluid[i]
*/ 
  
    // Mass flux
    // Initialisation does not matter because rhoPhi is reset after the
    // alpha1 solution before it is used in the U equation.
    surfaceScalarField rhoPhi
    (
        IOobject
        (
            "rhoPhi",
            runTime.timeName(),
            fluidRegions[i],
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        fvc::interpolate(rhoFluid[i])*phiFluid[i]
    );

    volScalarField dgdt
    (
        pos(alpha2)*fvc::div(phiFluid[i])/max(alpha2, scalar(0.0001))
    );

    // Construct interface from alpha1 distribution
    interfaceProperties interface(alpha1, UFluid[i], twoPhaseProperties);

    // Construct compressible turbulence model
    autoPtr<compressible::turbulenceModel> turbulence
    (
        compressible::turbulenceModel::New(rhoFluid[i], UFluid[i], rhoPhi, twoPhaseProperties)
    );

    Info<< "Creating field kinetic energy K\n" << endl;
    volScalarField K("K", 0.5*magSqr(UFluid[i]));  
    
//----------------------------------------



        /*
        Info<< "    Adding turbulence \n" << endl;
        turbulence.set
        (
            i,
            compressible::turbulenceModel::New
            (
                rhoFluid[i],
                UFluid[i],
                phiFluid[i],
                thermoFluid[i]
            ).ptr()
        );
        */

        Info<< "    Adding to ghFluid\n" << endl;
        ghFluid.set
        (
            i,
            new volScalarField("gh", gFluid[i] & fluidRegions[i].C())
        );

        Info<< "    Adding to ghfFluid\n" << endl;
        ghfFluid.set
        (
            i,
            new surfaceScalarField("ghf", gFluid[i] & fluidRegions[i].Cf())
        );

        p_rghFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "p_rgh",
                    runTime.timeName(),
                    fluidRegions[i],
                    IOobject::MUST_READ,
                    IOobject::AUTO_WRITE
                ),
                fluidRegions[i]
            )
        );

        // Force p_rgh to be consistent with p
        p_rghFluid[i] = thermoFluid[i].p() - rhoFluid[i]*ghFluid[i];

        radiation.set
        (
            i,
            radiation::radiationModel::New(thermoFluid[i].T())
        );

        initialMassFluid[i] = fvc::domainIntegrate(rhoFluid[i]).value();

        Info<< "    Adding to KFluid\n" << endl;
        KFluid.set
        (
            i,
            new volScalarField
            (
                "K",
                0.5*magSqr(UFluid[i])
            )
        );

        Info<< "    Adding to dpdtFluid\n" << endl;
        dpdtFluid.set
        (
            i,
            new volScalarField
            (
                IOobject
                (
                    "dpdt",
                    runTime.timeName(),
                    fluidRegions[i]
                ),
                fluidRegions[i],
                dimensionedScalar
                (
                    "dpdt",
                    thermoFluid[i].p().dimensions()/dimTime,
                    0
                )
            )
        );

        Info<< "    Adding fvOptions\n" << endl;
        fluidFvOptions.set
        (
            i,
            new fv::IOoptionList(fluidRegions[i])
        );
    }
