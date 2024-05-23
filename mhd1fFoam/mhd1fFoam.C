/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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
    mhdFoam

Description
    Solver for magnetohydrodynamics (MHD): incompressible, laminar flow of a
    conducting fluid under the influence of a magnetic field.

    An applied magnetic field H acts as a driving force,
    at present boundary conditions cannot be set via the
    electric field E or current density J. The fluid viscosity nu,
    conductivity sigma and permeability mu are read in as uniform
    constants.

    A fictitous magnetic flux pressure pH is introduced in order to
    compensate for discretisation errors and create a magnetic face flux
    field which is divergence free as required by Maxwell's equations.

    However, in this formulation discretisation error prevents the normal
    stresses in UB from cancelling with those from BU, but it is unknown
    whether this is a serious error.  A correction could be introduced
    whereby the normal stresses in the discretised BU term are replaced
    by those from the UB term, but this would violate the boundedness
    constraint presently observed in the present numerics which
    guarantees div(U) and div(H) are zero.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"
#include "pisoControl.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char *argv[])
{
    #include "postProcess.H"

    #include "setRootCaseLists.H"
    #include "createTime.H"
    #include "createMesh.H"
    #include "createControl.H"
    #include "createFields.H"
    #include "initContinuityErrs.H"

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
         Info<< "Getting cell centroid coordinates" <<endl;
         int cnt=0;
         forAll(U,cellI)
         {
             cnt++; //Info<< " cellI = " << cellI << endl;
         }
         int Length = cnt;  //sizeof(U)/sizeof(U[0]);
         Info<< " Length = " << Length << endl;
         double x[Length];
         double y[Length];
        // /* 
         for(int cellI=0; cellI<=Length; cellI++ )
	{
          x[cellI] = mesh.C()[cellI].component(0);
          y[cellI] = mesh.C()[cellI].component(1);
         // Info<< "Cell Number: " << cellI << endl; 
          //Info<< " x = " << x[cellI] << endl; 
          //Info<< " y = " << y[cellI] << nl << endl; 
         }
        // */
        Info<< "Makeing momentoun source Si field with cell centroid coordinates" <<endl;
        volScalarField Si("Si", rho); //Si=p;
        forAll(Si, cellI)
        {
           Si[cellI] = 10e4*Foam::exp(
          -10e-4*((0.1+x[cellI])*(0.1+
          x[cellI])+y[cellI]*y[cellI]));
        }
        
        
        
    Info<< nl << "Starting time loop" << endl;

    while (runTime.loop())
    {
        Info<< "Time = " << runTime.timeName() << nl << endl;

        #include "CourantNo.H"

        {
            
            fvVectorMatrix UEqn
            (
                fvm::ddt(U) 
              + fvm::div(phi, U)
              - fvc::div(phiB, 2.0*DBU*B)
              - fvm::laplacian(nu, U) 
              + fvc::grad(DBU*magSqr(B))
            );
            
            if (piso.momentumPredictor())
            {
                
                solve(UEqn == -fvc::grad(p)); //*rho
                //solve(fvm::ddt(rho) == mag(rho*U)); //fvm::div(1.0,rho*U)
            }

	   
	   
            // --- PISO loop
            while (piso.correct())
            {
                volScalarField rAU(1.0/UEqn.A());
                surfaceScalarField rAUf("rAUf", fvc::interpolate(rAU));
                volVectorField HbyA(constrainHbyA(rAU*UEqn.H(), U, p));
                surfaceScalarField phiHbyA
                (
                    "phiHbyA",
                    fvc::flux(HbyA)
                  + rAUf*fvc::ddtCorr(U, phi)
                );

                // Update the pressure BCs to ensure flux consistency
                constrainPressure(p, U, phiHbyA, rAUf);

                while (piso.correctNonOrthogonal())
                {
                    fvScalarMatrix pEqn
                    (
                        fvm::laplacian(rAUf, p) == fvc::div(phiHbyA)
                    );

                    pEqn.setReference(pRefCell, pRefValue);
                    pEqn.solve();

                    if (piso.finalNonOrthogonalIter())
                    {
                        phi = phiHbyA - pEqn.flux();
                    }
                }

                #include "continuityErrs.H"

                U = HbyA - rAU*fvc::grad(p);
                U.correctBoundaryConditions();
            }
        
        }
        // defining the time derivative of rho, continuity equation

        //fvScalarMatrix Cnt (fvm::ddt(rho) + fvm::div(rho,U)); 
        // solving the continuity equation
        //Cnt.relax();
       // Cnt.solve(); 

        // --- B-PISO loop
        while (bpiso.correct())
        {
            fvVectorMatrix BEqn
            (
                fvm::ddt(B)
              + fvm::div(phi, B)
              - fvc::div(phiB, U)
              - fvm::laplacian(DB, B)
            );

            BEqn.solve();

            volScalarField rAB(1.0/BEqn.A());
            surfaceScalarField rABf("rABf", fvc::interpolate(rAB));

            phiB = fvc::flux(B);

            while (bpiso.correctNonOrthogonal())
            {
                fvScalarMatrix pBEqn
                (
                    fvm::laplacian(rABf, pB) == fvc::div(phiB)
                );

                pBEqn.solve();

                if (bpiso.finalNonOrthogonalIter())
                {
                    phiB -= pBEqn.flux();
                }
            }

            #include "magneticFieldErr.H"
        }

        runTime.write();
    }

    Info<< "End\n" << endl;

    return 0;
}


// ************************************************************************* //
