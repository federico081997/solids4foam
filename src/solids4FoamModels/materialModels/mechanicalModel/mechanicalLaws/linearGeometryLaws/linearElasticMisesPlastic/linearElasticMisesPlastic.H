/*---------------------------------------------------------------------------*\
License
    This file is part of solids4foam.

    solids4foam is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    solids4foam is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with solids4foam.  If not, see <http://www.gnu.org/licenses/>.

Class
    linearElasticMisesPlastic

Description
    Elasto-plastic mechanical law as described by Simo & Hughes (1998)
    in Box 3.2.

    Linear (Hookean) elasticity with Mises/J2 plasticity.

    For elastic parameters, the user can specify:
        - Young's modulus (E) and Poisson's ratio (nu)
    or
        - Shear modulus (mu) and bulk modulus (K)

    More details found in:

    Simo & Hughes, Computational Inelasticity, 1998, Springer.

    P. Cardiff, Z. Tukovic, P. De Jaeger, M. Clancy and A. Ivankovic. A
    Lagrangian cell-centred finite volume method for metal forming simulation,
    doi=10.1002/nme.5345.

SourceFiles
    linearElasticMisesPlastic.C

Author
    Philip Cardiff, UCD. All rights reserved.

\*---------------------------------------------------------------------------*/

#ifndef linearElasticMisesPlastic_H
#define linearElasticMisesPlastic_H

#include "mechanicalLaw.H"
#include "surfaceMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "interpolationTable.H"
#include "pointFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                  Class linearElasticMisesPlastic Declaration
\*---------------------------------------------------------------------------*/

class linearElasticMisesPlastic
:
    public mechanicalLaw
{
    // Private data

        //- Point mesh
        const pointMesh& pMesh_;

        //- Shear modulus
        dimensionedScalar mu_;

        //- Bulk modulus
        dimensionedScalar K_;

        //- Young's modulus
        dimensionedScalar E_;

        //- Poisson's ratio
        dimensionedScalar nu_;

        //- Table of plastic strain versus post-yield stress
        interpolationTable<scalar> stressPlasticStrainSeries_;

        //- Yield stress
        volScalarField sigmaY_;

        //- Yield stress surface field
        surfaceScalarField sigmaYf_;

        //- Yield stress point field
        pointScalarField pSigmaY_;

        //- Incremental change in sigmaY
        volScalarField DSigmaY_;

        //- Incremental change in sigmaY surface field
        surfaceScalarField DSigmaYf_;

        //- Incremental change in sigmaY point field
        pointScalarField pDSigmaY_;

        //- Total plastic strain
        volSymmTensorField epsilonP_;

        //- Total plastic strain surface field
        surfaceSymmTensorField epsilonPf_;

        //- Total plastic strain point field
        pointSymmTensorField pEpsilonP_;

        //- Incremental change of plastic strain
        volSymmTensorField DEpsilonP_;

        //- Incremental change of plastic strain surface field
        surfaceSymmTensorField DEpsilonPf_;

        //- Equivalent plastic strain increment
        volScalarField DEpsilonPEq_;

        //- Equivalent plastic strain increment surface field
        surfaceScalarField DEpsilonPEqf_;

        //- Plastic multiplier increment - plastric strain scaled by sqrt(2/3)
        volScalarField DLambda_;

        //- Plastic multiplier increment - surface field
        surfaceScalarField DLambdaf_;

        //- Plastic multiplier increment - point field
        pointScalarField pDLambda_;

        //- Equivalent plastic strain
        volScalarField epsilonPEq_;

        //- Equivalent plastic strain surface field
        surfaceScalarField epsilonPEqf_;

        //- Equivalent plastic strain point field
        pointScalarField pEpsilonPEq_;

        //- Active yielding flag
        //     1.0 for active yielding
        //     0.0 otherwise
        volScalarField activeYield_;

        //- plasticN is the return direction to the yield surface
        volSymmTensorField plasticN_;

        //- plasticN surface field
        surfaceSymmTensorField plasticNf_;

        //- plasticN point field
        pointSymmTensorField pPlasticN_;

        //- An iterative procedure is used when the yield stress is a nonlinear
        //  function of plastic strain
        const bool nonLinearPlasticity_;

        //- Linear plastic modulus. It is only used if plasticity is linear,
        //  defined by two points on the stress-plastic strain curve
        scalar Hp_;

        //- Maximum allowed error in the plastic strain integration
        const scalar maxDeltaErr_;

        //- Tolerance for Newton loop
        static scalar LoopTol_;

        //- Maximum number of iterations for Newton loop
        static label MaxNewtonIter_;

        //- finiteDiff is the delta for finite difference differentiation
        static scalar finiteDiff_;

        //- Store sqrt(2/3) as it is used often
        static scalar sqrtTwoOverThree_;


    // Private Member Functions

        //- Disallow default bitwise copy construct
        linearElasticMisesPlastic(const linearElasticMisesPlastic&);

        //- Disallow default bitwise assignment
        void operator=(const linearElasticMisesPlastic&);

        //- Update plasticN, DLambda, DSigmaY and sigmaY
        void updatePlasticity
        (
            symmTensor& plasticN,          // Plastic return direction
            scalar& DLambda,               // Plastic multiplier increment
            scalar& DSigmaY,               // Increment of yield stress
            scalar& sigmaY,                // Yield stress
            const scalar sigmaYOld,        // Yield stress old time
            const scalar fTrial,           // Trial yield function
            const symmTensor& sTrial,      // Trial deviatoric stress
            const scalar epsilonPEqOld,    // Old equivalent plastic strain
            const scalar muBar,            // Scaled shear modulus
            const scalar maxMagDEpsilon    // Max strain increment magnitude
        ) const;

        //- Return the current yield stress
        scalar curYieldStress
        (
            const scalar curEpsilonPEq     // Current equivalent plastic strain
        ) const;

        //- Evaulate current value of the yield function
        scalar yieldFunction
        (
            const scalar epsilonPEqOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar DLambda,          // Plastic multiplier
            const scalar muBar             // Scaled shear modulus
        ) const;

        //- Iteratively calculate plastic multiplier increment (DLambda)
        //  and current yield stress using Newton's method
        void newtonLoop
        (
            scalar& DLambda,               // Plastic multiplier
            scalar& curSigmaY,             // Current yield stress
            const scalar epsilonPEqOld,    // Old equivalent plastic strain
            const scalar magSTrial,        // Deviatoric trial stress magnitude
            const scalar muBar,            // Scaled shear modulus
            const scalar maxMagDEpsilon    // Max strain increment magnitude
        ) const;

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            volScalarField& sigmaHyd,
            const volScalarField& trEpsilon
        );

        //- Calculate hydrostatic component of the stress tensor
        void calculateHydrostaticStress
        (
            surfaceScalarField& sigmaHyd,
            const surfaceScalarField& trEpsilon
        );

        //- Calculate stress given the strain
        void calculateStress
        (
            surfaceSymmTensorField& sigma,
            const surfaceTensorField& gradD,
            const surfaceScalarField& p
        ) const;


public:

    //- Runtime type information
    TypeName("linearElasticMisesPlastic");

    // Constructors

        //- Construct from dictionary
        linearElasticMisesPlastic
        (
            const word& name,
            const fvMesh& mesh,
            const dictionary& dict,
            const nonLinearGeometry::nonLinearType& nonLinGeom
        );

        //- Calculate stress at the primary points
        void calculatePStress
        (
            pointSymmTensorField& pSigma,
            const pointTensorField& pGradD
        ) const;


    // Destructor

        virtual ~linearElasticMisesPlastic();


    // Member Functions

        //- Return the implicit stiffness
        //  This is the diffusivity for the Laplacian term
        virtual tmp<volScalarField> impK() const;

        //- Return bulk modulus at the points
        virtual tmp<pointScalarField> pBulkModulus() const;

#ifdef OPENFOAM_COM
        //- Material tangent for Newton Raphson loop
        virtual tmp<Field<scalarSquareMatrix>> materialTangentField() const;
#endif

        //- Return the implicit stiffness as a diagTensor
#ifndef OPENFOAM_NOT_EXTEND
        virtual tmp<volDiagTensorField> impKdiagTensor() const;
#endif

        //- Update the volField stress
        virtual void correct(volSymmTensorField& sigma);

        //- Update the surfaceField stress
        virtual void correct(surfaceSymmTensorField& sigma);

        //- Return material residual i.e. a measured of how convergence of
        //  the material model
        virtual scalar residual();

        //- Update the yield stress: called at end of time-step
        virtual void updateTotalFields();

        //- Return the desired new time-step
        virtual scalar newDeltaT();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
