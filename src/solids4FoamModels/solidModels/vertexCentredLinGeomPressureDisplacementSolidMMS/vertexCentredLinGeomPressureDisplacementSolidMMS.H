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
    vertexCentredLinGeomPressureDisplacementSolidMMS

Description
    Mathematical model where linear geometry is assumed i.e. small strains and
    small rotations are assumed, and the total displacement is the primary
    unknown.

    Discretisation is performed using a vertex-centred finite volume approach
    where the governing equations are integrated over a dual mesh. The dual
    mesh is automatically constructed from the primary mesh using the
    meshDualiser.

Author
    Philip Cardiff, UCD.  All rights reserved.
    Federico Mazzanti, UCD. All rights reserved.

SourceFiles
    vertexCentredLinGeomPressureDisplacementSolidMMS.C

\*---------------------------------------------------------------------------*/

#ifndef vertexCentredLinGeomPressureDisplacementSolidMMS_H
#define vertexCentredLinGeomPressureDisplacementSolidMMS_H

#include "solidModel.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "uniformDimensionedFields.H"
#include "GeometricField.H"
#include "dualMechanicalModel.H"
#include "globalPointIndices.H"
#include "pointFieldFunctions.H"
#include "sparseMatrixExtended.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

/*---------------------------------------------------------------------------*\
        Class vertexCentredLinGeomPressureDisplacementSolidMMS Declaration
\*---------------------------------------------------------------------------*/

class vertexCentredLinGeomPressureDisplacementSolidMMS
:
    public solidModel
{
    // Private Data

        //- Mechanical model for dual mesh
        autoPtr<dualMechanicalModel> dualMechanicalPtr_;

        //- Flag to enable steadyState
        Switch steadyState_;

        //- Compact stencil flag
        Switch compactStencil_;

        //- zeta coefficient
        const scalar zeta_;

        //- Pressure smoothing coefficient
        const scalar pressureSmoothingFactor_;

        //- Bulk modulus at the primary points
        const pointScalarField pointK_;

        //- Flag indicating if the case is 2-D or 3-D
        const bool twoD_;

        //- Field of flags indicating the degrees of freedom where the value is
        //  set
        //  Currently the size of this list is equal to the number of points
        boolList fixedDofs_;

        //- Values for fixed degrees of freedom
        //  Currently the size of this list is equal to the number of points
        pointField fixedDofValues_;

        //- Fixed directions for fixed degrees of freedom
        //  Calculated as sqr(n) where n is a fixed direction, e.g. I for fully
        //  fixed
        symmTensorField fixedDofDirections_;

        //- Scale factor for the fixed degree of freedom equations to improve
        //  the matrix conditioning
        const scalar fixedDofScale_;

        //- Point pressure
        pointScalarField pointP_;

        //- Point velocity
        pointVectorField pointU_;

        //- Point acceleration field
        pointVectorField pointA_;

        //- Point density field
        pointScalarField pointRho_;

        //- Point volume field
        //  Dual mesh cell volumes associated with the primary mesh points
        pointScalarField pointVol_;

        //- Point gradD field
        //  gradD associated with the primary mesh points
        pointTensorField pointGradD_;

        //- Displacement gradient field at the dual mesh faces
        surfaceTensorField dualGradDf_;

        //- Stress field at the dual mesh faces
        surfaceSymmTensorField dualSigmaf_;

        //- Pressure at the dual mesh faces
        surfaceScalarField dualPf_;

        //- Cell pressure
        volScalarField volP_;

        //- Local-to-global point map and owner list
        globalPointIndices globalPointIndices_;


    // Private Member Functions

        //- Set fixed degree of freedom lists
        void setFixedDofs
        (
            const pointVectorField& pointD,
            boolList& fixedDofs,
            pointField& fixedDofValues,
            symmTensorField& fixedDofDirections
        ) const;

        //- Enforce extact tractions on traction boundaries
        void enforceTractionBoundaries
        (
            const pointVectorField& pointD,
            surfaceVectorField& dualTraction,
            const fvMesh& mesh, // Corresponding to pointD mesh
            const labelListList& pointToDualFaces
        ) const;

        //- Check convergence of Newton-Raphson loop
        bool converged
        (
            const label iCorr,
            scalar& initResidualD,
            scalar& initResidualP,
            const label nInterations,
            const pointVectorField& pointD,
            const pointScalarField& pointP,
            const Field<scalarRectangularMatrix>& pointDPcorr
        ) const;

        //- Calculate the residual of the momentum equation
        tmp<vectorField> residualD
        (
            const pointVectorField& pointD,
            const pointScalarField& pointP
        ) const;

        //- Calculate the residual of the pressure equation
        tmp<scalarField> residualP
        (
            const pointVectorField& pointD,
            const pointScalarField& pointP
        ) const;

        //- Calculate the matrix coefficients using finite difference
        void finiteDiffMatrix
        (
            sparseMatrixExtended& matrix,
            const vectorField& residualD,
            const scalarField& residualP,
            const pointVectorField& pointD,
            const pointScalarField& pointP
        );

        //- Disallow default bitwise copy construct
        vertexCentredLinGeomPressureDisplacementSolidMMS
        (
            const vertexCentredLinGeomPressureDisplacementSolidMMS&
        );

        //- Disallow default bitwise assignment
        void operator=(const vertexCentredLinGeomPressureDisplacementSolidMMS&);

protected:

    // Protected member functions

        //- Return nonlinear geometry enumerator
        virtual nonLinearGeometry::nonLinearType nonLinGeom() const
        {
            return nonLinearGeometry::LINEAR_GEOMETRY;
        }


public:

    //- Runtime type information
    TypeName("vertexCentredPressureDisplacementLinGeometryMMS");

    // Constructors

        //- Construct from components
        vertexCentredLinGeomPressureDisplacementSolidMMS
        (
            Time& runTime,
            const word& region = dynamicFvMesh::defaultRegion
        );

    // Destructor

        virtual ~vertexCentredLinGeomPressureDisplacementSolidMMS();


    // Member Functions

        // Access

            //- Each solidModel must indicate whether D or DD is the primary
            //  solution variable
            virtual volVectorField& solutionD()
            {
                notImplemented("For vertex-centred approach");

                // This model solves for D
                return D();
            }

        // Edit

            //- Evolve the solid solver and solve the mathematical model
            virtual bool evolve();

            //- Traction boundary surface normal gradient
            virtual tmp<vectorField> tractionBoundarySnGrad
            (
                const vectorField& traction,
                const scalarField& pressure,
                const fvPatch& patch
            ) const
            {
                notImplemented("For vertex-centred approach");

                // Keep compiler happy
                return tmp<vectorField>(new vectorField());
            }

            //- Set traction at specified patch
            virtual void setTraction
            (
                const label interfaceI,
                const label patchID,
                const vectorField& faceZoneTraction
            );

        // Write

            //- Write fields
            virtual void writeFields(const Time& runTime);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModel

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
