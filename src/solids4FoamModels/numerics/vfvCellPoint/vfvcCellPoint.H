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

Description
    Explicit functions for vertex-centred finite volume calculations: Vertex
    Finite Volume Calculus (vfvc).

SourceFiles
    vfvcCellPoint.C

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

\*---------------------------------------------------------------------------*/

#ifndef vfvcCellPoint_H
#define vfvcCellPoint_H

#include "volFields.H"
#include "surfaceFields.H"
#include "pointFields.H"
#include "vectorList.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vfvc
{
    // Gradient field
    // Returns volVectorField on primary mesh
    tmp<volVectorField> grad
    (
        const pointScalarField& pointT, // Primary mesh point field
        const fvMesh& mesh              // Primary mesh
    );

    // Gradient field
    // Returns volTensorField on primary mesh
    tmp<volTensorField> grad
    (
        const pointVectorField& pointD, // Primary mesh point displacement
        const fvMesh& mesh              // Primary mesh
    );

    // Gradient field
    // Returns pointTensorField on primary mesh
    tmp<pointTensorField> pGrad
    (
        const pointVectorField& pointD, // Primary mesh point displacement
        const fvMesh& mesh              // Primary mesh
    );

    // Gradient field
    // Returns surfaceVectorField on dual mesh
    // Only calculates the gradient for internal faces
    tmp<surfaceVectorField> fGrad
    (
        const pointScalarField& pointP, // Primary mesh point pressure
        const fvMesh& mesh,             // Primary mesh
        const fvMesh& dualMesh,         // Dual mesh
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug = false
    );

    // Gradient field
    // Returns surfaceTensorField on dual mesh
    // Only calculates the gradient for internal faces
    tmp<surfaceTensorField> fGrad
    (
        const pointVectorField& pointD, // Primary mesh point displacement
        const fvMesh& mesh,             // Primary mesh
        const fvMesh& dualMesh,         // Dual mesh
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug = false
    );

    // Second time derivative pointVectorField pointD
    tmp<vectorField> d2dt2
    (
        ITstream& d2dt2Scheme,
        const pointVectorField& pointD, // Displacement
        const pointVectorField& pointU, // Velocity
        const pointVectorField& pointA, // Acceleration
        const scalarField& pointRho,    // Density
        const scalarField& pointVol,    // Volumes
        const int debug                 // Debug switch
    );

    // First time derivative pointVectorField pointD
    tmp<vectorField> ddt
    (
        ITstream& ddtScheme,
        ITstream& d2dt2Scheme,
        const pointVectorField& pointP // Displacement
    );

    // Returns the laplacian of a scalar field on the dual mesh
    tmp<pointScalarField> laplacian
    (
        const pointScalarField& pointP,   // Pressure
        const fvMesh& mesh,               // Primary mesh
        const fvMesh& dualMesh,           // Dual mesh
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalar& zeta, // Fraction of compact edge direction used
        const bool debug
    );

    // Interpolate a pointScalarField to a volScalarField
    tmp<volScalarField> interpolate
    (
        const pointScalarField& pointP,
        const fvMesh& mesh
    );

    // Interpolate a pointScalarField to a surfaceScalarField
    // at the dual mesh faces
    tmp<surfaceScalarField> interpolate
    (
        const pointScalarField& pointP,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint
    );
    
} // End namespace vfvc

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
