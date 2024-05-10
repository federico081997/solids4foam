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
    Implicit functions for vertex-centred finite volume calculations: Vertex
    Finite Volume Matrix (vfvm).

    Populates the matrix for divergence of stress on a dual mesh, where the
    dual mesh cells corresponds to the primary mesh points.

    The Laplacian coefficients are also calculated.

    The d2dt2 operator is also defined.

    Note that these functions only set the matrix coefficients and they do not
    modify the right hand side.

SourceFiles
    vfvmCellPoint.C

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

\*---------------------------------------------------------------------------*/

#ifndef vfvmCellPoint_H
#define vfvmCellPoint_H

#include "volFields.H"
#include "pointFields.H"
#include "sparseMatrix.H"
#include "sparseMatrixExtended.H"
#include "RectangularMatrix.H"
#include "SparseMatrixTemplate.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace vfvm
{
    // Add coefficients to the matrix for the divergence of stress in the
    // small-strain displacement formulation. Each coefficient in the matrix
    // is a 3x3 tensor.
    void divSigma
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug = false
    );

    // Add coefficients to the matrix for the divergence of stress in the
    // large-strain displacement formulation. Each coefficient in the matrix
    // is a 3x3 tensor.
    void divSigma
    (
        sparseMatrix& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const Field<scalarRectangularMatrix>& geometricStiffnessField,
        const symmTensorField& sigma,
        const tensorField& dualGradDField,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug = false
    );


    // Add coefficients to the matrix for the divergence of stress in the
    // small-strain displacement-pressure formulation. Each coefficient in the
    // matrix is a 4x4 tensor.
    void divSigma
    (
        sparseMatrixExtended& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug
    );


    // Add coefficients to the matrix for the divergence of stress in the
    // large-strain displacement-pressure formulation. Each coefficient in the
    // matrix is a 4x4 tensor.
    void divSigma
    (
        sparseMatrixExtended& matrix,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const Field<scalarSquareMatrix>& materialTangentField,
        const Field<scalarRectangularMatrix>& geometricStiffnessField,
        const symmTensorField& sigma,
        const tensorField& dualGradDField,
        const scalar zeta, // Fraction of compact edge direction component used
        const bool debug = false
    );


    // Add coefficients to the matrix for the Laplacian of a scalar field with
    // a scalar diffusivity field. Each coefficient in the matrix is a scalar.
    void laplacian
    (
        sparseScalarMatrix& matrix,
        const Switch compactStencil,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalarField& diffusivity, // Diffusivity in the primary cells
        const bool debug = false
    );


    // Add coefficients to the matrix for the Laplacian of a scalar field with
    // a scalar diffusivity in the small-strain pressure-displacement
    // formulation. Each coefficient in the matrix is a 4x4 tensor.
    void laplacian
    (
        sparseMatrixExtended& matrix,
        const Switch compactStencil,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const scalar& diffusivity,
        const bool debug
    );


    // Add coefficients to the matrix for the Laplacian of a scalar field with
    // a scalar diffusivity in the large-strain pressure-displacement
    // formulation. Each coefficient in the matrix is a 4x4 tensor.
    void laplacian
    (
        sparseMatrixExtended& matrix,
        const Switch compactStencil,
        const fvMesh& mesh,
        const fvMesh& dualMesh,
        const labelList& dualFaceToCell,
        const labelList& dualCellToPoint,
        const tensorField& dualGradDField,
        const scalar& diffusivity,
        const bool debug
    );


    // Add coefficients to the matrix for the second time derivative. Each
    // coefficient in the matrix is a 3x3 tensor.
    void d2dt2
    (
        ITstream& d2dt2Scheme,
        const scalar& deltaT,           // Time-step
        const word& pointDname,
        sparseMatrix& matrix,
        const scalarField& pointRhoI,
        const scalarField& pointVolI,
        const int debug                 // Debug switch
    );


    // Add source coefficients to the matrix diagonal of the pressure
    // equation. Each coefficient in the matrix is a 4x4 tensor.
    void Sp
    (
        sparseMatrixExtended& matrix,
        const scalarField& pointVolI,
        const int debug
    );


    // Add the displacement coefficients of the pressure equation to the
    // matrix. Each coefficients in the matrix is a 4x4 tensor.
    void addPbar
    (
        sparseMatrixExtended& matrix,
        const fvMesh& dualMesh,
        const labelList& dualCellToPoint,
        const scalarField& pointVolI,
        const tensorField& pBarSensitivity,
        const int debug
    );


} // End namespace vfvc

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
