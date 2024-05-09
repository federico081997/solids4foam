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

\*---------------------------------------------------------------------------*/

#include "vertexCentredLinGeomPressureDisplacementSolid.H"
#include "addToRunTimeSelectionTable.H"
#include "vfvcCellPoint.H"
#include "vfvmCellPoint.H"
#include "fvcDiv.H"
#include "fixedValuePointPatchFields.H"
#include "solidTractionPointPatchVectorField.H"
#include "sparseMatrixExtendedTools.H"
#include "symmetryPointPatchFields.H"
#include "fixedDisplacementZeroShearPointPatchVectorField.H"
#ifdef USE_PETSC
    #include <petscksp.h>
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace solidModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(vertexCentredLinGeomPressureDisplacementSolid, 0);
addToRunTimeSelectionTable
(
    solidModel,
    vertexCentredLinGeomPressureDisplacementSolid,
    dictionary
);


// * * * * * * * * * * *  Private Member Functions * * * * * * * * * * * * * //

void vertexCentredLinGeomPressureDisplacementSolid::setFixedDofs
(
    const pointVectorField& pointD,
    boolList& fixedDofs,
    pointField& fixedDofValues,
    symmTensorField& fixedDofDirections
) const
{
    // Flag all fixed DOFs
    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            isA<fixedValuePointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointD.mesh().mesh().boundaryMesh()[patchI].meshPoints();

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];
                const vector& disp = pointD[pointID];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the new one
                    if
                    (
                        mag
                        (
                            fixedDofDirections[pointID]
                          & (fixedDofValues[pointID] - disp)
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                                "void vertexCentredLinGeomPressureDisplacementSolid::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // Set all directions as fixed, just in case it was
                    // previously marked as a symmetry point
                    fixedDofDirections[pointID] = symmTensor(I);
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = disp;
                    fixedDofDirections[pointID] = symmTensor(I);
                }
            }
        }
        else if
        (
            isA<symmetryPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const labelList& meshPoints =
                pointD.mesh().boundary()[patchI].meshPoints();
            const vectorField& pointNormals =
                pointD.mesh().boundary()[patchI].pointNormals();

            scalarField normalDisp(meshPoints.size(), 0.0);
            if
            (
                isA<fixedDisplacementZeroShearPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                )
            )
            {
                normalDisp =
                (
                    pointNormals
                  & pointD.boundaryField()[patchI].patchInternalField()
                );

                if (debug)
                {
                    Info<< "normalDisp = " << normalDisp << endl;
                }
            }

            forAll(meshPoints, pI)
            {
                const label pointID = meshPoints[pI];

                // Check if this point has already been fixed
                if (fixedDofs[pointID])
                {
                    // Check if the existing prescribed displacement is
                    // consistent with the current condition
                    if
                    (
                        mag
                        (
                            (pointNormals[pI] & fixedDofValues[pointID])
                          - normalDisp[pI]
                        ) > SMALL
                    )
                    {
                        FatalErrorIn
                        (
                            "void vertexCentredLinGeomPressureDisplacementSolid::setFixedDofs(...)"
                        )   << "Inconsistent displacements prescribed at point "
                            << "= " << pointD.mesh().mesh().points()[pointID]
                            << abort(FatalError);
                    }

                    // If the point is not fully fixed then make sure the normal
                    // direction is fixed
                    if (mag(fixedDofDirections[pointID] - symmTensor(I)) > 0)
                    {
                        // If the directions are orthogonal we can add them
                        const symmTensor curDir = sqr(pointNormals[pI]);
                        if (mag(fixedDofDirections[pointID] & curDir) > 0)
                        {
                            FatalError
                                << "Point " << pointID << " is fixed in two "
                                << "directions: this is only implemented for "
                                << "Cartesian axis directions"
                                << abort(FatalError);
                        }

                        fixedDofDirections[pointID] += curDir;
                    }
                }
                else
                {
                    fixedDofs[pointID] = true;
                    fixedDofValues[pointID] = normalDisp[pI]*pointNormals[pI];
                    fixedDofDirections[pointID] = sqr(pointNormals[pI]);
                }
            }
        }
    }
}


void vertexCentredLinGeomPressureDisplacementSolid::enforceTractionBoundaries
(
    const pointVectorField& pointD,
    surfaceVectorField& dualTraction,
    const fvMesh& mesh,
    const labelListList& pointToDualFaces
) const
{
    const pointMesh& pMesh = pointD.mesh();
    const fvMesh& dualMesh = dualTraction.mesh();

    forAll(pointD.boundaryField(), patchI)
    {
        if
        (
            isA<solidTractionPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            const solidTractionPointPatchVectorField& tracPatch =
                refCast<const solidTractionPointPatchVectorField>
                (
                    pointD.boundaryField()[patchI]
                );

            const labelList& meshPoints =
                mesh.boundaryMesh()[patchI].meshPoints();

            // Primary mesh point normals
            const vectorField& n =
                pMesh.boundary()[patchI].pointNormals();

            // Primary mesh point tractions
            const vectorField totalTraction
            (
                tracPatch.traction() - n*tracPatch.pressure()
            );

            // Create dual mesh faces traction field
            vectorField dualFaceTraction
            (
                dualMesh.boundaryMesh()[patchI].size(), vector::zero
            );

            // Multiple points map to each dual face so we will count them
            // and then divide the dualFaceTraction by this field so that it is
            // the average of all the points that map to it
            scalarField nPointsPerDualFace(dualFaceTraction.size(), 0.0);

            // Map from primary mesh point field to second mesh face field
            // using the pointToDualFaces map
            forAll(totalTraction, pI)
            {
                const label pointID = meshPoints[pI];
                const labelList& curDualFaces = pointToDualFaces[pointID];

                forAll(curDualFaces, dfI)
                {
                    const label dualFaceID = curDualFaces[dfI];

                    if (!dualMesh.isInternalFace(dualFaceID))
                    {
                        // Check which patch this dual face belongs to
                        const label dualPatchID =
                            dualMesh.boundaryMesh().whichPatch(dualFaceID);

                        if (dualPatchID == patchI)
                        {
                            // Find local face index
                            const label localDualFaceID =
                                    dualFaceID
                              - dualMesh.boundaryMesh()[dualPatchID].start();

                            // Set dual face traction
                            dualFaceTraction[localDualFaceID] +=
                                    totalTraction[pI];

                            // Update the count for this face
                            nPointsPerDualFace[localDualFaceID]++;
                        }
                    }
                }
            }

            if (gMin(nPointsPerDualFace) < 1)
            {
                FatalErrorIn
                (
                    "void vertexCentredLinGeomPressureDisplacementSolid::"
                    "enforceTractionBoundaries(...)"
                )   << "Problem setting tractions: gMin(nPointsPerDualFace) < 1"
                    << nl << "nPointsPerDualFace = "
                    << nPointsPerDualFace << abort(FatalError);
            }

            // Take the average
            dualFaceTraction /= nPointsPerDualFace;

            // Overwrite the dual patch face traction
#ifdef OPENFOAM_NOT_EXTEND
            dualTraction.boundaryFieldRef()[patchI] = dualFaceTraction;
#else
            dualTraction.boundaryField()[patchI] = dualFaceTraction;
#endif
        }
        else if
        (
            isA<symmetryPointPatchVectorField>(pointD.boundaryField()[patchI])
         || isA<fixedDisplacementZeroShearPointPatchVectorField>
            (
                pointD.boundaryField()[patchI]
            )
        )
        {
            // Set the dual patch face shear traction to zero
            const vectorField n(dualMesh.boundary()[patchI].nf());
#ifdef OPENFOAM_NOT_EXTEND
            dualTraction.boundaryFieldRef()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#else
            dualTraction.boundaryField()[patchI] =
                (sqr(n) & dualTraction.boundaryField()[patchI]);
#endif
        }
    }
}


bool vertexCentredLinGeomPressureDisplacementSolid::converged
(
    const label iCorr,
    scalar& initResidualD,
    scalar& initResidualP,
    const label nInterations,
    const pointVectorField& pointD,
    const pointScalarField& pointP,
    const Field<RectangularMatrix<scalar>>& pointDPcorr
) const
{
    scalar residualDAbs = 0;
    scalar residualPAbs = 0;

    // Calculate the residuals as the root mean square of the correction
    scalar maxDCorr = 0.0;
    scalar maxPCorr = 0.0;
    forAll(pointDPcorr, pointI)
    {
        // Displacement residual
        const vector curDCorr
        (
            pointDPcorr[pointI](0,0),
            pointDPcorr[pointI](1,0),
            pointDPcorr[pointI](2,0)
        );

        residualDAbs += magSqr(curDCorr);
        maxDCorr = max(maxDCorr, mag(curDCorr));

        // Pressure residual
        residualPAbs += sqr(pointDPcorr[pointI](3,0));
        maxPCorr = max(maxPCorr, mag(pointDPcorr[pointI](3,0)));
    }

    residualDAbs /= sqrt(residualDAbs/pointD.size());
    residualPAbs /= sqrt(residualPAbs/pointD.size());

    // Store initial residual
    if (iCorr == 0)
    {
        initResidualD = residualDAbs;
        initResidualP = residualPAbs;

        // If the initial residual is small then convergence has been achieved
        if (initResidualD < SMALL && initResidualP < SMALL)
        {
            Info<< "    Both displacement and pressure residuals are less "
                << "than 1e-15"
                << "    Converged" << endl;
            return true;
        }
        Info<< "    Initial displacement residual = " << initResidualD << endl;
        Info<< "    Initial pressure residual = " << initResidualP << endl;
    }

    // Define a normalised residual wrt the initial residual
    const scalar residualDNorm = residualDAbs/initResidualD;
    const scalar residualPNorm = residualPAbs/initResidualP;

    // Calculate the maximum displacement
    const scalar maxMagD = gMax(mag(pointD.primitiveField()));

    // Calculate the maximum pressure
    const scalar maxMagP = gMax(mag(pointP.primitiveField()));

    // Print information for the displacement
    Info<< "    Iter = " << iCorr
        << ", relRef = " << residualDNorm
        << ", resAbs = " << residualDAbs
        << ", nIters = " << nInterations
        << ", maxD = " << maxMagD
        << ", maxDCorr = " << maxDCorr << endl;

    // Print information for the pressure
    Info<< "    Iter = " << iCorr
        << ", relRef = " << residualPNorm
        << ", resAbs = " << residualPAbs
        << ", nIters = " << nInterations
        << ", maxP = " << maxMagP
        << ", maxPCorr = " << maxPCorr << endl;

    // Displacement tolerance
    const scalar DTol =
        solidModelDict().lookupOrDefault<scalar>("solutionDTolerance", 1e-11);

    // Pressure tolerance
    const scalar PTol =
        solidModelDict().lookupOrDefault<scalar>("solutionPTolerance", 1e-6);

    // Check for convergence
    if (residualDNorm < DTol && residualPNorm < PTol)
    {
        Info<< "    Converged" << endl;
        return true;
    }
    else if (iCorr >= nCorr() - 1)
    {
        if (nCorr() > 1)
        {
            Warning
                << "Max iterations reached within the momentum Newton-Raphson "
                "loop" << endl;
        }

        return true;
    }

    // Convergence has not been reached
    return false;
}


tmp<vectorField>
vertexCentredLinGeomPressureDisplacementSolid::residualD
(
    const pointVectorField& pointD,
    const pointScalarField& pointP
) const
{
    // Prepare the result
    tmp<vectorField> tresult(new vectorField(pointD.size(), vector::zero));
    vectorField& result = tresult.ref();

    // The momentum residual (residualD) vector is
    // F = div(sigma) + rho*g - rho*d2dt2(D)
    //   = div(s - p*I) + rho*g - rho*d2dt2(D)

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Point density field
    const scalarField& pointRhoI = pointRho_.internalField();

    // Dual face area vectors
    const surfaceVectorField dualSf = dualMesh().Sf();

    // Magnitude of dualSf
    const surfaceScalarField dualMagSf(mag(dualSf));

    // Dual face unit normals
    const surfaceVectorField dualN(dualSf/dualMagSf);

    // Calculate the Cauchy tractions on the dual faces
    surfaceVectorField dualTraction(dualN & dualSigmaf_);

    // Enforce exact tractions on traction boundaries
    enforceTractionBoundaries
    (
        pointD, dualTraction, mesh(), dualMeshMap().pointToDualFaces()
    );

    // Set coupled boundary (e.g. processor) traction fields to zero: this
    // ensures their global contribution is zero
    forAll(dualTraction.boundaryField(), patchI)
    {
        if (dualTraction.boundaryField()[patchI].coupled())
        {
            dualTraction.boundaryFieldRef()[patchI] = vector::zero;
        }
    }

    // Calculate the divergence of stress for the dual cells
    const vectorField dualDivSigma = fvc::div(dualTraction*dualMagSf);

    // Map dual cell field to primary mesh point field
    vectorField pointDivSigma(mesh().nPoints(), vector::zero);
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualDivSigma, dualCellI)
    {
        const label pointID = dualCellToPoint[dualCellI];
        pointDivSigma[pointID] = dualDivSigma[dualCellI];
    }

    // Add surface forces
    result += pointDivSigma*pointVolI;

    // Add gravity body forces
    result += pointRhoI*g().value()*pointVolI;

    // Add transient term
    result -= vfvc::d2dt2
    (
        mesh().d2dt2Scheme("d2dt2(pointD)"),
        pointD,
        pointU_,
        pointA_,
        pointRho_,
        pointVol_,
        int(bool(debug))
    );

    // Return the momentum residual field
    return tresult;
}


tmp<scalarField>
vertexCentredLinGeomPressureDisplacementSolid::residualP
(
    const pointVectorField& pointD,
    const pointScalarField& pointP
) const
{
    // Prepare the result
    tmp<scalarField> tresult(new scalarField(pointD.size(), 0));
    scalarField& result = tresult.ref();

    // The residual for the pressure equation is:
    // F = p - gamma*laplacian(p) - pBar(D)

    // Calculate gradD at the primary mesh points
    const pointTensorField pointGradD
    (
        vfvc::pGrad
        (
            pointD,
            mesh()
        )
    );

    // Calculate laplacian(p) field
    const pointScalarField laplacianP
    (
        vfvc::laplacian
        (
            pointP,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta_,
            int(bool(debug))
        )
    );

    // Calculate the pBar field
    const scalarField pBar(-pointK_.internalField()*tr(pointGradD));

    // Point volume field
    const scalarField& pointVolI = pointVol_.internalField();

    // Add pressure
    result += pointP*pointVolI;

    // Add gamma*laplacian(p)
    result -= pressureSmoothingFactor_*laplacianP;

    // Add pBar
    result -= pBar*pointVolI;

    // Return the residual field
    return tresult;
}


void vertexCentredLinGeomPressureDisplacementSolid::finiteDiffMatrix
(
    sparseMatrixExtended& matrix,
    const vectorField& momentumRes,
    const scalarField& pressureRes,
    const pointVectorField& pointD,
    const pointScalarField& pointP
)
{
    Info<< "Calculating the Jacobian using finite differences" << endl;

    // Small number used for perturbations
    const scalar relEps = 1e-8; // sqrt of computer tolerance
    const scalar typicalDisplacementValue =
        readScalar(solidModelDict().lookup("typicalDisplacementValue"));
    const scalar typicalPressureValue =
        readScalar(solidModelDict().lookup("typicalPressureValue"));
    const scalar epsD = relEps*max
        (
            average(mag(pointD.primitiveField())),
            typicalDisplacementValue
        );
    const scalar epsP = relEps*max
        (
            average(mag(pointP.primitiveField())),
            typicalPressureValue
        );

    Info<< "epsD = " << epsD << ", epsP = " << epsP << endl;

    // Create fields to be used for perturbations
    vectorField momentumResPerturb = momentumRes;
    scalarField pressureResPerturb = pressureRes;
    pointVectorField pointDPerturb("pointDPerturb", pointD);
    pointScalarField pointPPerturb("pointPPerturb", pointP);

    // Displacement coefficients
    forAll(pointD, blockRowI)
    {
        forAll(pointD, blockColI)
        {
            // For each component of pointD, sequentially apply a perturbation
            // and then calculate the resulting residuals
            for (label cmptI = 0; cmptI < vector::nComponents; cmptI++)
            {
                // Reset pointDPerturb and multiply by 1.0 to avoid it being
                // removed from the object registry
                pointDPerturb = 1.0*pointD;

                // Perturb this component of pointD
                pointDPerturb[blockColI].component(cmptI) =
                    pointD[blockColI].component(cmptI) + epsD;

                // Calculate residualD with this component perturbed
                momentumResPerturb = residualD(pointDPerturb, pointP);

                // Calculate residualP with this component perturbed
                pressureResPerturb = residualP(pointDPerturb, pointP);

                // Calculate each component
                const vector tangCmptD
                (
                    (momentumResPerturb[blockRowI] - momentumRes[blockRowI])
                   /epsD
                );
                const scalar tangCmptP
                (
                    (pressureResPerturb[blockRowI] - pressureRes[blockRowI])
                   /epsD
                );

                // Insert components
                matrix(blockRowI, blockColI)(0,cmptI) =
                    tangCmptD.component(vector::X);
                matrix(blockRowI, blockColI)(1,cmptI) =
                    tangCmptD.component(vector::Y);
                matrix(blockRowI, blockColI)(2,cmptI) =
                    tangCmptD.component(vector::Z);
                matrix(blockRowI, blockColI)(3,cmptI) = tangCmptP;
            }
        }
    }

    // Pressure coefficients
    forAll(pointP, blockRowI)
    {
        forAll(pointP, blockColI)
        {
            // Reset pointPPerturb and multiply by 1.0 to avoid it being
            // removed from the object registry
            pointPPerturb = 1.0*pointP;

            // Perturb pointP
            pointPPerturb[blockColI] = pointP[blockColI] + epsP;

            // Calculate residualD with this component perturbed
            momentumResPerturb = residualD(pointD, pointPPerturb);

            // Calculate residualP with this component perturbed
            pressureResPerturb = residualP(pointD, pointPPerturb);

            // Calculate the components
            const vector tangCmptD
            (
                (momentumResPerturb[blockRowI] - momentumRes[blockRowI])
               /epsP
            );
            const scalar tangCmptP
            (
                (pressureResPerturb[blockRowI] - pressureRes[blockRowI])
               /epsP
            );

            // Insert components
            matrix(blockRowI, blockColI)(0,3) = tangCmptD.component(vector::X);
            matrix(blockRowI, blockColI)(1,3) = tangCmptD.component(vector::Y);
            matrix(blockRowI, blockColI)(2,3) = tangCmptD.component(vector::Z);
            matrix(blockRowI, blockColI)(3,3) = tangCmptP;
        }
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomPressureDisplacementSolid::vertexCentredLinGeomPressureDisplacementSolid
(
        Time& runTime,
        const word& region
)
:
    solidModel(typeName, runTime, region),
    dualMechanicalPtr_
    (
        new dualMechanicalModel
        (
            dualMesh(),
            nonLinGeom(),
            incremental(),
            mechanical(),
            dualMeshMap().dualFaceToCell()
        )
    ),
    steadyState_(false),
    compactStencil_(false),
    zeta_(solidModelDict().lookupOrDefault<scalar>("zeta", 0.2)),
    pressureSmoothingFactor_
    (
        solidModelDict().lookupOrDefault<scalar>("pressureSmoothFactor", 0.0)
    ),
    pointK_(mechanical().pBulkModulus()),
    twoD_(sparseMatrixExtendedTools::checkTwoD(mesh())),
    fixedDofs_(mesh().nPoints(), false),
    fixedDofValues_(fixedDofs_.size(), vector::zero),
    fixedDofDirections_(fixedDofs_.size(), symmTensor::zero),
    fixedDofScale_
    (
        solidModelDict().lookupOrDefault<scalar>
        (
            "fixedDofScale",
            (
                average(mechanical().impK())
               *Foam::sqrt(gAverage(mesh().magSf()))
            ).value()
        )
    ),
    pointP_
    (
        IOobject
        (
            "pointP",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimPressure, 0)
    ),
    pointU_
    (
        IOobject
        (
            "pointU",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity, vector::zero)
    ),
    pointA_
    (
        IOobject
        (
            "pointA",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedVector("0", dimVelocity/dimTime, vector::zero)
    ),
    pointRho_
    (
        IOobject
        (
            "point(rho)",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimDensity, 0.0)
    ),
    pointVol_
    (
        IOobject
        (
            "pointVolumes",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimVolume, 0.0)
    ),
    pointGradD_
    (
        IOobject
        (
            "pGrad(D)",
            runTime.timeName(),
            mesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        pMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualGradDf_
    (
        IOobject
        (
            "grad(D)f",
            runTime.timeName(),
            dualMesh(),
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedTensor("zero", dimless, tensor::zero),
        "calculated"
    ),
    dualSigmaf_
    (
        IOobject
        (
            "sigmaf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedSymmTensor("zero", dimPressure, symmTensor::zero),
        "calculated"
    ),
    dualPf_
    (
        IOobject
        (
            "pf",
            runTime.timeName(),
            dualMesh(),
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        dualMesh(),
        dimensionedScalar("zero", dimPressure, 0),
        "calculated"
    ),
    volP_
    (
        IOobject
        (
            "p",
            runTime.timeName(),
            mesh(),
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh(),
        dimensionedScalar("zero", dimPressure, 0),
        "calculated"
    ),
    globalPointIndices_(mesh())
{
    // Create dual mesh and set write option
    dualMesh().objectRegistry::writeOpt() = IOobject::NO_WRITE;

    // pointD field must be defined
    pointDisRequired();

    // Set fixed degree of freedom list
    setFixedDofs(pointD(), fixedDofs_, fixedDofValues_, fixedDofDirections_);

    // Set point density field
    mechanical().volToPoint().interpolate(rho(), pointRho_);

    // Set the pointVol field
    // Map dualMesh cell volumes to the primary mesh points
#ifdef OPENFOAM_NOT_EXTEND
    scalarField& pointVolI = pointVol_;
#else
    scalarField& pointVolI = pointVol_.internalField();
#endif
    const scalarField& dualCellVol = dualMesh().V();
    const labelList& dualCellToPoint = dualMeshMap().dualCellToPoint();
    forAll(dualCellToPoint, dualCellI)
    {
        // Find point which maps to this dual cell
        const label pointID = dualCellToPoint[dualCellI];

        // Map the cell volume
        pointVolI[pointID] = dualCellVol[dualCellI];
    }

    // Store old time fields
    pointD().oldTime().storeOldTime();
    pointP_.oldTime().storeOldTime();
    pointU_.oldTime().storeOldTime();
    pointA_.storeOldTime();

    // Write fixed degree of freedom equation scale
    Info<< "fixedDofScale: " << fixedDofScale_ << endl;

    // Disable the writing of the unused fields
    D().writeOpt() = IOobject::NO_WRITE;
    D().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    DD().writeOpt() = IOobject::NO_WRITE;
    DD().oldTime().oldTime().writeOpt() = IOobject::NO_WRITE;
    U().writeOpt() = IOobject::NO_WRITE;
    pointDD().writeOpt() = IOobject::NO_WRITE;
}


// * * * * * * * * * * * * * * * *  Destructors  * * * * * * * * * * * * * * //

vertexCentredLinGeomPressureDisplacementSolid::
~vertexCentredLinGeomPressureDisplacementSolid()
{
#ifdef USE_PETSC
    if (Switch(solidModelDict().lookup("usePETSc")))
    {
        PetscFinalize();
    }
#endif
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool vertexCentredLinGeomPressureDisplacementSolid::evolve()
{
    Info<< "Evolving solid solver" << endl;

    // Print out the pressure smoothing coefficient
    Info<< "    Laplacian equation will be solved for pressure" << nl
        << "    pressureSmoothingFactor: "
        << pressureSmoothingFactor_
        << endl;

    // Initialise matrix where each coefficient is a 4x4 tensor
    sparseMatrixExtended matrix(sum(globalPointIndices_.stencilSize()));

    // Initialise matrix calculated using finite differences
    sparseMatrixExtended finiteDifferenceMatrix
    (
        sum(globalPointIndices_.stencilSize())
    );

    // Interpolate the the hydrostatic pressure to the dual faces
    dualPf_ = vfvc::interpolate
    (
        pointP_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint()
    );

    // Calculate gradD at the primary mesh points
    pointGradD_ = vfvc::pGrad
    (
        pointD(),
        mesh()
    );

    // Interpolate pointP to the primary cells
    volP_ = vfvc::interpolate
    (
        pointP_,
        mesh()
    );

    // Store material tangent field for the dual mesh faces
    Field<RectangularMatrix<scalar>> materialTangent
    (
        dualMechanicalPtr_().materialTangentFaceField()
    );

    // Calculate stress field at the dual faces
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Calculate stress for the primary cells
    mechanical().correct(sigma());

    // Global point index lists
    const boolList& ownedByThisProc = globalPointIndices_.ownedByThisProc();
    const labelList& localToGlobalPointMap =
        globalPointIndices_.localToGlobalPointMap();

    // Unknown pressure and displacement correction
    Field<RectangularMatrix<scalar>> pointDPcorr
    (
        pointD().internalField().size(),
        RectangularMatrix<scalar>(4,1,0)
    );

    // Newton-Raphson loop over momentum and pressure equations
    int iCorr = 0;
    scalar initResidualD = 0.0;
    scalar initResidualP = 0.0;
#ifdef OPENFOAM_NOT_EXTEND
    SolverPerformance<vector> solverPerf;
#else
    BlockSolverPerformance<vector> solverPerf;
#endif
    do
    {
        // Update gradD at the dual faces
        dualGradDf_ = vfvc::fGrad
        (
            pointD(),
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            zeta_,
            debug
        );

        // Update the pressure at the dual faces
        dualPf_ = vfvc::interpolate
        (
            pointP_,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint()
        );

        // Update gradD at the primary mesh points
        pointGradD_ = vfvc::pGrad
        (
            pointD(),
            mesh()
        );

        // Update the primary cell pressure
        volP_ = vfvc::interpolate
        (
            pointP_,
            mesh()
        );

        // Calculate stress at the dual faces
        dualMechanicalPtr_().correct(dualSigmaf_);

        // Initialise the source
        Field<RectangularMatrix<scalar>> source
        (
            mesh().nPoints(),
            RectangularMatrix<scalar>(4, 1, 0.0)
        );

        // Calculate d(pBar)/d(gradD) for the matrix coefficients
        const pointTensorField pBarSensitivity(-pointK_*tensor(I));

        pointP_.correctBoundaryConditions();
        pointD().correctBoundaryConditions();

        // Calculate the residuals of the momentum and pressure equations
        const vectorField momentumRes
        (
            residualD
            (
                pointD(),
                pointP_
            )
        );
        const scalarField pressureRes
        (
            residualP
            (
                pointD(),
                pointP_
            )
        );

        // Assemble the source
        forAll(momentumRes, pointI)
        {
            source[pointI](0,0) = -momentumRes[pointI].component(vector::X);
            source[pointI](1,0) = -momentumRes[pointI].component(vector::Y);
            source[pointI](2,0) = -momentumRes[pointI].component(vector::Z);
            source[pointI](3,0) = -pressureRes[pointI];
        }

        // Assemble the matrix once per outer iteration
        matrix.clear();

        // Update material tangent
        materialTangent = dualMechanicalPtr_().materialTangentFaceField();

        // Add div(sigma) pressure and displacement coefficients
        vfvm::divSigma
        (
            matrix,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            materialTangent,
            zeta_,
            debug
        );

        // Add laplacian coefficients to the pressure equation
        vfvm::laplacian
        (
            matrix,
            compactStencil_,
            mesh(),
            dualMesh(),
            dualMeshMap().dualFaceToCell(),
            dualMeshMap().dualCellToPoint(),
            pressureSmoothingFactor_,
            debug
        );

        // Add displacement coefficients of pressure equation
        vfvm::addPbar
        (
            matrix,
            mesh(),
            dualMeshMap().dualCellToPoint(),
            pointVol_,
            pBarSensitivity,
            debug
        );

        // Add source terms of pressure equation to the matrix
        vfvm::Sp
        (
            matrix,
            pointVol_,
            debug
        );

        // Calculate the matrix using finite difference
        if (Switch(solidModelDict().lookup("finiteDifferenceJacobian")))
        {
            vertexCentredLinGeomPressureDisplacementSolid::finiteDiffMatrix
            (
                finiteDifferenceMatrix,
                momentumRes,
                pressureRes,
                pointD(),
                pointP_
            );

            // Enforce the fixedDOFs on the linear system for the displacement
            sparseMatrixExtendedTools::enforceFixedDisplacementDof
            (
                finiteDifferenceMatrix,
                source,
                fixedDofs_,
                fixedDofDirections_,
                fixedDofValues_,
                fixedDofScale_
            );
        }

        if (debug > 1)
        {
            Info<< "Matrix before enforcing DOFs: " << endl << endl;
            matrix.print();
        }

        // Enforce fixed DOF on the linear system for the displacement
        sparseMatrixExtendedTools::enforceFixedDisplacementDof
        (
            matrix,
            source,
            fixedDofs_,
            fixedDofDirections_,
            fixedDofValues_,
            fixedDofScale_
        );

        if (debug > 1)
        {
            Info<< "Matrix after enforcing DOFs: " << endl << endl;
            matrix.print();
        }

        // Solve the linear system
        if (debug)
        {
            Info<< "bool vertexCentredLinGeomPressureDisplacementSolid::evolve(): "
                << " solving linear system: start" << endl;
        }

        Info<< "    Solving" << endl;

        if (Switch(solidModelDict().lookup("usePETSc")))
        {
#ifdef USE_PETSC
            fileName optionsFile(solidModelDict().lookup("optionsFile"));

            // Solve for displacement and pressure correction
            solverPerf = sparseMatrixExtendedTools::solveLinearSystemPETSc
            (
                matrix,
                source,
                pointDPcorr,
                twoD_,
                optionsFile,
                mesh().points(),
                ownedByThisProc,
                localToGlobalPointMap,
                globalPointIndices_.stencilSizeOwned(),
                globalPointIndices_.stencilSizeNotOwned(),
                solidModelDict().lookupOrDefault<bool>("debugPETSc", false)
            );
#else
            FatalErrorIn("vertexCentredLinGeomPressureDisplacementSolid::evolve()")
                << "PETSc not available. Please set the PETSC_DIR "
                << "environment variable and re-compile solids4foam"
                << abort(FatalError);
#endif
        }
        else
        {
            // Use Eigen SparseLU direct solver
            sparseMatrixExtendedTools::solveLinearSystemEigen
            (
                matrix, source, pointDPcorr, twoD_, true, debug
            );
        }

        if (debug)
        {
            Info<< "bool vertexCentredLinGeomPressureDisplacementSolid::evolve(): "
                << " solving linear system: end" << endl;
        }

        // Update point displacement field
        if (Switch(solidModelDict().lookup("lineSearch")))
        {
            notImplemented("Line search not implemented.")
        }
#ifdef OPENFOAM_NOT_EXTEND
        else if (mesh().relaxField(pointD().name()))
#else
        else if (mesh().solutionDict().relaxField(pointD().name()))
#endif
        {
            notImplemented("pointD or pointP relaxation not implemented.")
        }
        else
        {
            forAll(pointDPcorr, pointI)
            {
#ifdef OPENFOAM_NOT_EXTEND
                pointD().primitiveFieldRef()[pointI].component(0) +=
                    pointDPcorr[pointI](0,0);
                pointD().primitiveFieldRef()[pointI].component(1) +=
                    pointDPcorr[pointI](1,0);
                pointD().primitiveFieldRef()[pointI].component(2) +=
                    pointDPcorr[pointI](2,0);
                pointP_.primitiveFieldRef()[pointI] +=
                    pointDPcorr[pointI](3,0);
#else
                pointD().internalField()[pointI].component(0) +=
                    pointDPcorr[pointI](0,0);
                pointD().internalField()[pointI].component(1) +=
                    pointDPcorr[pointI](1,0);
                pointD().internalField()[pointI].component(2) +=
                    pointDPcorr[pointI](2,0);
                pointP_.internalField()[pointI] +=
                    pointDPcorr[pointI](3,0);
#endif
            }
        }
        pointD().correctBoundaryConditions();
        pointP_.correctBoundaryConditions();

        // Update point accelerations
        // Note: for NewmarkBeta, this needs to come before the pointU
        // update
#ifdef OPENFOAM_NOT_EXTEND
        pointA_.primitiveFieldRef() =
        vfvc::ddt
        (
            mesh().ddtScheme("ddt(pointU)"),
            mesh().d2dt2Scheme("d2dt2(pointD)"),
            pointU_
        );

        // Update point velocities
        pointU_.primitiveFieldRef() =
        vfvc::ddt
        (
            mesh().ddtScheme("ddt(pointD)"),
            mesh().d2dt2Scheme("d2dt2(pointD)"),
            pointD()
        );
#else
        pointA_.internalField() =
        vfvc::ddt
        (
            mesh().schemesDict().ddtScheme("ddt(pointU)"),
            mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
            pointU_
        );

        // Update point velocities
        pointU_.internalField() =
        vfvc::ddt
        (
            mesh().schemesDict().ddtScheme("ddt(pointD)"),
            mesh().schemesDict().d2dt2Scheme("d2dt2(pointD)"),
            pointD()
        );
#endif
    }
    while
    (
        !converged
        (
            iCorr,
            initResidualD,
            initResidualP,
            cmptMax(solverPerf.nIterations()),
            pointD(),
            pointP_,
            pointDPcorr
        ) && ++iCorr
    );

    // Update gradD at the dual faces
    dualGradDf_ = vfvc::fGrad
    (
        pointD(),
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint(),
        zeta_,
        debug
    );

    // Update the hydrostatic pressure at the dual faces
    dualPf_ = vfvc::interpolate
    (
        pointP_,
        mesh(),
        dualMesh(),
        dualMeshMap().dualFaceToCell(),
        dualMeshMap().dualCellToPoint()
    );

    // Update gradD at the primary mesh points
    pointGradD_ = vfvc::pGrad
    (
        pointD(),
        mesh()
    );

    // Update primary cell pressure
    volP_ = vfvc::interpolate
    (
        pointP_,
        mesh()
    );

    // Update the increment of displacement
    pointDD() = pointD() - pointD().oldTime();

    // Calculate cell gradient
    // This assumes a constant gradient within each primary mesh cell
    // This is a first-order approximation
    gradD() = vfvc::grad(pointD(), mesh());

    // Map primary cell gradD field to sub-meshes for multi-material cases
    if (mechanical().PtrList<mechanicalLaw>::size() > 1)
    {
        mechanical().mapGradToSubMeshes(gradD());
    }

    // Update dual face stress field
    dualMechanicalPtr_().correct(dualSigmaf_);

    // Update primary mesh cell stress field, assuming it is constant per
    // primary mesh cell
    // This stress will be first-order accurate

    mechanical().correct(sigma());

    return true;
}

void vertexCentredLinGeomPressureDisplacementSolid::setTraction
(
    const label interfaceI,
    const label patchID,
    const vectorField& faceZoneTraction
)
{
    // Get point field on patch
    const vectorField traction
    (
        globalPatches()[interfaceI].globalPointToPatch
        (
            globalPatches()[interfaceI].interpolator().faceToPointInterpolate
            (
                faceZoneTraction
            )()
        )
    );

    // Lookup point patch field
#ifdef OPENFOAM_NOT_EXTEND
    pointPatchVectorField& ptPatch = pointD().boundaryFieldRef()[patchID];
#else
    pointPatchVectorField& ptPatch = pointD().boundaryField()[patchID];
#endif

    if (isA<solidTractionPointPatchVectorField>(ptPatch))
    {
        solidTractionPointPatchVectorField& patchD =
            refCast<solidTractionPointPatchVectorField>(ptPatch);

        patchD.traction() = traction;
    }
    else
    {
        FatalErrorIn
        (
            "void Foam::vertexCentredLinGeomPressureDisplacementSolid::setTraction\n"
            "(\n"
            "    fvPatchVectorField& tractionPatch,\n"
            "    const vectorField& traction\n"
            ")"
        )   << "Boundary condition "
            << ptPatch.type()
            << " for point patch " << ptPatch.patch().name()
            << " should instead be type "
            << solidTractionPointPatchVectorField::typeName
            << abort(FatalError);
    }
}

void vertexCentredLinGeomPressureDisplacementSolid::writeFields
(
    const Time& runTime
)
{
    // Calculate gradD at the primary points using least squares: this
    // should be second-order accurate (... I think).
    const pointTensorField pGradD(vfvc::pGrad(pointD(), mesh()));

    // Calculate strain at the primary points based on pGradD
    // Note: the symm operator is not defined for pointTensorFields so we
    // will do it manually
    pointSymmTensorField pEpsilon
    (
        IOobject
        (
            "pEpsilon",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedSymmTensor("0", dimless, symmTensor::zero)
    );

#ifdef FOAMEXTEND
    pEpsilon.internalField() = symm(pGradD.internalField());
#else
    pEpsilon.primitiveFieldRef() = symm(pGradD.internalField());
#endif
    pEpsilon.write();

    // Equivalent strain at the points
    pointScalarField pEpsilonEq
    (
        IOobject
        (
            "pEpsilonEq",
            runTime.timeName(),
            runTime,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        pMesh(),
        dimensionedScalar("0", dimless, 0.0)
    );

#ifdef FOAMEXTEND
    pEpsilonEq.internalField() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#else
    pEpsilonEq.primitiveFieldRef() =
        sqrt((2.0/3.0)*magSqr(dev(pEpsilon.internalField())));
#endif
    pEpsilonEq.write();

    Info<< "Max pEpsilonEq = " << gMax(pEpsilonEq) << endl;

    solidModel::writeFields(runTime);
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace solidModels

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
