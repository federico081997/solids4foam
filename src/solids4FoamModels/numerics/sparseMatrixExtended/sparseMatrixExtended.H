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
    sparseMatrixExtended

Description
    A simple but convenient matrix for storing sparse matrices.
    Values are stored using a HashTable where the Hash is based on the
    (rowI, columns) key.

    It is assumed that every value is a 4x4 RectangularMatrix

    Used in coupled solid mechanics solvers in solids4foam.

Author
    Philip Cardiff, UCD.
    Federico Mazzanti, UCD.

SourceFiles
    sparseMatrixExtended.C

\*---------------------------------------------------------------------------*/

#ifndef sparseMatrixExtended_H
#define sparseMatrixExtended_H

#include "HashTable.H"
#include "FixedList.H"
#include "scalarMatrices.H"
#include "refCount.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    typedef HashTable
    <
        scalarRectangularMatrix,
        FixedList<label, 2>,
        FixedList<label, 2>::Hash<>
    > sparseMatrixExtendedData;

/*---------------------------------------------------------------------------*\
                     Class sparseMatrixExtended Declaration
\*---------------------------------------------------------------------------*/

class sparseMatrixExtended
:
    public refCount
{
    // Private data

        //- Hashtable (map) where the values are stored
        sparseMatrixExtendedData data_;


    // Private Member Functions

        //- Disallow default bitwise assignment
        void operator=(const sparseMatrixExtended&);

public:

    //- Runtime type information
    TypeName("sparseMatrixExtended");


    // Constructors

        //- Construct given approximate size
        sparseMatrixExtended(const label size = 128);

        //- Copy construct
        sparseMatrixExtended(const sparseMatrixExtended& sm);

        //- Construct from tmp\<sparseMatrixExtended\> deleting argument
        sparseMatrixExtended(const tmp<sparseMatrixExtended>& tsm);

    // Destructor

        virtual ~sparseMatrixExtended()
        {}


    // Member Functions

        // Access

            //- Const access to the data
            const sparseMatrixExtendedData& data() const
            {
                return data_;
            }

            //- Non-const access to the data
            sparseMatrixExtendedData& data()
            {
                return data_;
            }

            //- Count and return the number of block rows in the matrix
            label nBlockRows() const;

            //- Print out the matrix coefficients
            void print() const;

        // Modifiers

            //- Clear the matrix but do not change its capacity
            void clear()
            {
                data_.clear();
            }

        // Operators

            //- Non-const access to a value entry. If the entry does not exist
            //  then create it and initialise it to zero
            scalarRectangularMatrix& operator()
            (
                const label rowI,
                const label colI
            );

            //- Operator to add a sparseMatrixExtended A to this object
            void operator+=(const sparseMatrixExtended& A);

            //- Operator to subtract a sparseMatrixExtended A from this object
            void operator-=(const sparseMatrixExtended& A);

            //- Operator to add a sparseMatrixExtended A to this object and
            //  return the result
            tmp<sparseMatrixExtended> operator+
            (
                const sparseMatrixExtended& A
            ) const;

            //- Operator to multiply the diagonal coefficients of this
            //  sparseMatrixExtended by a scalar and return the result
            tmp<sparseMatrixExtended> operator*(const scalar s) const;

            //- Operator to multiply the diagonal coefficients of this
            //  sparseMatrixExtended by a scalarField and return the result
            //  It is assumed that the field index corresponds to the matrix row
            tmp<sparseMatrixExtended> operator*(const scalarField& sf) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
