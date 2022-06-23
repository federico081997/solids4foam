/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Class
    Foam::ListSumOp

Description
    Op function class to perform operations on a list.

Author
    Philip Cardiff, UCD.  All rights reserved

\*---------------------------------------------------------------------------*/

#ifndef ListSumOp_H
#define ListSumOp_H

#include "Field.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    template<class Type>
    struct ListSumOp
    {
        public:

        // Op function
        List<Type> operator()
        (
            const List<Type>& x, const List<Type>& y
        ) const
        {
            if (x.size() != y.size())
            {
                FatalErrorIn
                (
                    "List<Type> ListSumOp()"
                    "(const List<Type>& x, const List<Type>& y) const"
                )   << "The list sizes do not match! "
                    << "x: " << x.size() << ", y: " << y.size()
                    << abort(FatalError);
            }

            List<Type> result(x.size());

            forAll(result, i)
            {
                result[i] = x[i] + y[i];
            }

            return result;
        }
    };

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //