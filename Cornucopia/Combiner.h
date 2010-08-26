/*--
    Combiner.h  

    This file is part of the Cornucopia curve sketching library.
    Copyright (C) 2010 Ilya Baran (baran37@gmail.com)

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef COMBINER_H_INCLUDED
#define COMBINER_H_INCLUDED

#include "defs.h"
#include "Algorithm.h"
#include "VectorC.h"
#include "CurvePrimitive.h"

NAMESPACE_Cornu

SMART_FORW_DECL(PrimitiveSequence);

template<>
struct AlgorithmOutput<COMBINING> : public AlgorithmOutputBase
{
    PrimitiveSequenceConstPtr output;
};

template<>
class Algorithm<COMBINING> : public AlgorithmBaseTemplate<COMBINING>
{
public:
    //override
    std::string stageName() const { return "Combining"; }

private:
    friend class AlgorithmBase;
    static void _initialize();
};

END_NAMESPACE_Cornu

#endif //COMBINER_H_INCLUDED
