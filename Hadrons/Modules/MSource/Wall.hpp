/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: Hadrons/Modules/MSource/Wall.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>
Author: Vera Guelpers    <Vera.Guelpers@ed.ac.uk>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */

#ifndef Hadrons_MSource_WallSource_hpp_
#define Hadrons_MSource_WallSource_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Module.hpp>
#include <Hadrons/ModuleFactory.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Wall source
 -----------------------------
 * src_x = theta(x_3 - tA) * theta(tB - x_3) * exp(i x.mom) * exp(energy*t)
 
 * options:
 - tA: begin timeslice (integer)
 - tB: end timesilce (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 - energy: for exp(energy*t), (double)
 
 */

/******************************************************************************
 *                         Wall                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class WallPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(WallPar,
                                    unsigned int,   tA,
                                    unsigned int,   tB,
                                    std::string, mom,
				    double, energy);
};

template <typename FImpl>
class TWall: public Module<WallPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWall(const std::string name);
    // destructor
    virtual ~TWall(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    bool        hasPhase_{false};
    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(Wall, TWall<FIMPL>, MSource);

/******************************************************************************
 *                 TWall implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWall<FImpl>::TWall(const std::string name)
: Module<WallPar>(name)
, momphName_ (name + "_momph")
, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWall<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::setup(void)
{
    envCreateLat(PropagatorField, getName());
    envCache(Lattice<iScalar<vInteger>>, tName_, 1, envGetGrid(LatticeComplex));
    envCacheLat(LatticeComplex, momphName_);
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(LatticeComplex, "tcoor");
}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall<FImpl>::execute(void)
{    
    if (par().tA == par().tB)
    {
	LOG(Message) << "Generating wall source at t = " << par().tA
                 << " with momentum " << par().mom << std::endl;
    }
    else
    {
	LOG(Message) << "Generating wall source for " << par().tA << " <= t <= " << par().tB
                 << " with momentum " << par().mom << std::endl;
    }
    
    auto  &src = envGet(PropagatorField, getName());
    auto  &ph  = envGet(LatticeComplex, momphName_);
    auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);
    
    if (!hasPhase_)
    {
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = zero;
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getDim(mu))*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        hasPhase_ = true;
    }

    envGetTmp(LatticeComplex, tcoor);
    LatticeCoordinate(tcoor, Tp);

    src = 1.;
    src = where((t >= par().tA) and (t <= par().tB), src*ph*exp(tcoor*par().energy), 0.*src);
}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_WallSource_hpp_
