/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules/MSource/Wall5d.hpp

Copyright (C) 2015-2018

Author: Antonin Portelli <antonin.portelli@me.com>
Author: Lanny91 <andrew.lawson@gmail.com>

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

#ifndef Hadrons_MSource_Wall5dSource_hpp_
#define Hadrons_MSource_Wall5dSource_hpp_

#include <Grid/Hadrons/Global.hpp>
#include <Grid/Hadrons/Module.hpp>
#include <Grid/Hadrons/ModuleFactory.hpp>

#include <Grid/Hadrons/Modules/MFermion/GaugeProp.hpp>

BEGIN_HADRONS_NAMESPACE

/*
 
 Wall5d source
 -----------------------------
 * src_x = delta(x_3 - tW) * exp(i x.mom)
 
 * options:
 - tW: source timeslice (integer)
 - mom: momentum insertion, space-separated float sequence (e.g ".1 .2 1. 0.")
 
 */

/******************************************************************************
 *                         Wall5d                                               *
 ******************************************************************************/
BEGIN_MODULE_NAMESPACE(MSource)

class Wall5dPar: Serializable
{
public:
    GRID_SERIALIZABLE_CLASS_MEMBERS(Wall5dPar,
                                    unsigned int, tW,
                                    std::string, mom,
                                    bool, fiveD,
				    int, Ls);
};

template <typename FImpl>
class TWall5d: public Module<Wall5dPar>
{
public:
    FERM_TYPE_ALIASES(FImpl,);
public:
    // constructor
    TWall5d(const std::string name);
    // destructor
    virtual ~TWall5d(void) {};
    // dependency relation
    virtual std::vector<std::string> getInput(void);
    virtual std::vector<std::string> getOutput(void);
protected:
    // setup
    virtual void setup(void);
    // execution
    virtual void execute(void);
private:
    unsigned int Ls_;
//private:
//    bool        hasPhase_{false};
//    std::string momphName_, tName_;
};

MODULE_REGISTER_TMP(Wall5d, TWall5d<FIMPL>, MSource);

/******************************************************************************
 *                 TWall5d implementation                                       *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename FImpl>
TWall5d<FImpl>::TWall5d(const std::string name)
: Module<Wall5dPar>(name)
//, momphName_ (name + "_momph")
//, tName_ (name + "_t")
{}

// dependencies/products ///////////////////////////////////////////////////////
template <typename FImpl>
std::vector<std::string> TWall5d<FImpl>::getInput(void)
{
    std::vector<std::string> in = {};
    
    return in;
}

template <typename FImpl>
std::vector<std::string> TWall5d<FImpl>::getOutput(void)
{
    std::vector<std::string> out = {getName()};
    
    return out;
}

// setup ///////////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall5d<FImpl>::setup(void)
{
if (par().fiveD) Ls_ = par().Ls; //env().getObjectLs(par().solver);
    else Ls_ = 1; //later have only one wall src module
    envCreateLat(PropagatorField, getName(), Ls_);
    envTmpLat(Lattice<iScalar<vInteger>>, "t");
    envTmpLat(LatticeComplex, "ph");
    envTmpLat(LatticeComplex, "coor");
    envTmpLat(PropagatorField, "fourdsrc");
    envTmpLat(FermionField, "tmpferm");
    envTmpLat(FermionField, "tmpfermfive", Ls_);
    envTmpLat(PropagatorField, "srcbuf");

}

// execution ///////////////////////////////////////////////////////////////////
template <typename FImpl>
void TWall5d<FImpl>::execute(void)
{    
    LOG(Message) << "Generating Wall5d source at t = " << par().tW 
                 << " with momentum " << par().mom << std::endl;
    
    auto  &src = envGet(PropagatorField, getName());
//    auto  &ph  = envGet(LatticeComplex, momphName_);
    envGetTmp(Lattice<iScalar<vInteger>>, t);
    envGetTmp(LatticeComplex, ph);
    envGetTmp(PropagatorField, fourdsrc);
    envGetTmp(FermionField, tmpferm);
    envGetTmp(FermionField, tmpfermfive);
    envGetTmp(PropagatorField, srcbuf);
    //auto  &t   = envGet(Lattice<iScalar<vInteger>>, tName_);

    
    //if (!hasPhase_)
    //{
        Complex           i(0.0,1.0);
        std::vector<Real> p;

        envGetTmp(LatticeComplex, coor);
        p  = strToVec<Real>(par().mom);
        ph = zero;
        for(unsigned int mu = 0; mu < env().getNd(); mu++)
        {
            LatticeCoordinate(coor, mu);
            ph = ph + (p[mu]/env().getGrid()->_fdimensions[mu])*coor;
        }
        ph = exp((Real)(2*M_PI)*i*ph);
        LatticeCoordinate(t, Tp);
        //hasPhase_ = true;
    //}



    fourdsrc = 1.;
    fourdsrc = where((t == par().tW), fourdsrc*ph, 0.*fourdsrc);



    for(unsigned int s=0; s< Ls_; ++s){
	InsertSlice(fourdsrc, src, s, 0);
    }


}

END_MODULE_NAMESPACE

END_HADRONS_NAMESPACE

#endif // Hadrons_MSource_Wall5dSource_hpp_
