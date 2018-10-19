/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
 Copyright (C) 2015-2018
 
 Author: Antonin Portelli <antonin.portelli@me.com>
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
 
 See the full license in the file "LICENSE" in the top level distribution
 directory.
 *******************************************************************************/

#include <Hadrons/Application.hpp>
#include <Hadrons/Modules.hpp>

using namespace Grid;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    // initialization //////////////////////////////////////////////////////////
    Grid_init(&argc, &argv);
    HadronsLogError.Active(GridLogError.isActive());
    HadronsLogWarning.Active(GridLogWarning.isActive());
    HadronsLogMessage.Active(GridLogMessage.isActive());
    HadronsLogIterative.Active(GridLogIterative.isActive());
    HadronsLogDebug.Active(GridLogDebug.isActive());
    LOG(Message) << "Grid initialized" << std::endl;


    // run setup ///////////////////////////////////////////////////////////////
    Application              application;
    std::vector<std::string> flavour = {"h"}; //{"l", "s", "c1", "c2", "c3"};
    std::vector<double>      mass    = {.2}; //{.01, .04, .2  , .25 , .3  };

    //lepton parameters
    std::vector<std::string> lepton_flavour    = {"mu"};
    std::vector<double>      lepton_mass    = {.2};
    std::vector<std::string> lepton_twist = {"0.0 0.0 0.0 0.5"};

    std::vector<double>      lepton_energy (lepton_flavour.size(),.0);
    std::vector<Real> tw;
    for (unsigned int i = 0; i < lepton_mass.size(); ++i)
    {
        tw  = strToVec<Real>(lepton_twist[i]);
        lepton_energy[i] = lepton_mass[i]*lepton_mass[i];
	for(unsigned int mu = 0; mu < 3; mu++){
	    lepton_energy[i] += 2*M_PI*tw[mu]*2*M_PI*tw[mu]/(GridDefaultLatt()[mu]*GridDefaultLatt()[mu]);
	}
	lepton_energy[i] = std::sqrt(lepton_energy[i]);
    }


    unsigned int  nt    = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;
    globalPar.trajCounter.end   = 1520;
    globalPar.trajCounter.step  = 20;
    globalPar.runId             = "test";
    application.setPar(globalPar);
    // gauge field
    application.createModule<MGauge::Unit>("gauge");
    // unit gauge field for lepton 
    application.createModule<MGauge::Unit>("free_gauge");
    // pt source
    MSource::Point::Par ptPar;
    ptPar.position = "0 0 0 0";
    application.createModule<MSource::Point>("pt", ptPar);
    // sink
    MSink::Point::Par sinkPar;
    sinkPar.mom = "0 0 0";
    application.createModule<MSink::ScalarPoint>("sink", sinkPar);
    
    // set fermion boundary conditions to be periodic space, antiperiodic time.
    std::string boundary = "1 1 1 -1";

    //stochastic photon field
    MGauge::StochEm::Par photonPar;
    photonPar.gauge = PhotonR::Gauge::feynman;
    photonPar.zmScheme = PhotonR::ZmScheme::qedL;
    application.createModule<MGauge::StochEm>("ph_field", photonPar);


    for (unsigned int i = 0; i < lepton_mass.size(); ++i)
    {

        // Wall source for lepton
        MSource::Wall::Par wallPar;
        wallPar.tA       = 0;
        wallPar.tB       = nt-1;
        wallPar.mom       = "0. 0. 0. 0.";
        wallPar.energy = lepton_energy[i];
        application.createModule<MSource::Wall>("wall", wallPar);


        // lepton action
        MAction::DWF::Par actionPar_lep;
        actionPar_lep.gauge = "free_gauge";
        actionPar_lep.Ls    = 8;
        actionPar_lep.M5    = 1.8;
        actionPar_lep.mass  = lepton_mass[i];
        actionPar_lep.boundary = "1 1 1 1";//boundary conditions via twist in freeprop
        application.createModule<MAction::DWF>("free_DWF_" + lepton_flavour[i], actionPar_lep);


	//Aslash insertion at lepton
        MSource::SeqAslash::Par seqKl2Par;
        seqKl2Par.q         = "wall";
        seqKl2Par.tA        = 0;
        seqKl2Par.tB        = nt-1;
        seqKl2Par.photon    = "ph_field";
        seqKl2Par.mom       = "0. 0. 0. 0.";
        application.createModule<MSource::SeqAslash>("lVA" + lepton_flavour[i], seqKl2Par);



    	// free propagators for lepton
    	MFermion::FreeProp::Par freeKl2Par;
	freeKl2Par.source = "lVA" + lepton_flavour[i];
    	freeKl2Par.action = "free_DWF_" + lepton_flavour[i]; 
	freeKl2Par.twist = lepton_twist[i];
    	freeKl2Par.mass = lepton_mass[i];
    	application.createModule<MFermion::FreeProp>("lVA_" + lepton_flavour[i],
							 freeKl2Par);
    }


    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // quark actions
        MAction::DWF::Par actionPar;
        actionPar.gauge = "gauge";
        actionPar.Ls    = 8;
        actionPar.M5    = 1.8;
        actionPar.mass  = mass[i];
        actionPar.boundary = boundary;
        application.createModule<MAction::DWF>("DWF_" + flavour[i], actionPar);

        
        // solvers
        MSolver::RBPrecCG::Par solverPar;
        solverPar.action       = "DWF_" + flavour[i];
        solverPar.residual     = 1.0e-8;
        solverPar.maxIteration = 10000;
        application.createModule<MSolver::RBPrecCG>("CG_" + flavour[i],
                                                    solverPar);
        
        // propagators
        MFermion::GaugeProp::Par quarkPar;
        quarkPar.solver = "CG_" + flavour[i];
        quarkPar.source = "pt";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i],
							 quarkPar);

	//seq sources with Aslash insertion
        MSource::SeqAslash::Par seqPar;
        seqPar.q         = "Qpt_" + flavour[i];
        seqPar.tA        = 0;
        seqPar.tB        = nt-1;
        seqPar.photon    = "ph_field";
        seqPar.mom       = "0. 0. 0. 0.";
        application.createModule<MSource::SeqAslash>("Qpt_" + flavour[i] 
						    + "_seq_lV_ph", seqPar);
        // seq propagator with Aslash insertion
        MFermion::GaugeProp::Par quarkPar_seq;
        quarkPar_seq.solver = "CG_" + flavour[i];
        quarkPar_seq.source = "Qpt_" + flavour[i] + "_seq_lV_ph";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i] 
						+ "_seq_lV_ph_" + flavour[i], 
							quarkPar_seq);


    }

 
    //Kl2 contraction
    MContraction::WeakMesonDecayKl2::Par kl2Par;
    kl2Par.q1      = "Qpt_" + flavour[0] + "_seq_lV_ph_" + flavour[0];
    kl2Par.q2      = "Qpt_" + flavour[0];
    kl2Par.lepton      = "lVA_" + lepton_flavour[0];
    kl2Par.output  = "QED_mat/l_weakdecay_" + flavour[0] + flavour[0] + "_to_"  + lepton_flavour[0];
    application.createModule<MContraction::WeakMesonDecayKl2>("final",
                                                      kl2Par);

   
    // execution
    application.saveParameterFile("QED_kl2.xml");
    application.run();
    
    // epilogue
    LOG(Message) << "Grid is finalizing now" << std::endl;
    Grid_finalize();
    
    return EXIT_SUCCESS;
}
