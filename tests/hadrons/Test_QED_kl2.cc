/*******************************************************************************
 Grid physics library, www.github.com/paboyle/Grid
 
 Source file: tests/hadrons/Test_hadrons_spectrum.cc
 
 Copyright (C) 2015-2018
 
 Author: Antonin Portelli <antonin.portelli@me.com>
 Author: Vera Guelpers    <v.m.guelpers@soton.ac.uk>
 
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

#include <Grid/Hadrons/Application.hpp>
#include <Grid/Hadrons/Modules.hpp>

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
    std::vector<std::string> lepton_flavour    = {"mu"};
    std::vector<double>      lepton_mass    = {.2};

    unsigned int  nt    = GridDefaultLatt()[Tp];
    
    // global parameters
    Application::GlobalPar globalPar;
    globalPar.trajCounter.start = 1500;//5061;//5056;
    globalPar.trajCounter.end   = 1501;//5162;//5061;
    globalPar.trajCounter.step  = 1;
    globalPar.seed              = "1 2 3 4";//"1 2 3 10";//"1 2 3 9";
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
    std::string boundary = "1 1 1 -1";//"1 1 1 -1";

    //stochastic photon field
    MGauge::StochEm::Par photonPar;
    photonPar.gauge = PhotonR::Gauge::feynman;
    photonPar.zmScheme = PhotonR::ZmScheme::qedL;
    application.createModule<MGauge::StochEm>("ph_field", photonPar);





    // Wall source for lepton
    MSource::Wall5d::Par wallPar;
    wallPar.tW = 4;//26;
    wallPar.mom       = "0. 0. 0. 0.";
    wallPar.fiveD = true;
    wallPar.Ls = 8;
    application.createModule<MSource::Wall5d>("wall", wallPar);




    for (unsigned int i = 0; i < lepton_mass.size(); ++i)
    {
        // actions
        MAction::DWF::Par actionPar_lep;
        actionPar_lep.gauge = "free_gauge";
        actionPar_lep.Ls    = 8;
        actionPar_lep.M5    = 1.8;
        actionPar_lep.mass  = lepton_mass[i];
        actionPar_lep.boundary = "1 1 1 1";//boundary conditions via twist in freeprop
        application.createModule<MAction::DWF>("free_DWF_" + lepton_flavour[i], actionPar_lep);


    	MSource::SeqConserved::Par seqKl2Par;
    	seqKl2Par.q         = "wall";
    	seqKl2Par.action    = "free_DWF_" + lepton_flavour[i]; 
    	seqKl2Par.tA        = 0;
    	seqKl2Par.tB        = nt-1;
    	seqKl2Par.curr_type = Current::Vector;
    	seqKl2Par.mu_min	   = 0;
    	seqKl2Par.mu_max	   = 3;
    	seqKl2Par.mom       = "0. 0. 0. 0.";
    	seqKl2Par.photon	   = "ph_field";
    	application.createModule<MSource::SeqConserved>("VA",seqKl2Par);



    	// free propagators
    	MFermion::FreeProp::Par freeKl2Par;
    	freeKl2Par.source = "VA";
    	freeKl2Par.action = "free_DWF_" + lepton_flavour[i]; 
    	freeKl2Par.twist = "0 0 0 0.5";
    	freeKl2Par.mass = lepton_mass[i];
    	application.createModule<MFermion::FreeProp>("VA_" + lepton_flavour[i],
							 freeKl2Par);

    }




    for (unsigned int i = 0; i < flavour.size(); ++i)
    {
        // actions
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




	//seq sources with conserved vector and photon insertion
        MSource::SeqConserved::Par seqPar_V;
        seqPar_V.q         = "Qpt_" + flavour[i] + "_5d";
        seqPar_V.action    = "DWF_" + flavour[i];
        seqPar_V.tA        = 0;
        seqPar_V.tB        = nt-1;
        seqPar_V.curr_type = Current::Vector;
	seqPar_V.mu_min	   = 0;
	seqPar_V.mu_max	   = 3;
        seqPar_V.mom       = "0. 0. 0. 0.";
	seqPar_V.photon	   = "ph_field";
        application.createModule<MSource::SeqConserved>("Qpt_" + flavour[i] 
						    + "_seq_V_ph", seqPar_V);
        //seq propagator with conserved vector and photon insertion
        MFermion::GaugeProp::Par quarkPar_seq_V;
        quarkPar_seq_V.solver = "CG_" + flavour[i];
        quarkPar_seq_V.source = "Qpt_" + flavour[i] + "_seq_V_ph";
        application.createModule<MFermion::GaugeProp>("Qpt_" + flavour[i] 
						+ "_seq_V_ph_" + flavour[i], 
							quarkPar_seq_V);


    }




//LATER: LOOP OVER QUARK/LEPTON FLAVOURS


    //Kl2 hadron contraction
    MContraction::WeakMesonDecayKl2::Par kl2Par;
    kl2Par.q1      = "Qpt_" + flavour[0] + "_seq_V_ph_" + flavour[0];
    kl2Par.q2      = "Qpt_" + flavour[0];
    kl2Par.lepton      = "VA_" + lepton_flavour[0];
    //kl2Par.gammas  = "(Gamma5 Gamma5)";//don't need gammas, at least not two of them
    kl2Par.output  = "QED_mat_test/weakdecay_" + flavour[0] + flavour[0] + "_to_"  + lepton_flavour[0];
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
