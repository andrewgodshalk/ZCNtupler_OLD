/*-----------------------------------------------------------------------------
ZCNtupler.cc
  Created:  2014-03-05 by Andrew Godshalk (godshalk@buffalo.edu)
  Modified: 2014-03-05 by Andrew Godshalk

Step-2 Ntupler, working from PATTuples created by the VHbb group.

Takes a standard python input file for configuration information
    - Default: "zc_ntupler.py"

Modifications:
2014-03-05 - Used to make plots from Step1 PATTuples. No selection as of yet.

-----------------------------------------------------------------------------*/


#include <TFile.h>
#include <TH1F.h>
#include <TROOT.h>
#include <TString.h>
#include <TSystem.h>
#include <TTree.h>

#include <iostream>
#include <vector>

#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"

#include "ZCHelperFunctions.h"

using namespace std;


int main( int argc, char* argv[] )
{
    // Initialize ROOT and FWLite Libraries
    gROOT->Reset();
    gSystem->Load("libFWCoreFWLite");
    gSystem->Load("libDataFormatsFWLite");
    AutoLibraryLoader::enable();

    // ---- Declare Working Variables ----
    // Timestamps
    string startTimeForDisplay = timestamp();
    string startTimeForOutput  = file_timestamp();

    // File info
    TString outputFilename = TString("zcNtuple_") + file_timestamp() + TString("_");
    TFile *outputFile, *inputFile;
    TTree *zjCandTree;

    // Input data information
    int maxEvents, skipEvents,
        runMin,    runMax;
    vector<edm::LuminosityBlockRange> jsonVector;    // Vector containing luminosity range to process.
    vector<string> inputFilenames;

    // Variables to map to tree
    int nJets_, nMuons_, nElecs_;
    EventInfo event_;
    JetInfo jetInfo_;
    LeptonInfo muonInfo_, elecInfo_;
    
    // Counts
    int nEventsLoopedOver = 0;      // For keeping track of how many events program has looped over, for skipping and maxEvent counts.
    TH1F *nEventsProcessed;

    // ---- Do some high-tech computering!! ----
    cout << "\n  ====== ZCNtupler (" << startTimeForDisplay << ") ======" << endl;
    
    // Read in parameter sets from python file
    if ( argc < 2 ) { cout << "\n  -->ERROR: Please input a python configuration file.\n  ====== ZCNtupler (ERROR) ======\n" << endl; return 0; }
    PythonProcessDesc builder(argv[1]);
    const edm::ParameterSet& inputInfo        = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("inputInformation" );
    const edm::ParameterSet& outputInfo       = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("outputInformation");
    const edm::ParameterSet& analysisCriteria = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("analysisCriteria" );

    // Read luminosity information if it exists
    if ( inputInfo.exists("lumisToProcess") )
    {
        vector<edm::LuminosityBlockRange> const & lumisTemp = inputInfo.getUntrackedParameter<vector<edm::LuminosityBlockRange> > ("lumisToProcess");
        jsonVector.resize( lumisTemp.size() );
        copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
    }

    // Read other input parameters
    maxEvents  = inputInfo.getParameter<int>("maxEvents" );
    skipEvents = inputInfo.getParameter<int>("skipEvents");
    runMin     = inputInfo.getParameter<int>("runMin"    );
    runMax     = inputInfo.getParameter<int>("runMax"    );
    inputFilenames = inputInfo.getParameter<vector<string> >("fileNames");

    // Set up output file, counts, and tree.
    outputFilename += TString(outputInfo.getParameter<string>("label")) + TString(".root");
    cout << "  Creating output file: " << outputFilename << endl;
    outputFile = new TFile(outputFilename, "RECREATE");

    nEventsProcessed = new TH1F("nEventsProcessed", "# Events Processed", 1, 0, 2 );

    zjCandTree = new TTree("zjCandidates", "Z+j Candidate Events");
    zjCandTree->Branch("event",  &event_,  "run/I:lumi/I:event/I:json/I" );

    zjCandTree->Branch("nJets",  &nJets_,  "nJets_/I"  );
    zjCandTree->Branch("jet_pt",		jetInfo_.pt,					"pt[nJets_]/F");
    zjCandTree->Branch("jet_eta",		jetInfo_.eta,				"eta[nJets_]/F");
    zjCandTree->Branch("jet_phi",		jetInfo_.phi,				"phi[nJets_]/F");
    zjCandTree->Branch("jet_e",			jetInfo_.e,					"e[nJets_]/F");
    zjCandTree->Branch("jet_csv",		jetInfo_.csv, 				"csv[nJets_]/F");
    zjCandTree->Branch("jet_csv_nominal",	jetInfo_.csv_nominal, 		"csv_nominal[nJets_]/F");
    zjCandTree->Branch("jet_csv_upBC",		jetInfo_.csv_upBC, 			"csv_upBC[nJets_]/F");
    zjCandTree->Branch("jet_csv_downBC",	jetInfo_.csv_downBC, 		"csv_downBC[nJets_]/F");
    zjCandTree->Branch("jet_csv_upL",		jetInfo_.csv_upL, 			"csv_upL[nJets_]/F");
    zjCandTree->Branch("jet_csv_downL",		jetInfo_.csv_downL, 			"csv_downL[nJets_]/F");
    zjCandTree->Branch("jet_csv_nominal4p",	jetInfo_.csv_nominal4p, 		"csv_nominal4p[nJets_]/F");
    zjCandTree->Branch("jet_csv_upBC4p",	jetInfo_.csv_upBC4p, 		"csv_upBC4p[nJets_]/F");
    zjCandTree->Branch("jet_csv_downBC4p",	jetInfo_.csv_downBC4p, 		"csv_downBC4p[nJets_]/F");
    zjCandTree->Branch("jet_csv_nominal1Bin",   jetInfo_.csv_nominal1Bin, 	"csv_nominal1Bin[nJets_]/F");
    zjCandTree->Branch("jet_csvivf",		jetInfo_.csvivf, 			"csvivf[nJets_]/F");
    zjCandTree->Branch("jet_cmva",		jetInfo_.cmva, 				"cmva[nJets_]/F");
    zjCandTree->Branch("jet_cosTheta",		jetInfo_.cosTheta, 			"cosTheta[nJets_]/F");
    zjCandTree->Branch("jet_numTracksSV",	jetInfo_.numTracksSV, 		"numTracksSV[nJets_]/I");
    zjCandTree->Branch("jet_chf",		jetInfo_.chf, 				"chf[nJets_]/F");
    zjCandTree->Branch("jet_nhf",		jetInfo_.nhf, 				"nhf[nJets_]/F");
    zjCandTree->Branch("jet_cef",		jetInfo_.cef, 				"cef[nJets_]/F");
    zjCandTree->Branch("jet_nef",		jetInfo_.nef, 				"nef[nJets_]/F");
    zjCandTree->Branch("jet_nch",		jetInfo_.nch, 				"nch[nJets_]/F");
    zjCandTree->Branch("jet_nconstituents",	jetInfo_.nconstituents,		"nconstituents[nJets_]");
    zjCandTree->Branch("jet_flavour",		jetInfo_.flavour, 			"flavour[nJets_]/F");
    zjCandTree->Branch("jet_isSemiLept",	jetInfo_.isSemiLept, 		"isSemiLept[nJets_]/I");
    zjCandTree->Branch("jet_isSemiLeptMCtruth",	jetInfo_.isSemiLeptMCtruth, 	"isSemiLeptMCtruth[nJets_]/I");
    zjCandTree->Branch("jet_SoftLeptpdgId",	jetInfo_.SoftLeptpdgId,  	"SoftLeptpdgId[nJets_]/I");
    zjCandTree->Branch("jet_SoftLeptIdlooseMu", jetInfo_.SoftLeptIdlooseMu,  	"SoftLeptIdlooseMu[nJets_]/I");
    zjCandTree->Branch("jet_SoftLeptId95", 	jetInfo_.SoftLeptId95, 		"SoftLeptId95[nJets_]/I");
    zjCandTree->Branch("jet_SoftLeptPt", 	jetInfo_.SoftLeptPt, 		"SoftLeptPt[nJets_]/F");
    zjCandTree->Branch("jet_SoftLeptdR", 	jetInfo_.SoftLeptdR, 		"SoftLeptdR[nJets_]/F");
    zjCandTree->Branch("jet_SoftLeptptRel", 	jetInfo_.SoftLeptptRel, 		"SoftLeptptRel[nJets_]/F");
    zjCandTree->Branch("jet_SoftLeptRelCombIso",jetInfo_.SoftLeptRelCombIso,  	"SoftLeptRelCombIso[nJets_]/F");
    zjCandTree->Branch("jet_puJetIdL", 		jetInfo_.puJetIdL,  			"puJetIdL[nJets_]/F");
    zjCandTree->Branch("jet_puJetIdM", 		jetInfo_.puJetIdM,  			"puJetIdM[nJets_]/F");
    zjCandTree->Branch("jet_puJetIdT", 		jetInfo_.puJetIdT,  			"puJetIdT[nJets_]/F");
    zjCandTree->Branch("jet_puJetIdMva", 	jetInfo_.puJetIdMva,  		"puJetIdMva[nJets_]/F");
    zjCandTree->Branch("jet_charge",		jetInfo_.charge, 			"charge[nJets_]/F");
    zjCandTree->Branch("jet_genPt",		jetInfo_.genPt, 				"genPt[nJets_]/F");
    zjCandTree->Branch("jet_genEta",		jetInfo_.genEta, 			"genEta[nJets_]/F");
    zjCandTree->Branch("jet_genPhi",		jetInfo_.genPhi, 			"genPhi[nJets_]/F");
    zjCandTree->Branch("jet_JECUnc",		jetInfo_.JECUnc, 			"JECUnc[nJets_]/F");
    zjCandTree->Branch("jet_vtxMass",		jetInfo_.vtxMass, 			"vtxMass[nJets_]/F");
    zjCandTree->Branch("jet_vtx3dL",		jetInfo_.vtx3dL, 			"vtx3dL[nJets_]/F");
    zjCandTree->Branch("jet_vtx3deL",		jetInfo_.vtx3deL,	 		"vtx3deL[nJets_]/F");
    zjCandTree->Branch("jet_id",		jetInfo_.id, 				"id[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVL",		jetInfo_.SF_CSVL, 			"SF_CSVL[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVM",		jetInfo_.SF_CSVM, 			"SF_CSVM[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVT",		jetInfo_.SF_CSVT, 			"SF_CSVT[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVLerr",	jetInfo_.SF_CSVLerr, 		"SF_CSVLerr[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVMerr",	jetInfo_.SF_CSVMerr, 		"SF_CSVMerr[nJets_]/b");
    zjCandTree->Branch("jet_SF_CSVTerr",	jetInfo_.SF_CSVTerr, 		"SF_CSVTerr[nJets_]/b");

    zjCandTree->Branch("nMuons", &nMuons_, "nMuons_/I" );
    zjCandTree->Branch("muon_mass",		    muonInfo_.mass,			"mass[nMuons_]/F");
    zjCandTree->Branch("muon_pt",		    muonInfo_.pt,			"pt[nMuons_]/F");
    zjCandTree->Branch("muon_eta",		    muonInfo_.eta,			"eta[nMuons_]");
    zjCandTree->Branch("muon_phi",		    muonInfo_.phi,			"phi[nMuons_]/F");
    zjCandTree->Branch("muon_aodCombRelIso",	    muonInfo_.aodCombRelIso,	"aodCombRelIso[nMuons_]/F");
    zjCandTree->Branch("muon_pfCombRelIso",	    muonInfo_.pfCombRelIso,	"pfCombRelIso[nMuons_]/F");
    zjCandTree->Branch("muon_photonIso",	    muonInfo_.photonIso,		"photonIso[nMuons_]/F");
    zjCandTree->Branch("muon_neutralHadIso",	    muonInfo_.neutralHadIso,	"neutralHadIso[nMuons_]/F");
    zjCandTree->Branch("muon_chargedHadIso",	    muonInfo_.chargedHadIso,	"chargedHadIso[nMuons_]/F");
    zjCandTree->Branch("muon_chargedPUIso",	    muonInfo_.chargedPUIso,	"chargedPUIso[nMuons_]/F");
    zjCandTree->Branch("muon_particleIso",	    muonInfo_.particleIso,	"particleIso[nMuons_]/F");
    zjCandTree->Branch("muon_dxy",		    muonInfo_.dxy,			"dxy[nMuons_]/F");
    zjCandTree->Branch("muon_dz",		    muonInfo_.dz,			"dz[nMuons_]/F");
    zjCandTree->Branch("muon_type",		    muonInfo_.type,			"type[nMuons_]/I");
    zjCandTree->Branch("muon_id80",		    muonInfo_.id80,			"id80[nMuons_]/F");
    zjCandTree->Branch("muon_id95",		    muonInfo_.id95,			"id95[nMuons_]/F");
    zjCandTree->Branch("muon_vbtf",		    muonInfo_.vbtf,			"vbtf[nMuons_]/F");
    zjCandTree->Branch("muon_id80NoIso",	    muonInfo_.id80NoIso,		"id80NoIso[nMuons_]/F");
    zjCandTree->Branch("muon_genPt",		    muonInfo_.genPt,			"genPt[nMuons_]/F");
    zjCandTree->Branch("muon_genEta",		    muonInfo_.genEta,		"genEta[nMuons_]");
    zjCandTree->Branch("muon_genPhi",		    muonInfo_.genPhi,		"genPhi[nMuons_]/F");
    zjCandTree->Branch("muon_charge",		    muonInfo_.charge,		"charge[nMuons_]/F");
    zjCandTree->Branch("muon_pfCorrIso",            muonInfo_.pfCorrIso,		"pfCorrIso[nMuons_]/F");
    zjCandTree->Branch("muon_pfCorrIsoHCP",         muonInfo_.pfCorrIsoHCP,	"pfCorrIsoHCP[nMuons_]/F");
    zjCandTree->Branch("muon_id2012tight",          muonInfo_.id2012tight,	"id2012tight[nMuons_]/F");
    zjCandTree->Branch("muon_idMVAnotrig",          muonInfo_.idMVAnotrig,	"idMVAnotrig[nMuons_]/F");
    zjCandTree->Branch("muon_idMVAtrig",	    muonInfo_.idMVAtrig,		"idMVAtrig[nMuons_]/F");
    zjCandTree->Branch("muon_idMVApresel",	    muonInfo_.idMVApresel,	"idMVApresel[nMuons_]/F");
    zjCandTree->Branch("muon_innerHits",	    muonInfo_.innerHits,		"innerHits[nMuons_]/F");
    zjCandTree->Branch("muon_photonIsoDoubleCount", muonInfo_.photonIsoDoubleCount, "photonIsoDoubleCount[nMuons_]/F");
    zjCandTree->Branch("muon_wpHWW",		    muonInfo_.wpHWW,		    "wpHWW[nMuons_]/F");
    zjCandTree->Branch("muon_wp95",		    muonInfo_.wp95,	            "wp95[nMuons_]/F");
    zjCandTree->Branch("muon_wp90",		    muonInfo_.wp90,		    "wp90[nMuons_]/F");
    zjCandTree->Branch("muon_wp85",		    muonInfo_.wp85,		    "wp85[nMuons_]/F");
    zjCandTree->Branch("muon_wp80",		    muonInfo_.wp80,		    "wp80[nMuons_]/F");
    zjCandTree->Branch("muon_wp70",		    muonInfo_.wp70,		    "wp70[nMuons_]/F");

    zjCandTree->Branch("nElecs", &nElecs_, "nElecs_/I" );
    zjCandTree->Branch("elec_mass",			elecInfo_.mass,			"mass[nElecs_]/F");
    zjCandTree->Branch("elec_pt",			elecInfo_.pt,			"pt[nElecs_]/F");
    zjCandTree->Branch("elec_eta",			elecInfo_.eta,			"eta[nElecs_]");
    zjCandTree->Branch("elec_phi",			elecInfo_.phi,			"phi[nElecs_]/F");
    zjCandTree->Branch("elec_aodCombRelIso",	elecInfo_.aodCombRelIso,	"aodCombRelIso[nElecs_]/F");
    zjCandTree->Branch("elec_pfCombRelIso",	elecInfo_.pfCombRelIso,	"pfCombRelIso[nElecs_]/F");
    zjCandTree->Branch("elec_photonIso",		elecInfo_.photonIso,		"photonIso[nElecs_]/F");
    zjCandTree->Branch("elec_neutralHadIso",	elecInfo_.neutralHadIso,	"neutralHadIso[nElecs_]/F");
    zjCandTree->Branch("elec_chargedHadIso",	elecInfo_.chargedHadIso,	"chargedHadIso[nElecs_]/F");
    zjCandTree->Branch("elec_chargedPUIso",	elecInfo_.chargedPUIso,	"chargedPUIso[nElecs_]/F");
    zjCandTree->Branch("elec_particleIso",	elecInfo_.particleIso,	"particleIso[nElecs_]/F");
    zjCandTree->Branch("elec_dxy",			elecInfo_.dxy,			"dxy[nElecs_]/F");
    zjCandTree->Branch("elec_dz",			elecInfo_.dz,			"dz[nElecs_]/F");
    zjCandTree->Branch("elec_type",			elecInfo_.type,			"type[nElecs_]/I");
    zjCandTree->Branch("elec_id80",			elecInfo_.id80,			"id80[nElecs_]/F");
    zjCandTree->Branch("elec_id95",			elecInfo_.id95,			"id95[nElecs_]/F");
    zjCandTree->Branch("elec_vbtf",			elecInfo_.vbtf,			"vbtf[nElecs_]/F");
    zjCandTree->Branch("elec_id80NoIso",		elecInfo_.id80NoIso,		"id80NoIso[nElecs_]/F");
    zjCandTree->Branch("elec_genPt",			elecInfo_.genPt,			"genPt[nElecs_]/F");
    zjCandTree->Branch("elec_genEta",		elecInfo_.genEta,		"genEta[nElecs_]");
    zjCandTree->Branch("elec_genPhi",		elecInfo_.genPhi,		"genPhi[nElecs_]/F");
    zjCandTree->Branch("elec_charge",		elecInfo_.charge,		"charge[nElecs_]/F");
    zjCandTree->Branch("elec_pfCorrIso",		elecInfo_.pfCorrIso,		"pfCorrIso[nElecs_]/F");
    zjCandTree->Branch("elec_pfCorrIsoHCP",	elecInfo_.pfCorrIsoHCP,	"pfCorrIsoHCP[nElecs_]/F");
    zjCandTree->Branch("elec_id2012tight",	elecInfo_.id2012tight,	"id2012tight[nElecs_]/F");
    zjCandTree->Branch("elec_idMVAnotrig",	elecInfo_.idMVAnotrig,	"idMVAnotrig[nElecs_]/F");
    zjCandTree->Branch("elec_idMVAtrig",		elecInfo_.idMVAtrig,		"idMVAtrig[nElecs_]/F");
    zjCandTree->Branch("elec_idMVApresel",	elecInfo_.idMVApresel,	"idMVApresel[nElecs_]/F");
    zjCandTree->Branch("elec_innerHits",		elecInfo_.innerHits,		"innerHits[nElecs_]/F");
    zjCandTree->Branch("elec_photonIsoDoubleCount",elecInfo_.photonIsoDoubleCount,"photonIsoDoubleCount[nElecs_]/F");
    zjCandTree->Branch("elec_wpHWW",			elecInfo_.wpHWW,			"wpHWW[nElecs_]/F");
    zjCandTree->Branch("elec_wp95",			elecInfo_.wp95,			"wp95[nElecs_]/F");
    zjCandTree->Branch("elec_wp90",			elecInfo_.wp90,			"wp90[nElecs_]/F");
    zjCandTree->Branch("elec_wp85",			elecInfo_.wp85,			"wp85[nElecs_]/F");
    zjCandTree->Branch("elec_wp80",			elecInfo_.wp80,			"wp80[nElecs_]/F");
    zjCandTree->Branch("elec_wp70",			elecInfo_.wp70,			"wp70[nElecs_]/F");




    // --- Loop over all files and events --- //
    for( unsigned int curFile=0; curFile<inputFilenames.size(); curFile++)
    {
        // Open input files
        cout << "  Opening file " << curFile << ": " << inputFilenames[curFile] << endl;
        inputFile = TFile::Open(inputFilenames[curFile].c_str());
        if(inputFile==0) { cout << "  ERROR: Unable to open file: " << inputFilenames[curFile] << endl; continue;}

        // Loop over events
        fwlite::Event curEvent(inputFile);
        for(curEvent.toBegin(); !curEvent.atEnd(); ++curEvent, nEventsLoopedOver++)
        {
            // Extract event, run, lumi information
            event_.run   = curEvent.id().run();
            event_.lumi  = curEvent.id().luminosityBlock();
            event_.event = curEvent.id().event();
            event_.json  = jsonContainsEvent(jsonVector, curEvent);

            // Check to see:
            //   - if enough events have been skipped
            //   - if the maximum have been processed
            //   - if the event is out of the run range
            if(nEventsLoopedOver < skipEvents) continue;
            if(maxEvents>=0 && nEventsLoopedOver > maxEvents+skipEvents) break;
            if(runMin>0 && event_.run < runMin) continue;
            if(runMax>0 && event_.run > runMax) continue;

            nEventsProcessed->Fill(1.);
            
            // Get HbbAnalyzer Object from Step1 PATTuble
            const char * lab = "HbbAnalyzerNew";
            fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
            vhbbAuxHandle.getByLabel(curEvent,lab,0,0);
            const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();

            // Get VHbbEvent object from file
            fwlite::Handle< VHbbEvent > vhbbHandle;
            vhbbHandle.getByLabel(curEvent,"HbbAnalyzerNew");
            VHbbEvent vhbbEvent = *vhbbHandle.product();        // RAW EVENT. Actual Ntupler has JEC corrections, etc.

            // Extract Jet, Muon information from VHbbEvent object
            vector<VHbbEvent::MuonInfo>     &rawMuons = vhbbEvent.muInfo;
            vector<VHbbEvent::ElectronInfo> &rawElecs = vhbbEvent.eleInfo;
            vector<VHbbEvent::SimpleJet>    &rawJets  = vhbbEvent.simpleJets2;

            jetInfo_.reset(); muonInfo_.reset(); elecInfo_.reset();
            nJets_  = rawJets.size();
            nMuons_ = rawMuons.size();
            nElecs_ = rawElecs.size();
            for(int i=0; i<nJets_  && i<MAXJ; i++)  jetInfo_.set(rawJets[i],i);
            for(int i=0; i<nMuons_ && i<MAXL; i++) muonInfo_.set(rawMuons[i],i,13,aux);
            for(int i=0; i<nElecs_ && i<MAXL; i++) elecInfo_.set(rawElecs[i],i,11,aux);

            zjCandTree->Fill();
        }

        inputFile->Close();
    }


    // Write to output and close
    //outputFile->cd();
    //zjCandTree->Write();
    outputFile->Write();
    outputFile->Close();
    
    cout << "\n  ZCNtupler Completed Successfully at " << timestamp() << "."
         << "\n  ====== ZCNtupler (COMPLETE) ======\n"   << endl;


return 0;
}



