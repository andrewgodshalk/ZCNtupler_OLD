/*-----------------------------------------------------------------------------
ZCHelperFunctions
  Created:  2014-03-05 by Andrew Godshalk (godshalk@buffalo.edu)
  Modified: 2014-03-05 by Andrew Godshalk

Part of a Step-2 Ntupler, working from PATTuples created by the VHbb group.
Contains functions and structs used in the Ntupler.

Functions:
    string timestamp()
    string file_timestamp()
    bool jsonContainsEvent()

Structs:
    EventInfo
    LeptonInfo
    xxxJetInfo  - made into independent h-file, into a class.

Modifications:
2014-03-05 - Included timestamp functions from previous work on Lorentz Angle.

-----------------------------------------------------------------------------*/


#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagReshaping.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include <sstream>
#include <string>
#include <vector>

#define MAXJ 130
#define MAXL 110


using namespace std;

// --- [FUNCTION DECLARATIONS] --- //
  string timestamp(char d = '-', char b = ' ', char t = ':');
    // Outputs string w/ date and time, by default in YYYY-MM-DD HH:MM:SS format, for display
    // Input variables are for separation of date values, data and time, and time values.
  string file_timestamp();
    // Outputs string w/ date and time, in YYYY-MM-DD_HHMMSS format, for filenames
  bool jsonContainsEvent (const vector<edm::LuminosityBlockRange> &jsonVec, const edm::EventBase &event);
    // Checks JSON vector to see if it contains the event

// --- [GLOBAL VARIABLES] --- //
  BTagShapeInterface * nominalShape=0;
  BTagShapeInterface * downBCShape=0;
  BTagShapeInterface * upBCShape=0;
  BTagShapeInterface * downLShape=0;
  BTagShapeInterface * upLShape=0; 
  BTagShapeInterface * nominalShape4p=0;
  BTagShapeInterface * downBCShape4p=0;
  BTagShapeInterface * upBCShape4p=0;
  BTagShapeInterface * nominalShape1Bin=0;

// --- [STRUCT DEFINITIONS] --- //
  struct EventInfo { int run, lumi, event, json; };
  
  struct LeptonInfo
  {
    void reset()
    {
      for(int i =0; i < MAXL;i++)
      { 
        mass[i]=-99; pt[i]=-99; eta[i]=-99; phi[i]=-99; aodCombRelIso[i]=-99; pfCombRelIso[i]=-99; photonIso[i]=-99; neutralHadIso[i]=-99; chargedHadIso[i]=-99; chargedPUIso[i]=-99; particleIso[i]=-99; dxy[i]=-99; dz[i]=-99; type[i]=-99;  genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99;  
        id80[i]=-99; id95[i]=-99; vbtf[i]=-99; id80NoIso[i]=-99;
        charge[i]=-99;wp70[i]=-99; wp80[i]=-99;wp85[i]=-99;wp90[i]=-99;wp95[i]=-99;wpHWW[i]=-99;
        pfCorrIso[i]=-99.; id2012tight[i]=-99; idMVAnotrig[i]=-99; idMVAtrig[i]=-99;innerHits[i]=-99.,pfCorrIsoHCP[i]=-99.;
  
        decayModeFinding[i] =-99.;
        byLooseCombinedIsolationDeltaBetaCorr[i] =-99.;
        againstMuonTight[i] =-99.;
        againstElectronLoose[i] =-99.;
        againstElectronMedium[i] =-99.;
        againstElectronMVA[i] =-99.;
        NsignalPFChargedHadrCands[i] =-99;
        NsignalPFGammaCands[i] =-99;
        leadPFChargedHadrCandPt[i] = -99.;
        byLooseIsolation[i]  = -99.;
        byMediumIsolation[i] = -99.;
        byTightIsolation[i] = -99.;
        byLooseCombinedIsolationDeltaBetaCorr3Hits[i] = -99.;
        byMediumCombinedIsolationDeltaBetaCorr3Hits[i] = -99.;
        byTightCombinedIsolationDeltaBetaCorr3Hits[i] = -99.;
        againstElectronMVA3raw[i] = -99.;
        againstElectronMVA3category[i] = -99.;
        againstElectronLooseMVA3[i] = -99.;
        againstElectronMediumMVA3[i] = -99.;
        againstElectronTightMVA3[i] = -99.;
        againstElectronVTightMVA3[i] = -99.;
        againstElectronDeadECAL[i] = -99.;
        byLooseIsolationMVA[i] = -99.;
        byMediumIsolationMVA[i] = -99.;
        byTightIsolationMVA[i] = -99.;
        byLooseIsolationMVA2[i] = -99.;
        byMediumIsolationMVA2[i] = -99.;
        byTightIsolationMVA2[i] = -99.;
        againstMuonLoose2[i] = -99.;
        againstMuonMedium2[i] = -99.;
        againstMuonTight2[i] = -99.;
      }
    }
  
    template <class Input> void set(const Input & i, int j,int t,const VHbbEventAuxInfo & aux)
    {
      type[j]=t;
      pt[j]=i.p4.Pt(); 
      mass[j]=i.p4.M();
      eta[j]=i.p4.Eta();
      phi[j]=i.p4.Phi();
      aodCombRelIso[j]=(i.hIso+i.eIso+i.tIso)/i.p4.Pt();
      pfCombRelIso[j]=(i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt();
      photonIso[j]=i.pfPhoIso;
      neutralHadIso[j]=i.pfNeuIso;
      chargedHadIso[j]=i.pfChaIso;
      chargedPUIso[j]=i.pfChaPUIso;
      charge[j]=i.charge;
      if(i.mcFourMomentum.Pt() > 0)
      { 
        genPt[j]=i.mcFourMomentum.Pt();
        genEta[j]=i.mcFourMomentum.Eta();
        genPhi[j]=i.mcFourMomentum.Phi();
      }
      setSpecific(i,j,aux);
    }
    
    template <class Input> void setSpecific(const Input & i, int j,const VHbbEventAuxInfo & aux) {}      
  
    float mass[MAXL];  //MT in case of W
    float pt[MAXL];
    float eta[MAXL];
    float phi[MAXL];
    float aodCombRelIso[MAXL];
    float pfCombRelIso[MAXL];
    float photonIso[MAXL];
    float neutralHadIso[MAXL];
    float chargedHadIso[MAXL];
    float chargedPUIso[MAXL];
    float particleIso[MAXL];
    float genPt[MAXL];
    float genEta[MAXL];
    float genPhi[MAXL];
    float dxy[MAXL];
    float dz[MAXL];
    int type[MAXL];
    float id80[MAXL];
    float id95[MAXL];
    float vbtf[MAXL];
    float id80NoIso[MAXL];
    float charge[MAXL];
    float pfCorrIso[MAXL];
    float pfCorrIsoHCP[MAXL];
    float id2012tight[MAXL];
    float idMVAnotrig[MAXL];
    float idMVAtrig[MAXL];
    float idMVApresel[MAXL];
    float wp70[MAXL]; 
    float wp80[MAXL]; 
    float wp85[MAXL]; 
    float wp90[MAXL]; 
    float wp95[MAXL]; 
    float wpHWW[MAXL];
    float innerHits[MAXL]; 
    float photonIsoDoubleCount[MAXL];
  
    float decayModeFinding[MAXL],byLooseCombinedIsolationDeltaBetaCorr[MAXL],againstMuonTight[MAXL],againstElectronLoose[MAXL],againstElectronMedium[MAXL],againstElectronMVA[MAXL];
    int NsignalPFChargedHadrCands[MAXL], NsignalPFGammaCands[MAXL];
    float leadPFChargedHadrCandPt[MAXL], byLooseIsolation[MAXL], byMediumIsolation[MAXL], byTightIsolation[MAXL]; 
    float byLooseCombinedIsolationDeltaBetaCorr3Hits[MAXL],  byMediumCombinedIsolationDeltaBetaCorr3Hits[MAXL], byTightCombinedIsolationDeltaBetaCorr3Hits[MAXL], againstElectronMVA3raw[MAXL], againstElectronMVA3category[MAXL], againstElectronLooseMVA3[MAXL], againstElectronMediumMVA3[MAXL], againstElectronTightMVA3[MAXL], againstElectronVTightMVA3[MAXL], againstElectronDeadECAL[MAXL], byLooseIsolationMVA[MAXL], byMediumIsolationMVA[MAXL], byTightIsolationMVA[MAXL], byLooseIsolationMVA2[MAXL], byMediumIsolationMVA2[MAXL], byTightIsolationMVA2[MAXL], againstMuonLoose2[MAXL], againstMuonMedium2[MAXL], againstMuonTight2[MAXL];
  };
  
  struct JetInfo
  {
    // set() - Copy info from SimpleJet j onto JetInfo index i
    void set(const VHbbEvent::SimpleJet & j, int i) 
    {
      pt[i]=j.p4.Pt();
      eta[i]=j.p4.Eta();
      phi[i]=j.p4.Phi();
      e[i]=j.p4.E();
      csv[i]=j.csv;
      if(nominalShape)
      {
        csv_nominal[i]=nominalShape->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_upBC[i]=upBCShape->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_downBC[i]=downBCShape->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_upL[i]=upLShape->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_downL[i]=downLShape->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_nominal4p[i]=nominalShape4p->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_upBC4p[i]=upBCShape4p->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_downBC4p[i]=downBCShape4p->reshape(eta[i],pt[i],j.csv,j.flavour);
        csv_nominal1Bin[i]=nominalShape1Bin->reshape(eta[i],pt[i],j.csv,j.flavour);
      }
      else
      {
        csv_nominal[i]=csv[i];
        csv_downBC[i]=csv[i];
        csv_upBC[i]=csv[i];
        csv_downL[i]=csv[i];
        csv_upL[i]=csv[i];
        csv_nominal4p[i]=csv[i];
        csv_downBC4p[i]=csv[i];
        csv_upBC4p[i]=csv[i];
        csv_nominal1Bin[i]=csv[i];
      }
      csvivf[i]=j.csvivf;
      cmva[i]=j.cmva;
      numTracksSV[i] = j.vtxNTracks;
      vtxMass[i]= j.vtxMass;
      vtxPt[i]= j.vtxP4.Pt();
      if(j.vtxP4.Pt() > 0) vtxEta[i]= j.vtxP4.Eta(); else  vtxEta[i]=-99;
      vtxPhi[i]= j.vtxP4.Phi();
      vtxE[i]= j.vtxP4.E();
      vtx3dL[i] = j.vtx3dL;
      vtx3deL[i] = j.vtx3deL;
      chf[i]=j.chargedHadronEFraction;
      nhf[i]  =j.neutralHadronEFraction;
      cef[i]  =j.chargedEmEFraction;
      nef[i]  =j.neutralEmEFraction;
      nconstituents[i]  =j.nConstituents;
      nch[i]=j.ntracks;
      SF_CSVL[i]=j.SF_CSVL;
      SF_CSVM[i]=j.SF_CSVM;
      SF_CSVT[i]=j.SF_CSVT;
      SF_CSVLerr[i]=j.SF_CSVLerr;
      SF_CSVMerr[i]=j.SF_CSVMerr;
      SF_CSVTerr[i]=j.SF_CSVTerr;
  
      flavour[i]=j.flavour;
      isSemiLeptMCtruth[i]=j.isSemiLeptMCtruth;
      isSemiLept[i]=j.isSemiLept;
      SoftLeptpdgId[i] = j.SoftLeptpdgId;
      SoftLeptIdlooseMu[i] = j.SoftLeptIdlooseMu;
      SoftLeptId95[i] = j.SoftLeptId95;
      SoftLeptPt[i] = j.SoftLeptPt;
      SoftLeptdR[i] = j.SoftLeptdR;  
      SoftLeptptRel[i] = j.SoftLeptptRel;
      SoftLeptRelCombIso[i] = j.SoftLeptRelCombIso;
      if(j.bestMCp4.Pt() > 0)
      {
            genPt[i]=j.bestMCp4.Pt();
            genEta[i]=j.bestMCp4.Eta();
            genPhi[i]=j.bestMCp4.Phi();
      }
      JECUnc[i]=j.jecunc;
      id[i]=jetId(i);
      ptRaw[i]=j.ptRaw;
      ptLeadTrack[i]=j.ptLeadTrack; 
      puJetIdL[i]=j.puJetIdL; 
      puJetIdM[i]=j.puJetIdM; 
      puJetIdT[i]=j.puJetIdT; 
      puJetIdMva[i]=j.puJetIdMva; 
      charge[i]=j.charge;     
      jetArea[i]=j.jetArea;
    }
  
    bool jetId(int i)
    {
      if(nhf[i] > 0.99) return false;
      if(nef[i] > 0.99) return false;
      if(nconstituents[i]  <= 1) return false;
      if(fabs(eta[i])<2.5)
      {
        if(cef[i] > 0.99) return false;
        if(chf[i] == 0) return false;
        if(nch[i]== 0) return false;
      }
      return true;
    }
  
    void reset()
    {
      for(int i=0;i<MAXJ;i++)
      {
        pt[i]=-99; eta[i]=-99; phi[i]=-99;e[i]=-99;csv[i]=-99;csv_nominal[i]=-99.;csv_upBC[i]=-99.;csv_downBC[i]=-99.;csv_upL[i]=-99.;csv_downL[i]=-99.;csv_nominal1Bin[i]=-99.;
        csvivf[i]=-99; cmva[i]=-99;
        cosTheta[i]=-99; numTracksSV[i]=-99; chf[i]=-99; nhf[i]=-99; cef[i]=-99; nef[i]=-99; nch[i]=-99; nconstituents[i]=-99; flavour[i]=-99; isSemiLeptMCtruth[i]=-99; isSemiLept[i]=-99;      
        SoftLeptpdgId[i] = -99; SoftLeptIdlooseMu[i] = -99;  SoftLeptId95[i] =  -99;   SoftLeptPt[i] = -99;  SoftLeptdR[i] = -99;   SoftLeptptRel[i] = -99; SoftLeptRelCombIso[i] = -99;  
        genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99; JECUnc[i]=-99; ptRaw[i]=-99.; ptLeadTrack[i]=-99.; puJetIdL[i]=-99; puJetIdM[i]=-99; puJetIdT[i]=-99; puJetIdMva[i]=-99; charge[i]=-99; jetArea[i]=-99; selectedTauDR[i] = -99.;
      }
    }
  
    float pt[MAXJ];
    float eta[MAXJ];
    float phi[MAXJ];
    float e[MAXJ];
    float csv[MAXJ];
    float csv_nominal[MAXJ];
    float csv_upBC[MAXJ];
    float csv_downBC[MAXJ];
    float csv_upL[MAXJ];
    float csv_downL[MAXJ];
    float csv_nominal4p[MAXJ];
    float csv_upBC4p[MAXJ];
    float csv_downBC4p[MAXJ];
    float csv_nominal1Bin[MAXJ];
    float csvivf[MAXJ];
    float cmva[MAXJ];
    float cosTheta[MAXJ];
    int numTracksSV[MAXJ];
    float chf[MAXJ];
    float nhf[MAXJ];
    float cef[MAXJ];
    float nef[MAXJ];
    float nch[MAXJ];
    float nconstituents[MAXJ];
    float flavour[MAXJ];
    int isSemiLept[MAXJ];
    int isSemiLeptMCtruth[MAXJ];
    int SoftLeptpdgId[MAXJ] ;
    int SoftLeptIdlooseMu[MAXJ] ;  
    int SoftLeptId95[MAXJ]  ;   
    float SoftLeptPt[MAXJ] ;  
    float SoftLeptdR[MAXJ] ;   
    float SoftLeptptRel[MAXJ] ; 
    float SoftLeptRelCombIso[MAXJ];
    float genPt[MAXJ];
    float genEta[MAXJ];
    float genPhi[MAXJ];
    float JECUnc[MAXJ];
    float vtxMass[MAXJ];
    float vtxPt[MAXJ];
    float vtxEta[MAXJ];
    float vtxPhi[MAXJ];
    float vtxE[MAXJ];
    float vtx3dL [MAXJ];
    float vtx3deL[MAXJ];
    bool id[MAXJ];
    float SF_CSVL[MAXJ];
    float SF_CSVM[MAXJ];
    float SF_CSVT[MAXJ]; 
    float SF_CSVLerr[MAXJ];
    float SF_CSVMerr[MAXJ];
    float SF_CSVTerr[MAXJ];  
    float ptRaw[MAXJ];
    float ptLeadTrack[MAXJ];
    float puJetIdL[MAXJ];
    float puJetIdM[MAXJ];
    float puJetIdT[MAXJ];
    float puJetIdMva[MAXJ];
    float charge[MAXJ];
    float jetArea[MAXJ];
    float selectedTauDR[MAXJ];
  };

// --- [FUNCTION DEFINITIONS] --- //

  string timestamp( char d, char b, char t )
  {
    stringstream output;
    TDatime now;
    int date = now.GetDate();
    int time = now.GetTime();
    if( time < 0 )
    {
      time-= 240000;
      date+= 1;
    }
  
    output << date/10000 << d << (date%10000)/1000 << (date%1000)/100 << d << (date%100)/10 << date%10;
    if( b!='0' && t!='0' ) output << b << time/10000 << t << (time%10000)/1000 << (time%1000)/100 << t << (time%100)/10 << time%10;

  return output.str();
  }


  string file_timestamp()
  {
    stringstream output;
    TDatime now;
    int date = now.GetDate();
    int time = now.GetTime();
    if( time < 0 )
    {
      time-= 240000;
      date+= 1;
    }
  
    output << date/100000 << (date%100000)/10000 << '-' << (date%10000)/1000 << (date%1000)/100 << '-' << (date%100)/10 << date%10;
    output << '_' << time/100000 << (time%100000)/10000 << (time%10000)/1000 << (time%1000)/100 << (time%100)/10 << time%10;
  return output.str();
  }


  bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec, const edm::EventBase &event)
  {
   // if the jsonVec is empty, then no JSON file was provided so all events should pass    
    if (jsonVec.empty()) return true;
   // Create a pointer to an edm function to see if an event is in a Lumi Block.
    bool (*funcPtr) (edm::LuminosityBlockRange const &, edm::LuminosityBlockID const &) = &edm::contains;
    edm::LuminosityBlockID lumiID (event.id().run(), event.id().luminosityBlock());
    vector< edm::LuminosityBlockRange >::const_iterator iter =  find_if (jsonVec.begin(), jsonVec.end(), boost::bind(funcPtr, _1, lumiID) );
    return jsonVec.end() != iter;
  }






