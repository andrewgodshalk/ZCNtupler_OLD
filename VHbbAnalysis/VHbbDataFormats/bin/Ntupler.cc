//Dummy line
#include <TH1F.h>
#include <TH3F.h>
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
//#include "PhysicsTools/Utilities/interface/Lumi3DReWeighting.h"
#include <TH2F.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TSystem.h>
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "FWCore/FWLite/interface/AutoLibraryLoader.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "PhysicsTools/FWLite/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "DataFormats/Math/interface/deltaPhi.h"
#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/HbbCandidateFinderAlgo.h" 
#include "VHbbAnalysis/VHbbDataFormats/src/HbbCandidateFinderAlgo.cc"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEvent.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbEventAuxInfo.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/VHbbCandidate.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerReader.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/TopMassReco.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/JECFWLite.h"

//for IVF
#include "RecoBTag/SecondaryVertex/interface/SecondaryVertex.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/Measurement1D.h"
#include "DataFormats/GeometryVector/interface/GlobalVector.h"
#include "DataFormats/GeometryVector/interface/VectorUtil.h"
#include <DataFormats/GeometrySurface/interface/Surface.h>
#include "Math/SMatrix.h"

//for LHE info
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"

//Move class definition to Ntupler.h ?
//#include "VHbbAnalysis/VHbbDataFormats/interface/Ntupler.h"

#include "ZSV/BAnalysis/interface/SimBHadron.h"
//btagging
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagWeight.h"
#include "VHbbAnalysis/VHbbDataFormats/interface/BTagReshaping.h"
//trigger
#include "VHbbAnalysis/VHbbDataFormats/interface/TriggerWeight.h"

//Gen particles
#include "DataFormats/Candidate/interface/Candidate.h"

//GenJet
#include "SimDataFormats/JetMatching/interface/JetFlavour.h"
#include "SimDataFormats/JetMatching/interface/JetFlavourMatching.h"

#include "DataFormats/Common/interface/View.h"
#include "DataFormats/Common/interface/Ref.h"
#include "DataFormats/Common/interface/getRef.h"

#include "DataFormats/JetReco/interface/Jet.h"
//#include "DataFormats/JetReco/interface/GenJetCollection.h"

#include <sstream>
#include <string>
#include <fstream>
#include <iostream>
#include <algorithm>
#include <vector>       // std::vector

#define MAXPDF 100
#define MAXJ 130
#define MAXL 110
#define MAXB 110
#define MAXT 160
#define nMetUnc 24 
#define MAXGENOBJ 20

//eleEnDown/Up, muEn, tauEn, JES, JER, Unclustered for type1 MET [0-11] and than type1p2 MET [12-23]

struct CompareDeltaR
{
  CompareDeltaR(TLorentzVector p4dir_): p4dir(p4dir_) {}
  bool operator()( const VHbbEvent::SimpleJet& j1, const VHbbEvent::SimpleJet& j2 ) const
  {
    return j1.p4.DeltaR(p4dir) > j2.p4.DeltaR(p4dir);
  }
  TLorentzVector p4dir;
};


const GlobalVector flightDirection(const TVector3 pv, const reco::Vertex &sv)
{
  GlobalVector fdir(sv.position().X() - pv.X(),
                    sv.position().Y() - pv.Y(),
                    sv.position().Z() - pv.Z());
  return fdir;
}


struct RLE
{
  RLE(){}
  RLE(
    unsigned long run_,
    unsigned long lumi_,
    unsigned long event_
  ) : run(run_),lumi(lumi_),event(event_) {}
  unsigned long run;
  unsigned long lumi;
  unsigned long event;
};


bool evcomp (RLE i,RLE j)
{
  if(i.run  != j.run)  return (i.run  < j.run );
  if(i.lumi != j.lumi) return (i.lumi < j.lumi);
  return (i.event<j.event);

}


float weightNLOEWKsignal(float pt)
{
 if(pt < 50) return 1;
 return 0.94-(0.2-0.068)/400.*(pt-50.);
}


float weightNNLOQCDsignal(float pt,int njets)
{
/*      =      0.99532   +/-   0.0705809   
p1                        = -0.000262157   +/-   0.000254908 
*/

/*
p0                        =       1.0064   +/-   0.353707    
p1                        =  0.000288831   +/-   0.00210877  
*/
 //  SF for jet veto
if(njets <=2)   return 0.99532-pt*0.000262157;
  

return 1.0064+pt*0.0002888;

}


void  readBadEvents(const char * filename, vector<RLE> & badEvents)
{
  string str ;
  RLE temp;
  ifstream myfile (filename);
  if (myfile.is_open())
  {
    while ( myfile.good() )
    {
      getline (myfile,str);
      replace( str.begin(), str.end(), ':', ' ' ) ;
      istringstream stm(str) ;
      int a, b, c ;
      stm >> a >> b >> c ;
      temp.run=a; temp.event=c; temp.lumi=b;
      badEvents.push_back(temp);
  //    cout << a << '\t' << b << '\t' << c << '\n' ; 
   }
 }

/*d::cout << "Presort" << std::endl;
std::cout <<  std::binary_search (badEvents.begin(), badEvents.end(), RLE(1,1,1),evcomp) << endl;
std::cout <<  std::binary_search (badEvents.begin(), badEvents.end(), RLE(190456,59,3971011),evcomp) << endl;
std::cout << "sort" << std::endl;*/
std::sort (badEvents.begin(), badEvents.end(), evcomp);
std::cout << "sorted" << std::endl;
/*d::cout <<  std::binary_search (badEvents.begin(), badEvents.end(), RLE(1,1,1),evcomp) << endl;
std::cout <<  std::binary_search (badEvents.begin(), badEvents.end(), RLE(190456,59,3971011),evcomp) << endl;
*/
 myfile.close();
}


bool jsonContainsEvent (const std::vector< edm::LuminosityBlockRange > &jsonVec,
                        const edm::EventBase &event)
{
  // if the jsonVec is empty, then no JSON file was provided so all
  // events should pass
  if (jsonVec.empty())
    {
      return true;
    }
  bool (* funcPtr) (edm::LuminosityBlockRange const &,
		    edm::LuminosityBlockID const &) = &edm::contains;
  edm::LuminosityBlockID lumiID (event.id().run(), 
				 event.id().luminosityBlock());
  std::vector< edm::LuminosityBlockRange >::const_iterator iter = 
    std::find_if (jsonVec.begin(), jsonVec.end(),
		  boost::bind(funcPtr, _1, lumiID) );
  return jsonVec.end() != iter;

}


float resolutionBias(float eta)
{
// return 0;//Nominal!
/* if(eta< 1.1) return 0.05;
 if(eta< 2.5) return 0.10;
 if(eta< 5) return 0.30;*/
//new numbers from:	https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
//==These are for 7 TeV==
//   if(eta< 0.5) return 0.052;
//   if(eta< 1.1) return 0.057;
//   if(eta< 1.7) return 0.096;
//   if(eta< 2.3) return 0.134;
//   if(eta< 5) return 0.28;
//==NOTE : update for 8 TeV data (Duong 10-11-2015)
//put here c-1 with c from https://twiki.cern.ch/twiki/bin/view/CMS/JetResolution
   if(eta< 0.5) return 0.079;
   if(eta< 1.1) return 0.099;
   if(eta< 1.7) return 0.121;
   if(eta< 2.3) return 0.208;
   if(eta< 2.8) return 0.254;
   if(eta< 3.2) return 0.395;
   if(eta<5) return 0.056;

 return 0;
}


#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
namespace LHAPDF {
  void initPDFSet(int nset, int setid, int member=0);
  void initPDFSet(int nset, const std::string& filename, int member=0);
  int numberPDF(int nset);
  void usePDFMember(int nset, int member);
  double xfx(int nset, double x, double Q, int fl);
  double getXmin(int nset, int member);
  double getXmax(int nset, int member);
  double getQ2min(int nset, int member);
  double getQ2max(int nset, int member);
  void extrapolate(bool extrapolate=true);
}


typedef struct
{
  void set(const SimBHadron & sbhc, int i){
    mass[i] = sbhc.mass();
    pt[i] = sbhc.pt();
    eta[i] = sbhc.eta();
    phi[i] = sbhc.phi();
    vtx_x[i] = sbhc.decPosition.x();
    vtx_y[i] = sbhc.decPosition.y();
    vtx_z[i] = sbhc.decPosition.z();
    pdgId[i] = sbhc.pdgId();
    status[i] = sbhc.status();
  };
  void reset(){
    for(int i=0; i < MAXB; ++i){
      mass[i] = -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; vtx_x[i] = -99; vtx_y[i] = -99; vtx_z[i] = -99; pdgId[i] = -99; status[i] = -99;
    }
  };
  float mass[MAXB];
  float pt[MAXB];
  float eta[MAXB];
  float phi[MAXB];
  float vtx_x[MAXB];
  float vtx_y[MAXB];
  float vtx_z[MAXB];
  int pdgId[MAXB];
  int status[MAXB];
//   int quarkStatus[MAXB];
//   int brotherStatus[MAXB];
//   int otherId[MAXB];
//   bool etaOk[MAXB];
//   bool simOk[MAXB];
//   bool trackOk[MAXB];
//   bool cutOk[MAXB];
//   bool cutNewOk[MAXB];
//   bool mcMatchOk[MAXB];
//   bool matchOk[MAXB];
}
SimBHadronInfo;


typedef struct
{
  void set( const reco::SecondaryVertex & recoSv, const TVector3 recoPv, int isv){
    pt[isv]   = recoSv.p4().Pt();
    eta[isv]  = flightDirection(recoPv,recoSv).eta();
    phi[isv]  = flightDirection(recoPv,recoSv).phi();
    massBcand[isv] = recoSv.p4().M();
    massSv[isv] = recoSv.p4().M();
    dist3D[isv] = recoSv.dist3d().value();  
    distSig3D[isv] = recoSv.dist3d().significance(); 
    dist2D[isv] = recoSv.dist2d().value();  
    distSig2D[isv] = recoSv.dist2d().significance(); 
    dist3D_norm[isv] = recoSv.dist3d().value()/recoSv.p4().Gamma(); 
  };
  void reset(){
    for(int i = 0; i < MAXB; ++i){
      massBcand[i] = -99; massSv[i]= -99; pt[i] = -99; eta[i] = -99; phi[i] = -99; dist3D[i] = -99; distSig3D[i] = -99; dist2D[i] = -99; distSig2D[i] = -99; dist3D_norm[i] = -99;
    }
  };
  float massBcand[MAXB];
  float massSv[MAXB];
  float pt[MAXB]; 
  float eta[MAXB]; 
  float phi[MAXB];
  float dist3D[MAXB]; 
  float distSig3D[MAXB];
  float dist2D[MAXB]; 
  float distSig2D[MAXB];
  float dist3D_norm[MAXB];
}
IVFInfo;


typedef struct 
{
  int HiggsFlag;
  float mass; 
  float pt;
  float eta;
  float phi;
  float dR;
  float dPhi;
  float dEta;
}
HiggsInfo;

 
typedef struct
{
  int FatHiggsFlag; 
  float mass;
  float pt;
  float eta;
  float phi;
  float filteredmass;
  float filteredpt; 
  float filteredeta;
  float filteredphi;  
//  float dR;
//  float dPhi;
//  float dEta;
}
FatHiggsInfo;

 
typedef struct 
{
  float mass; 
  float pt;
  float eta;
  float phi;
  float status;
  float charge;
  float momid;
}
genParticleInfo;


typedef struct 
{
  float bmass; 
  float bpt;
  float beta;
  float bphi;
  float bstatus;
  float wdau1mass; 
  float wdau1pt;
  float wdau1eta;
  float wdau1phi;
  float wdau1id;
  float wdau2mass; 
  float wdau2pt;
  float wdau2eta;
  float wdau2phi;
  float wdau2id;  
}
genTopInfo;


typedef struct
{
  bool HiggsCSVtkSharing;
  bool HiggsIPtkSharing;
  bool HiggsSVtkSharing;
  bool FatHiggsCSVtkSharing;
  bool FatHiggsIPtkSharing;
  bool FatHiggsSVtkSharing;
}
TrackSharingInfo;


typedef struct 
{
  float mass;  //MT in case of W
  float pt;
  float eta;
  float phi;
}
TrackInfo;

  
struct  LeptonInfo
{
  void reset()
  {
    for(int i =0; i < MAXL;i++)
	{ 
      mass[i]=-99; pt[i]=-99; eta[i]=-99; phi[i]=-99; aodCombRelIso[i]=-99; pfCombRelIso[i]=-99; photonIso[i]=-99; neutralHadIso[i]=-99; chargedHadIso[i]=-99; chargedPUIso[i]=-99; particleIso[i]=-99; dxy[i]=-99; dz[i]=-99; type[i]=-99;  genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99;  
      id80[i]=-99; id95[i]=-99; vbtf[i]=-99; id80NoIso[i]=-99;
      charge[i]=-99;wp70[i]=-99; wp80[i]=-99;wp85[i]=-99;wp90[i]=-99;wp95[i]=-99;wpHWW[i]=-99;
      pfCorrIso[i]=-99.; id2012tight[i]=-99; idMVAnotrig[i]=-99; idMVAtrig[i]=-99;innerHits[i]=-99.,pfCorrIsoHCP[i]=-99.;
      tIso      [i]=-99.;
      eIso      [i]=-99.;
      hIso      [i]=-99.;
      pfChaIso  [i]=-99.;
      pfChaPUIso[i]=-99.;
      pfPhoIso  [i]=-99.;
      pfNeuIso  [i]=-99.;

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
      doBelongToVhCand[i] = 0 ;
    }
  }

  template <class Input> void set(const Input & i, int j,int t, const VHbbEventAuxInfo & aux, int doBelongToVhCandIn)
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
    tIso      [j]=i.tIso;
    eIso      [j]=i.eIso;
    hIso      [j]=i.hIso;
    pfChaIso  [j]=i.pfChaIso;
    pfChaPUIso[j]=i.pfChaPUIso;
    pfPhoIso  [j]=i.pfPhoIso;
    pfNeuIso  [j]=i.pfNeuIso;
    doBelongToVhCand[j] = doBelongToVhCandIn ;

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
  float tIso      [MAXL];
  float eIso      [MAXL];
  float hIso      [MAXL];
  float pfChaIso  [MAXL];
  float pfChaPUIso[MAXL];
  float pfPhoIso  [MAXL];
  float pfNeuIso  [MAXL];

  float decayModeFinding[MAXL],byLooseCombinedIsolationDeltaBetaCorr[MAXL],againstMuonTight[MAXL],againstElectronLoose[MAXL],againstElectronMedium[MAXL],againstElectronMVA[MAXL];
  int NsignalPFChargedHadrCands[MAXL], NsignalPFGammaCands[MAXL];
  float leadPFChargedHadrCandPt[MAXL], byLooseIsolation[MAXL], byMediumIsolation[MAXL], byTightIsolation[MAXL]; 
  float byLooseCombinedIsolationDeltaBetaCorr3Hits[MAXL],  byMediumCombinedIsolationDeltaBetaCorr3Hits[MAXL], byTightCombinedIsolationDeltaBetaCorr3Hits[MAXL], againstElectronMVA3raw[MAXL], againstElectronMVA3category[MAXL], againstElectronLooseMVA3[MAXL], againstElectronMediumMVA3[MAXL], againstElectronTightMVA3[MAXL], againstElectronVTightMVA3[MAXL], againstElectronDeadECAL[MAXL], byLooseIsolationMVA[MAXL], byMediumIsolationMVA[MAXL], byTightIsolationMVA[MAXL], byLooseIsolationMVA2[MAXL], byMediumIsolationMVA2[MAXL], byTightIsolationMVA2[MAXL], againstMuonLoose2[MAXL], againstMuonMedium2[MAXL], againstMuonTight2[MAXL];
  int doBelongToVhCand[MAXL] ;
};

  
template <> void LeptonInfo::setSpecific<VHbbEvent::ElectronInfo>(const VHbbEvent::ElectronInfo & i, int j,const VHbbEventAuxInfo & aux)
{
  photonIsoDoubleCount[j]=i.pfPhoIsoDoubleCounted;
  id80[j]=i.id80;
  id95[j]=i.id95;
  id80NoIso[j]=(i.innerHits ==0 && !(fabs(i.convDist)<0.02 && fabs(i.convDcot)<0.02) && ((i.isEB && i.sihih<0.01 && fabs(i.Dphi)<0.06 && fabs(i.Deta)<0.004) || (i.isEE && i.sihih<0.03 && fabs(i.Dphi)<0.03  && fabs(i.Deta)<0.007)));
  innerHits[j]=i.innerHits;
  float mincor=0.0;
  float minrho=0.0;
  float rho = std::max(aux.puInfo.rho25Iso,minrho);
  float eta=i.p4.Eta();
  float areagamma=0.5;
  float areaNH=0.5;
  float areaComb=0.5;

  // new EA for electrons to use for Moriond 13
  if(fabs(eta) <= 1.0 ) {areaComb=0.21;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areaComb=0.21;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areaComb=0.11;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areaComb=0.14;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areaComb=0.18;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areaComb=0.19;}
  if(fabs(eta) > 2.4  ) {areaComb=0.26;}
  /*
  if(fabs(eta) <= 1.0 ) {areagamma=0.14; areaNH=0.044; areaComb=0.18;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.13; areaNH=0.065; areaComb=0.20;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.079; areaNH=0.068; areaComb=0.15;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.13; areaNH=0.057; areaComb=0.19;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.15; areaNH=0.058; areaComb=0.21;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.16; areaNH=0.061; areaComb=0.22;}
  if(fabs(eta) > 2.4  ) {areagamma=0.18; areaNH=0.11; areaComb=0.29;}
  */
  /*
  if(fabs(eta) <= 1.0 ) {areagamma=0.081; areaNH=0.024; areaComb=0.10;}
  if(fabs(eta) > 1.0 &&  fabs(eta) <= 1.479 ) {areagamma=0.084; areaNH=0.037; areaComb=0.12;}
  if(fabs(eta) > 1.479 &&  fabs(eta) <= 2.0 ) {areagamma=0.048; areaNH=0.037; areaComb=0.085;}
  if(fabs(eta) > 2.0 &&  fabs(eta) <= 2.2 ) {areagamma=0.089; areaNH=0.023; areaComb=0.11;}
  if(fabs(eta) > 2.2 &&  fabs(eta) <= 2.3 ) {areagamma=0.092; areaNH=0.023; areaComb=0.12;}
  if(fabs(eta) > 2.3 &&  fabs(eta) <= 2.4 ) {areagamma=0.097; areaNH=0.021; areaComb=0.12;}
  if(fabs(eta) > 2.4  ) {areagamma=0.11; areaNH=0.021; areaComb=0.13;}
  */

  float pho=i.pfPhoIso;
  if(i.innerHits>0) pho-=i.pfPhoIsoDoubleCounted;

  pfCorrIsoHCP[j] = (i.pfChaIso+ std::max(pho-rho*areagamma,mincor )+std::max(i.pfNeuIso-rho*areaNH,mincor))/i.p4.Pt(); 
  // new definition for Moriond13 
  pfCorrIso[j] = (i.pfChaIso+ std::max(pho+i.pfNeuIso-rho*areaComb,mincor))/i.p4.Pt();
  
  id2012tight[j] = fabs(i.dxy) < 0.02  &&fabs(i.dz) < 0.1  &&(
	(i.isEE  &&fabs(i.Deta) < 0.005 &&fabs(i.Dphi) < 0.02 &&i.sihih < 0.03  &&i.HoE < 0.10  &&fabs(i.fMVAVar_IoEmIoP) < 0.05) ||
    (i.isEB  &&fabs(i.Deta) < 0.004 &&fabs(i.Dphi) < 0.03 &&i.sihih < 0.01  &&i.HoE < 0.12  &&fabs(i.fMVAVar_IoEmIoP) < 0.05
  ));

  bool mvaPreSel =
  (
    (    i.isEE 
	  //&& fabs(electrons[it].Deta) < 0.009
	  //&& fabs(electrons[it].Dphi) < 0.1
	  && i.sihih < 0.03
	  && i.HoE < 0.10
	  && i.innerHits == 0
	  && i.tIso/i.p4.Pt() < 0.2
	  && i.eIso/i.p4.Pt() < 0.2
	  && i.hIso/i.p4.Pt() < 0.2
	) 
   ||
    (	 i.isEB
	  //&& fabs(electrons[it].Deta) < 0.007 &&
	  //&& fabs(electrons[it].Dphi) < 0.015 &&
	  && i.sihih < 0.01
	  && i.HoE < 0.12
	  && i.innerHits == 0
	  && i.tIso/i.p4.Pt() < 0.2
	  && i.eIso/i.p4.Pt() < 0.2
	  && i.hIso/i.p4.Pt() < 0.2) 
  );

  float id=i.mvaOutTrig;
  if(!mvaPreSel) id=0;
  float iso=pfCorrIso[j];
  wp70[j] =((fabs(eta) < 0.8 && id>0.977 && iso < 0.093) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.956 && iso < 0.095) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.966 && iso < 0.171));
  wp80[j] =((fabs(eta) < 0.8 && id>0.913 && iso < 0.105) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.964 && iso < 0.178) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.899 && iso < 0.150));
  wp85[j] =((fabs(eta) < 0.8 && id>0.929 && iso < 0.135) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.931 && iso < 0.159) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.805 && iso < 0.155));
  wp90[j] =((fabs(eta) < 0.8 && id>0.877 && iso < 0.177) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.794 && iso < 0.180) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.846 && iso < 0.244));
  wp95[j] =((fabs(eta) < 0.8 && id>0.858 && iso < 0.253) ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.425 && iso < 0.225) || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.759 && iso < 0.308));
  wpHWW[j]=((fabs(eta) < 0.8 && id>0.94  && iso < 0.15)  ||  (fabs(eta) >= 0.8 && fabs(eta) < 1.479 && id>0.85 && iso  < 0.15)  || (fabs(eta) >= 1.479 && fabs(eta) < 2.5 && id>0.92 && iso < 0.15));

  idMVAnotrig[j]=i.mvaOut;
  idMVAtrig[j]=i.mvaOutTrig;
  idMVApresel[j]=mvaPreSel;
}


template <> void LeptonInfo::setSpecific<VHbbEvent::MuonInfo>(const VHbbEvent::MuonInfo & i, int j,const VHbbEventAuxInfo & aux)
{
  dxy[j]=i.ipDb;
  dz[j]=i.zPVPt;
  vbtf[j]=( i.globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nHits > 10 && i.cat & 0x1 && i.cat & 0x2 && i.nMatches >=2 && i.ipDb<.2 &&
          (i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt()<.15  && fabs(i.p4.Eta())<2.4 && i.p4.Pt()>20 ) ;
 
  float mincor=0.0;
  float minrho=0.0;
  float NoverCh=0.5;
  float rhoN = std::max(aux.puInfo.rhoNeutral,minrho);
  float eta=i.p4.Eta();
  float area=0.5;
  if(fabs(eta)>0.0 && fabs(eta) <= 1.0) {area=0.674;}
  if(fabs(eta)>1.0 && fabs(eta) <= 1.5) {area=0.565;}
  if(fabs(eta)>1.5 && fabs(eta) <= 2.0) {area=0.442;}
  if(fabs(eta)>2.0 && fabs(eta) <= 2.2) {area=0.515;}
  if(fabs(eta)>2.2 && fabs(eta) <= 2.3) {area=0.821;}
  if(fabs(eta)>2.3 && fabs(eta) <= 2.4) {area=0.660;}
  pfCorrIsoHCP[j] = (i.pfChaIso+ std::max(i.pfPhoIso+i.pfNeuIso-rhoN*area,mincor))/i.p4.Pt();
  // Moriond13: using dbeta corrections
  pfCorrIso[j] = (i.pfChaIso+ std::max(i.pfPhoIso+i.pfNeuIso-NoverCh*i.pfChaPUIso,mincor))/i.p4.Pt();

  id2012tight[j]= i.isPF && i.globChi2<10 && i.nPixelHits>= 1 && i.globNHits != 0 && i.nValidLayers > 5 && (i.cat & 0x2) && i.nMatches >=2 && i.ipDb<.2;
}


template <> void LeptonInfo::setSpecific<VHbbEvent::TauInfo>(const VHbbEvent::TauInfo &i, int j,const VHbbEventAuxInfo & aux)
{
  //std::cout << "In tau set specific" << std::endl;

  pt[j]=i.p4.Pt(); 
  mass[j]=i.p4.M();
  eta[j]=i.p4.Eta();
  phi[j]=i.p4.Phi();
  pfCombRelIso[j]=(i.pfChaIso+i.pfPhoIso+i.pfNeuIso)/i.p4.Pt();
  charge[j]=i.charge;
  if(i.mcFourMomentum.Pt() > 0)
  { 
    genPt[j]=i.mcFourMomentum.Pt();
    genEta[j]=i.mcFourMomentum.Eta();
    genPhi[j]=i.mcFourMomentum.Phi();
  }

  decayModeFinding[j] = i.decayModeFinding;
  byLooseCombinedIsolationDeltaBetaCorr[j] =i.byLooseCombinedIsolationDeltaBetaCorr;
  againstMuonTight[j] =i.againstMuonTight;
  againstElectronLoose[j] =i.againstElectronLoose;
  againstElectronMedium[j] =i.againstElectronMedium;
  againstElectronMVA[j] =i.againstElectronMVA;
  NsignalPFChargedHadrCands[j] =i.NsignalPFChargedHadrCands;
  NsignalPFGammaCands[j] =i.NsignalPFGammaCands;
  leadPFChargedHadrCandPt[j] = i.leadPFChargedHadrCandPt;
  byLooseIsolation[j] = i.byLooseIsolation;
  byMediumIsolation[j] = i.byMediumIsolation;
  byTightIsolation[j] = i.byTightIsolation;

  byLooseCombinedIsolationDeltaBetaCorr3Hits[j] = i.byLooseCombinedIsolationDeltaBetaCorr3Hits;
  byMediumCombinedIsolationDeltaBetaCorr3Hits[j] = i.byMediumCombinedIsolationDeltaBetaCorr3Hits;
  byTightCombinedIsolationDeltaBetaCorr3Hits[j] = i.byTightCombinedIsolationDeltaBetaCorr3Hits;
  againstElectronMVA3raw[j] = i.againstElectronMVA3raw;
  againstElectronMVA3category[j] = i.againstElectronMVA3category;
  againstElectronLooseMVA3[j] = i.againstElectronLooseMVA3;
  againstElectronMediumMVA3[j] = i.againstElectronMediumMVA3;
  againstElectronTightMVA3[j] = i.againstElectronTightMVA3;
  againstElectronVTightMVA3[j] = i.againstElectronVTightMVA3;
  againstElectronDeadECAL[j] = i.againstElectronDeadECAL;
  byLooseIsolationMVA[j] = i.byLooseIsolationMVA;
  byMediumIsolationMVA[j] = i.byMediumIsolationMVA;
  byTightIsolationMVA[j] = i.byTightIsolationMVA;
  byLooseIsolationMVA2[j] = i.byLooseIsolationMVA2;
  byMediumIsolationMVA2[j] = i.byMediumIsolationMVA2;
  byTightIsolationMVA2[j] = i.byTightIsolationMVA2;
  againstMuonLoose2[j] = i.againstMuonLoose2;
  againstMuonMedium2[j] = i.againstMuonMedium2;
  againstMuonTight2[j] = i.againstMuonTight2;

}


typedef struct 
{
  void reset() {
    et = -100 ;
    sumet = -100 ;
    sig = -100 ;
    phi = -100 ;
  }
  float et; 
  float sumet;   
  float sig;
  float phi;
}
METInfo;


typedef struct 
{
  float  et[nMetUnc]; float phi[nMetUnc]; float sumet[nMetUnc]; 
  template <class Input> void set(const Input & i, int j)
  {
    et[j]=i.p4.Pt(); 
    phi[j]=i.p4.Phi();
    sumet[j]=i.sumEt;
  }      
}
METUncInfo ;

  
typedef struct 
{
  float mht;
  float ht;  
  float sig;
  float phi;
}
MHTInfo;


typedef struct 
{
  float mass;
  float pt;
  float wMass;
}
TopInfo;


typedef struct 
{
  int run;
  int lumi;
  int event;
  int json;
}
EventInfo;


BTagShapeInterface * nominalShape=0;
BTagShapeInterface * downBCShape=0;
BTagShapeInterface * upBCShape=0;
BTagShapeInterface * downLShape=0;
BTagShapeInterface * upLShape=0; 
BTagShapeInterface * nominalShape4p=0;
BTagShapeInterface * downBCShape4p=0;
BTagShapeInterface * upBCShape4p=0;
BTagShapeInterface * nominalShape1Bin=0;


//------- JetInfo Object -------//
// Contains info pertaining to a set of jets.
//
typedef struct 
{
  // set() - Copy info from SimpleJet j onto JetInfo index i
  void set(const VHbbEvent::SimpleJet & j, int i, const TLorentzVector& jetWithJEC, const TLorentzVector& jetBestMC, int bestMCidIn) 
  {
    pt[i]=j.p4.Pt();
    eta[i]=j.p4.Eta();
    phi[i]=j.p4.Phi();
    e[i]=j.p4.E();
    pt_withJEC[i]=jetWithJEC.Pt();
    eta_withJEC[i]=jetWithJEC.Eta();
    phi_withJEC[i]=jetWithJEC.Phi();
    e_withJEC[i]=jetWithJEC.E();
    pt_bestMC[i]=jetBestMC.Pt();
    if (jetBestMC.P() > 0.001) {
      eta_bestMC[i]=jetBestMC.Eta();
      phi_bestMC[i]=jetBestMC.Phi();
    }
    else {
      eta_bestMC[i]= -99 ;
      phi_bestMC[i]= -99 ;
    }
    e_bestMC[i]=jetBestMC.E();
    bestMCid[i]=bestMCidIn ;
    csv[i]=j.csv;
    jp[i]=j.jp;
    jpb[i]=j.jpb;    
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
    vtxProb[i] = j.vtxProb;
    ssvhe[i] = j.ssvhe  ;
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

    vtxPosition_x[i]=j.vtxPosition.x();
    vtxPosition_y[i]=j.vtxPosition.y();
    vtxPosition_z[i]=j.vtxPosition.z();
    tche[i]=j.tche;
    tchp[i]=j.tchp;
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
      pt[i]=-99; eta[i]=-99; phi[i]=-99;e[i]=-99;
      pt_withJEC[i]=-99; eta_withJEC[i]=-99; phi_withJEC[i]=-99;e_withJEC[i]=-99;
      pt_bestMC[i]=-99; eta_bestMC[i]=-99; phi_bestMC[i]=-99;e_bestMC[i]=-99;
      bestMCid[i]=-99 ;
      jp[i]=-99;jpb[i]=-99;csv[i]=-99;csv_nominal[i]=-99.;csv_upBC[i]=-99.;csv_downBC[i]=-99.;csv_upL[i]=-99.;csv_downL[i]=-99.;csv_nominal1Bin[i]=-99.;
      csvivf[i]=-99; cmva[i]=-99;
      cosTheta[i]=-99; numTracksSV[i]=-99; chf[i]=-99; nhf[i]=-99; cef[i]=-99; nef[i]=-99; nch[i]=-99; nconstituents[i]=-99; flavour[i]=-99; isSemiLeptMCtruth[i]=-99; isSemiLept[i]=-99;      
      SoftLeptpdgId[i] = -99; SoftLeptIdlooseMu[i] = -99;  SoftLeptId95[i] =  -99;   SoftLeptPt[i] = -99;  SoftLeptdR[i] = -99;   SoftLeptptRel[i] = -99; SoftLeptRelCombIso[i] = -99;  
      genPt[i]=-99; genEta[i]=-99; genPhi[i]=-99; JECUnc[i]=-99; ptRaw[i]=-99.; ptLeadTrack[i]=-99.; puJetIdL[i]=-99; puJetIdM[i]=-99; puJetIdT[i]=-99; puJetIdMva[i]=-99; charge[i]=-99; jetArea[i]=-99; selectedTauDR[i] = -99.;
      vtxPosition_x[i]=-99;
      vtxPosition_y[i]=-99;
      vtxPosition_z[i]=-99;
      tche[i]=-99;
      tchp[i]=-99;
      vtxProb[i]=-99;
      ssvhe  [i]=-99;
	}
  }

  float pt[MAXJ];
  float eta[MAXJ];
  float phi[MAXJ];
  float e[MAXJ];
  float pt_withJEC[MAXJ] ;
  float eta_withJEC[MAXJ];
  float phi_withJEC[MAXJ];
  float e_withJEC[MAXJ];
  float pt_bestMC[MAXJ] ;
  float eta_bestMC[MAXJ];
  float phi_bestMC[MAXJ];
  float e_bestMC[MAXJ];
  int   bestMCid[MAXJ];
  float jp[MAXJ];
  float jpb[MAXJ];
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
  float vtxProb[MAXJ];
  float ssvhe  [MAXJ];
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
  float vtxPosition_x[MAXJ];
  float vtxPosition_y[MAXJ];
  float vtxPosition_z[MAXJ];
  float tche[MAXJ];
  float tchp[MAXJ];
}
JetInfo;

//===NOTE : only run on new PAT with genParticles and genJet (Duong 10-20-2015)===
typedef struct 
{
  void reset() {
    for (int i = 0; i < MAXGENOBJ; i++) {
    	pt[i] = -100 ;
        eta[i] = -100 ;
        phi[i] = -100 ;
        charge[i] = -100 ;
        flavor[i] = -100 ;
        status[i] = -100 ;
    }
  }
  void set(reco::GenParticleRefVector::const_iterator& itr, int i) {
    pt[i] = (*itr)->pt() ;
    eta[i] = (*itr)->eta() ;
    phi[i] = (*itr)->phi() ;
    charge[i] = (*itr)->charge() ;
    flavor[i] = (*itr)->pdgId() ;
    status[i] = (*itr)->status() ;
  }
  //void set(reco::JetFlavourMatching& jet, unsigned i) {
  void set(reco::JetFlavourMatchingCollection::const_iterator& jet, unsigned i) {
    edm::RefToBase<reco::Jet> aJet  = (*jet).first ;
    const reco::JetFlavour aFlav = (*jet).second ;
    pt[i] = aJet.get()->pt() ;
    eta[i] = aJet.get()->eta() ;
    phi[i] = aJet.get()->phi() ;
    charge[i] = aJet.get()->charge() ;
    flavor[i] = aFlav.getFlavour() ;
    status[i] = -100 ;
  }
  float pt[MAXGENOBJ] ;
  float eta[MAXGENOBJ] ;
  float phi[MAXGENOBJ] ;
  int charge[MAXGENOBJ] ;
  int flavor[MAXGENOBJ] ;
  int status[MAXGENOBJ] ;
}
GenObjectInfo ; //used to store two leptons from Z, FSR photons and jets

int FindStatus1Lepton(reco::GenParticleRefVector::const_iterator lepIn, reco::GenParticleRefVector::const_iterator& lepOut, vector<reco::GenParticleRefVector::const_iterator>& phoItrs) {
  int status = (*lepIn)->status() ;
  int pdgId = (*lepIn)->pdgId() ;
  if (status == 1) {
    lepOut = lepIn ;
    return 1 ; //found 
  }
  else if (status == 2) {
    //===loop over daughter to find photon and status 1 lepton===
    //+ if status 1 lepton with the same flavor found return the leptOut==
    //+ if status 1 lepton with different flavor found for tau case return the leptOut==
    //+ if status 2 with same flavor found feed it again to FindStatus1Lepton===
    //+ other cases return 0
    const reco::GenParticleRefVector& daughterRefs = (*lepIn)->daughterRefVector();
    bool nSta1(0) ;
    bool nSta2(0) ;
    bool nHadTau(0) ;
    reco::GenParticleRefVector::const_iterator itrTmp ;
    for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
      int pdgIdTmp = (*idr)->pdgId() ;
      int statusTmp = (*idr)->status() ;
      //fill photons
      if (pdgIdTmp == 22 && statusTmp == 1) phoItrs.push_back(idr) ;
      //fill electron and muon, not tau case
      if (abs(pdgId) != 15 && (pdgIdTmp == pdgId) && statusTmp == 1) {
        lepOut = idr ;
        nSta1 += 1 ;
      }
      //fill electron and muon, tau case
      if (abs(pdgId) == 15) {
        if ((abs(pdgIdTmp) == 11 || abs(pdgIdTmp) == 13) && statusTmp == 1) {
          lepOut = idr ;
          nSta1 += 1 ;
        }
        if (abs(pdgIdTmp) > 18) nHadTau += 1 ; //hadronic tau decay found
      } 
      if (pdgIdTmp == pdgId && statusTmp == 2) { //status 2 with same flavor found
        itrTmp = idr ;
        nSta2 +=1 ;
      }
    } //end loop over daughter
    if (nHadTau == 0 && ((nSta1 == 0 && nSta2 == 0) || (nSta1 == 1 && nSta2 == 1) || nSta1 > 1 || nSta2 > 1)) { 
      cout << "\n Warning: found " << nSta1 << " status  1 particles and " << nSta2 << " status 2 particles. Odd case, exit " << pdgId << "  " << status ; 
      for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
        int pdgIdTmp = (*idr)->pdgId() ;
        int statusTmp = (*idr)->status() ;
        cout << "\n Daughter " << pdgIdTmp << "  " << statusTmp ;
      }
      return 0 ; 
    }
    if (nSta1 == 1) return 1 ; //==note lepOut and photon already set inside particle loop
    if (nSta2 == 1) {
      return FindStatus1Lepton(itrTmp, lepOut, phoItrs) ;
    }
  } //end else if (status ==2
  else {
    cout << "\n Warning: It is not status 1 or status 2 particle, exit" ;
    return 0 ;
  }
  return 0 ;
}

///------------THE MAIN FUNCTION!!!--------------///
int main(int argc, char* argv[]) 
{
//  gROOT->Reset();

  TTree *_outTree;
  IVFInfo IVF;
  SimBHadronInfo SimBs;
  float rho,rho25,rhoN;
  int nPVs;
  int nPV1s; //nPV outside of VH cand
  int nTruePUs ;  
  METInfo MET; //pfmetType1corr
  METInfo fakeMET;
  METInfo METnoPU;
  METInfo METnoPUCh;
  METInfo METtype1corr;
  METInfo METtype1p2corr;
  METInfo METnoPUtype1corr;
  METInfo METnoPUtype1p2corr;
  METInfo METtype1diff;
  MHTInfo MHT;
  METUncInfo  metUnc;
  METInfo aMET ; //MET outside VH candidate pfmetType1corr
  METInfo apfMET ; //pfMET outside VH candidate same as Andrew
  TopInfo top;
  EventInfo EVENT;
  // JetInfo jet1,jet2, addJet1, addJet2;
  // lepton1,lepton2;


  JetInfo hJets, aJets, fathFilterJets, aJetsFat;
  LeptonInfo vLeptons, aLeptons, vLeptonsTaus;
  int naJets=0, nhJets=0, nfathFilterJets=0, naJetsFat=0;

  //------- INFO ADDED 2014-02-19 ---------//
  JetInfo    allJets;
  LeptonInfo allMuons, allElectrons;
  int        nallJets = 0;
  int        nallMuons = 0, nallElectrons = 0;
  METInfo    pfMET;
  
  //====NOTE : some gen information stuffs (Duong 10-21-2015) ====
  int zdecayMode(0) ;
  GenObjectInfo genLeps, genPhos ; //==NOTE lep and photon status = 1
  int nGenLep(0), nGenPho(0) ; 
  GenObjectInfo genSta3objs ; //===Status 3 leptons (from Z) and partons===
  int nGenSta3obj(0) ;
  GenObjectInfo genAk5Jets ;
  int nGenAk5Jet(0) ;
  GenObjectInfo genPatPFjets ;
  int nGenPatPFjet(0) ;

  HiggsInfo H,SVH,SimBsH;
  FatHiggsInfo FatH;
  genParticleInfo genZ, genZstar, genWstar, genW,  genH, genB, genBbar; //add here the fatjet higgs
  genTopInfo genTop, genTbar;
  TrackInfo V;
  TrackInfo VTau;
  int nvlep=0,nalep=0,nvlepTau=0;
  float lheV_pt=0; 			  //for the Madgraph sample stitching
  float lheHT=0; 			  //for the Madgraph sample stitching
  float lheNj=0; 			  //for the Madgraph sample stitching
  TrackSharingInfo TkSharing; // track sharing info;
  unsigned int nPdf=0; 		  // number of pdf sets stored
  float PDFweight[MAXPDF];    // for pdf reweighting (only madgraph)

  float HVdPhi,HVMass,HMETdPhi,VMt,deltaPullAngle,deltaPullAngleAK7,deltaPullAngle2,deltaPullAngle2AK7,gendrcc,gendrbb, genZpt, genWpt, genHpt, weightTrig, weightTrigMay,weightTrigV4, weightTrigMET, weightTrigOrMu30, minDeltaPhijetMET,  jetPt_minDeltaPhijetMET , PUweight, PUweightP,PUweightM, PUweightAB, PUweight2011B,PUweight1DObs;
  float PU0,PUp1,PUm1,weightMCProd;
  bool isBadHcalEvent=false;
  float weightEleRecoAndId,weightEleTrigJetMETPart, weightEleTrigElePart,weightEleTrigEleAugPart;
  float weightTrigMET80, weightTrigMET100, weightTrig2CJet20, weightTrigMET150, weightTrigMET802CJet, weightTrigMET1002CJet, weightTrigMETLP;

  float weightTrig2012A, weightTrig2012ADiMuon, weightTrig2012ADiEle, weightTrig2012ASingleMuon, weightTrig2012AMuonPlusWCandPt, weightTrig2012ASingleEle;  
  float weightTrig2012AB, weightTrig2012ABDiMuon, weightTrig2012ABDiEle, weightTrig2012ABSingleMuon, weightTrig2012ABMuonPlusWCandPt, weightTrig2012ABSingleEle;  
  float weightTrig2012, weightTrig2012DiMuon, weightTrig2012DiEle, weightTrig2012SingleMuon, weightTrig2012MuonPlusWCandPt, weightTrig2012SingleEle;  

  float weightTrig2012DiJet30MHT80,weightTrig2012PFMET150,weightTrig2012SumpT100MET100;
  float weightTrig2012APFMET150orDiJetMET, weightTrig2012BPFMET150orDiJetMET, weightTrig2012CPFMET150orDiJetMET; 

  float weightSignalEWK, weightSignalQCD;
 
  int tauPlusMode, tauMinusMode;
  int WplusMode,WminusMode;
  int Vtype,VtypeWithTau,nSvs=0;
  int nSimBs=0;
  int numJets,numBJets,eventFlav;
  // bool isMET80_CJ80, ispfMHT150, isMET80_2CJ20,isMET65_2CJ20, isJETID,isIsoMu17;
  bool triggerFlags[500],hbhe,ecalFlag,totalKinematics,  cschaloFlag,  hcallaserFlag,   trackingfailureFlag , eebadscFlag ;

  float btag1T2CSF=1.,btag2TSF=1.,btag1TSF=1.,btagA0CSF=1., btagA0TSF=1., btag2CSF=1., btag1TA1C=1.;


  // ----------------------------------------------------------------------
  // First Part: 
  //
  //  * enable the AutoLibraryLoader 
  //  * book the histograms of interest 
  //  * open the input file
  // ----------------------------------------------------------------------

  // load framework libraries
  gSystem->Load("libFWCoreFWLite");
  gSystem->Load("libDataFormatsFWLite");
  AutoLibraryLoader::enable();

  //init the PDF set. Need to check for each sample beacuse pdfset and pdfmember info missing in AODSIM.
  //here the list of supported set http://lhapdf.hepforge.org/manual 
  //you can aslo call it by name
  //initPDFSet(int nset,const std::string &filename, int member=0) // If you want to use more than one (but not more than 3) at the same time, initialize them giving nset.

  
//-------------------------------------------------
//  PDF LOADING
//-------------------------------------------------

  //Initialise member in PDF set file filename. If filename contains a "/" character, it will be used as a path, otherwise it will be assumed to be a PDF file in the LHAPDF PDFsets directory.
  //LHAPDF::initPDFSet(1,_pdfset, _pdfmember);
  // (for Madgraph always (1,10042,0) :  mctqe6l) // CTEQ6l (LO fit/NLO alphas)
  LHAPDF::initPDFSet(1, 10042, 0); // this should always be the first
  
  std::vector<std::string> pdfNames;
  //from pdf prescription from LHC PDF group http://arxiv.org/pdf/1101.0538v1.pdf
  //pdfNames.push_back("NNPDF20_100.LHgrid"); //NNPDF20 (Neutral Network 100/1000 sets)
  //pdfNames.push_back("cteq6mE.LHgrid"); //
  pdfNames.push_back("cteq66.LHgrid"); //
  //pdfNames.push_back("MSTW2008nlo68cl.LHgrid"); // MSTW2008 (NLO central)
  //pdfNames.push_back("MSTW2008lo68cl.LHgrid"); // MSTW2008 (LO central)
  //Other possible examples:
  //pdfNames.push_back("cteq6lg.LHgrid"); // CTEQ61 (central value)
  //pdfNames.push_back("cteq61.LHgrid"); // CTEQ61 (central value)
  //pdfNames.push_back("CT10.LHgrid"); // CT10 (central value) 
  //pdfNames.push_back("CT10nlo.LHgrid"); // CT10nlo (new fit - central value)
  //pdfNames.push_back("CT10nnlo.LHgrid"); // CT10nnlo (new fit - central value)
  //pdfNames.push_back("cteq5m1.LHgrid"); //CTEQ5m1 (updated CTEQ5m)
  //pdfNames.push_back("MRST2001lo.LHgrid.gz");// MRST2001lo (LO fit)
  //pdfNames.push_back("MRST2001nlo.LHgrid"); // MRST2001nlo (Standard MSbar)
  //pdfNames.push_back("MRST2001nnlo.LHgrid.gz"); // MRST2001nnlo (NNLO fit)
  //pdfNames.push_back("NNPDF21_nnlo_as_0114.LHgrid"); // NNPDF23_nnlo_as_ (Neutral Network NNLO 0+100 replicas)
  //pdfNames.push_back("NNPDF21_100.LHgrid"); //NNPDF23_nlo_as_ (Neutral Network NLO 0+100 replicas)
  //pdfNames.push_back("MSTW2008lo68cl.LHgrid"); // MSTW2008 (LO central)
  //pdfNames.push_back("MSTW2008nnlo68cl.LHgrid"); // MSTW2008 (NNLO central)
  
  //initializing all the pdf sets
  nPdf=1; // one count the nominal PDF
  for (unsigned int setpdf=2; setpdf <= pdfNames.size()+1; setpdf++)
  {
    std::cout << "pdf set " << setpdf << std::endl;
    //LHAPDF::initPDFSet(setpdf, pdfNames[setpdf-2], 0);  // the member should always be zero (if you do not know what you are doing)
    LHAPDF::initPDFSet(setpdf, pdfNames[setpdf-2]); 	  // the member should always be zero (if you do not know what you are doing)
    //std::cout << "Initialized setpdf " << setpdf << ": " << pdfNames[setpdf-2] << std::endl;
    nPdf+=LHAPDF::numberPDF(setpdf); 					  // count all the members of all the pdfsets
  }
  std::cout << "pdf memebers " << nPdf << std::endl;




//-------------------------------------------------
//  LOAD CONFIGURATION
//-------------------------------------------------

  std::vector<VHbbCandidate> * candZlocal = new std::vector<VHbbCandidate>;
  std::vector<VHbbCandidate> * candWlocal = new std::vector<VHbbCandidate>;

  // parse arguments
  if ( argc < 2 ) return 0;

  // get the python configuration
  PythonProcessDesc builder(argv[1]);
  const edm::ParameterSet& in  = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteInput" );
  const edm::ParameterSet& out = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("fwliteOutput");
  const edm::ParameterSet& ana = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer");

  std::vector<edm::LuminosityBlockRange> jsonVector;
  if ( in.exists("lumisToProcess") ) 
  {
    std::vector<edm::LuminosityBlockRange> const & lumisTemp =
      in.getUntrackedParameter<std::vector<edm::LuminosityBlockRange> > ("lumisToProcess");
    jsonVector.resize( lumisTemp.size() );
    copy( lumisTemp.begin(), lumisTemp.end(), jsonVector.begin() );
  }
  
  // now get each parameter
  int maxEvents_( in.getParameter<int>("maxEvents") );
  int skipEvents_( in.getParameter<int>("skipEvents") );
  int runMin_( in.getParameter<int>("runMin") );
  int runMax_( in.getParameter<int>("runMax") );
  //unsigned int outputEvery_( in.getParameter<unsigned int>("outputEvery") );
  std::string outputFile_( out.getParameter<std::string>("fileName" ) );
  std::vector<std::string> triggers( ana.getParameter<std::vector<std::string> >("triggers") );
  double btagThr =  ana.getParameter<double>("bJetCountThreshold" );
  bool fromCandidate = ana.getParameter<bool>("readFromCandidates");
  bool useHighestPtHiggsZ = ana.getParameter<bool>("useHighestPtHiggsZ");
  bool useHighestPtHiggsW = ana.getParameter<bool>("useHighestPtHiggsW");
  HbbCandidateFinderAlgo * algoZ = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdZ"),useHighestPtHiggsZ);
  HbbCandidateFinderAlgo * algoW = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), ana.getParameter<double>("jetPtThresholdW"),useHighestPtHiggsW );
  //HbbCandidateFinderAlgo * algoRecoverLowPt = new HbbCandidateFinderAlgo(ana.getParameter<bool>("verbose"), 15, true);

  TriggerWeight triggerWeight(ana);

  const edm::ParameterSet& anaAB = builder.processDesc()->getProcessPSet()->getParameter<edm::ParameterSet>("Analyzer2012ABOnly");
  TriggerWeight triggerWeightAB(anaAB);

  BTagWeight btag(2); // 2 operating points "Custom" = 0.5 and "Tight = 0.898"
  BTagSampleEfficiency btagEff( ana.getParameter<std::string>("btagEffFileName" ).c_str() ); 

  std::vector<std::string> inputFiles_( in.getParameter<std::vector<std::string> >("fileNames") );
  //  std::string inputFile( in.getParameter<std::string> ("fileName") );

  std::string PUmcfileName_ = in.getParameter<std::string> ("PUmcfileName") ;
  std::string PUmcfileName2011B_ = in.getParameter<std::string> ("PUmcfileName2011B") ;

  std::string PUdatafileNameAB = in.getParameter<std::string> ("PUdatafileNameAB") ;
  std::string PUdatafileName_ = in.getParameter<std::string> ("PUdatafileName") ;
  std::string PUdatafileName_plusOneSigma = in.getParameter<std::string> ("PUdatafileNamePlus") ;
  std::string PUdatafileName_minusOneSigma = in.getParameter<std::string> ("PUdatafileNameMinus") ;
  std::string PUdatafileName2011B_ = in.getParameter<std::string> ("PUdatafileName2011B") ;
  std::string Weight3DfileName_ = in.getParameter<std::string> ("Weight3DfileName") ;
  
  JECFWLite jec(ana.getParameter<std::string>("jecFolder"));

  bool isMC_( ana.getParameter<bool>("isMC") );  
  TriggerReader trigger(false);
  TriggerReader patFilters(false);

  vector<RLE> badEvents;
  if(!isMC_)
  {
    readBadEvents(in.getParameter<std::string> ("badEventsFileName").c_str(), badEvents);
  }
  if(isMC_) 
  {
    nominalShape = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 0.0); 
    upBCShape = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 1.0, 0.0); 
    downBCShape = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), -1.0, 0.0); 
    upLShape = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 1.0); 
    downLShape = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, -1.0); 
    nominalShape4p = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 0.0,true,1.003,1.003); 
    upBCShape4p = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 0.0,true,1.001,1.001); 
    downBCShape4p = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 0.0,true,1.005,1.005); 
    nominalShape1Bin = new BTagShapeInterface(ana.getParameter<std::string>("csvDiscr").c_str(), 0.0, 0.0,true,1.001,1.001,1); 
  }
  edm::LumiReWeighting   lumiWeights,lumiWeightsPl,lumiWeightsMi,lumiWeightsAB;
  edm::LumiReWeighting   lumiWeights1DObs;
  //edm::Lumi3DReWeighting   lumiWeights2011B;
  if(isMC_)
  {
    lumiWeights = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_ , "pileup", "pileup");
    lumiWeightsPl = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_plusOneSigma , "pileup", "pileup");
    lumiWeightsMi = edm::LumiReWeighting(PUmcfileName_,PUdatafileName_minusOneSigma , "pileup", "pileup");
    lumiWeightsAB = edm::LumiReWeighting(PUmcfileName_,PUdatafileNameAB , "pileup", "pileup");
    lumiWeights1DObs = edm::LumiReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");

    //lumiWeights2011B = edm::Lumi3DReWeighting(PUmcfileName2011B_,PUdatafileName2011B_ , "pileup", "pileup");
    //if(Weight3DfileName_!="")
    //{
    //	lumiWeights2011B.weight3D_init(Weight3DfileName_.c_str());
    //}
    //else
    //{
    //	lumiWeights2011B.weight3D_init(1.0); // generate the weights the fisrt time;
    //}
  }

  bool doFillMoreMCtruth = ana.getParameter<bool> ("doFillMoreMCtruth") ; 

  //
  // CREATE OUTPUT FILE, ADD A TREE, ADD BRANCHES TO TREE.
  // 

  //TFile *_outPUFile	= new TFile((outputFile_+"_PU").c_str(), "recreate");	
  TFile *_outFile = new TFile(outputFile_.c_str(), "recreate");	
  TH1F * count = new TH1F("Count","Count", 1,0,2 );
  TH1F * countWithPU = new TH1F("CountWithPU","CountWithPU", 1,0,2 );
  TH1F * countWithPUP = new TH1F("CountWithPUP","CountWithPU plus one sigma", 1,0,2 );
  TH1F * countWithPUM = new TH1F("CountWithPUM","CountWithPU minus one sigma", 1,0,2 );
  TH1F * countWithPUAB = new TH1F("CountWithPUAB","CountWithPU 2012AB", 1,0,2 );

  TH1F * countWithPU2011B = new TH1F("CountWithPU2011B","CountWithPU2011B", 1,0,2 );
  TH1F * coutnWithMCProd = new TH1F("CountWithMCProd","CountWithMCProd", 1,0,2 );
  TH1F * countWithPUMCProd = new TH1F("CountWithPUMCProd","CountWithPUMCProd", 1,0,2 );
  TH1F * countWithSignalQCDcorrections = new TH1F("countWithSignalQCDcorrections","countWithSignalQCDcorrections", 1,0,2 );

  TH3F * input3DPU = new TH3F("Input3DPU","Input3DPU", 36,-0.5,35.5,36,-0.5,35.5, 36,-0.5,35.5 );
  
  TH1F * pu = new TH1F("pileup","",51,-0.5,50.5);
  
  TH1F* hEvent = new TH1F("Event", "", 3, 0, 3) ; 

  _outTree = new TTree("tree", "myTree");
  
  _outTree->Branch("H",  	&H,  	   "HiggsFlag/I:mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");
  _outTree->Branch("V",		&V,  	   "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("VTau",  	&VTau,     "mass/F:pt/F:eta:phi/F");
  _outTree->Branch("FatH",  	&FatH,     "FatHiggsFlag/I:mass/F:pt/F:eta:phi/F:filteredmass/F:filteredpt/F:filteredeta/F:filteredphi/F");
  _outTree->Branch("lheV_pt",  	&lheV_pt,  "lheV_pt/F");
  _outTree->Branch("lheHT",  	&lheHT,    "lheHT/F");
  _outTree->Branch("lheNj",  	&lheNj,    "lheNj/F");
  _outTree->Branch("nPdf",  	&nPdf,	   "nPdf/I");
  _outTree->Branch("PDFweight", PDFweight, "PDFweight[nPdf]/F");

  _outTree->Branch("genZ",     &genZ,     "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genZstar", &genZstar, "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genW",     &genW,     "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genWstar", &genWstar, "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genH",     &genH,     "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genB",     &genB,     "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genBbar",  &genBbar,  "mass/F:pt/F:eta:phi/F:status/F:charge:momid/F");
  _outTree->Branch("genTop",   &genTop,   "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");
  _outTree->Branch("genTbar",  &genTbar,  "bmass/F:bpt/F:beta:bphi/F:bstatus/F:wdau1mass/F:wdau1pt/F:wdau1eta:wdau1phi/F:wdau1id/F:wdau2mass/F:wdau2pt/F:wdau2eta:wdau2phi/F:wdau2id/F");

  _outTree->Branch("TkSharing", &TkSharing, "HiggsCSVtkSharing/b:HiggsIPtkSharing:HiggsSVtkSharing:FatHiggsCSVtkSharing:FatHiggsIPtkSharing:FatHiggsSVtkSharing");


  _outTree->Branch("weightMCProd"             ,  &weightMCProd                 ,   "weightMCProd/F");
  _outTree->Branch("isBadHcalEvent"             ,  &isBadHcalEvent                 ,   "isBadHcalEvent/b");

  _outTree->Branch("nhJets",	      &nhJets,		"nhJets/I");
  _outTree->Branch("nfathFilterJets", &nfathFilterJets, "nfathFilterJets/I");

//===NOTE : don't save hJet and fathFilterJet=== Duong (10-9-2015)
/*  _outTree->Branch("hJet_pt",hJets.pt ,"pt[nhJets]/F");
  _outTree->Branch("hJet_eta",hJets.eta ,"eta[nhJets]/F");
  _outTree->Branch("hJet_phi",hJets.phi ,"phi[nhJets]/F");
  _outTree->Branch("hJet_e",hJets.e ,"e[nhJets]/F");
  _outTree->Branch("hJet_csv",hJets.csv ,"csv[nhJets]/F");
  _outTree->Branch("hJet_csv_nominal",hJets.csv_nominal ,"csv_nominal[nhJets]/F");
  _outTree->Branch("hJet_csv_upBC",hJets.csv_upBC ,"csv_upBC[nhJets]/F");
  _outTree->Branch("hJet_csv_downBC",hJets.csv_downBC ,"csv_downBC[nhJets]/F");
  _outTree->Branch("hJet_csv_upL",hJets.csv_upL ,"csv_upL[nhJets]/F");
  _outTree->Branch("hJet_csv_downL",hJets.csv_downL ,"csv_downL[nhJets]/F");
  _outTree->Branch("hJet_csv_nominal4p",hJets.csv_nominal4p ,"csv_nominal4p[nhJets]/F");
  _outTree->Branch("hJet_csv_upBC4p",hJets.csv_upBC4p ,"csv_upBC4p[nhJets]/F");
  _outTree->Branch("hJet_csv_downBC4p",hJets.csv_downBC4p ,"csv_downBC4p[nhJets]/F");
  _outTree->Branch("hJet_csv_nominal1Bin",hJets.csv_nominal1Bin ,"csv_nominal1Bin[nhJets]/F");
  _outTree->Branch("hJet_csvivf",hJets.csvivf ,"csvivf[nhJets]/F");
  _outTree->Branch("hJet_cmva",hJets.cmva ,"cmva[nhJets]/F");
  _outTree->Branch("hJet_cosTheta",hJets.cosTheta ,"cosTheta[nhJets]/F");
  _outTree->Branch("hJet_numTracksSV",hJets.numTracksSV ,"numTracksSV[nhJets]/I");
  _outTree->Branch("hJet_chf",hJets.chf ,"chf[nhJets]/F");
  _outTree->Branch("hJet_nhf",hJets.nhf ,"nhf[nhJets]/F");
  _outTree->Branch("hJet_cef",hJets.cef ,"cef[nhJets]/F");
  _outTree->Branch("hJet_nef",hJets.nef ,"nef[nhJets]/F");
  _outTree->Branch("hJet_nch",hJets.nch ,"nch[nhJets]/F");
  _outTree->Branch("hJet_nconstituents",hJets.nconstituents ,"nconstituents[nhJets]");
  _outTree->Branch("hJet_flavour",hJets.flavour ,"flavour[nhJets]/F");
  _outTree->Branch("hJet_isSemiLept",hJets.isSemiLept ,"isSemiLept[nhJets]/I");
  _outTree->Branch("hJet_isSemiLeptMCtruth",hJets.isSemiLeptMCtruth ,"isSemiLeptMCtruth[nhJets]/I");
  _outTree->Branch("hJet_SoftLeptpdgId", hJets.SoftLeptpdgId , "SoftLeptpdgId[nhJets]/I");
  _outTree->Branch("hJet_SoftLeptIdlooseMu", hJets.SoftLeptIdlooseMu , "SoftLeptIdlooseMu[nhJets]/I");
  _outTree->Branch("hJet_SoftLeptId95", hJets.SoftLeptId95 , "SoftLeptId95[nhJets]/I");
  _outTree->Branch("hJet_SoftLeptPt", hJets.SoftLeptPt , "SoftLeptPt[nhJets]/F");
  _outTree->Branch("hJet_SoftLeptdR", hJets.SoftLeptdR , "SoftLeptdR[nhJets]/F");
  _outTree->Branch("hJet_SoftLeptptRel", hJets.SoftLeptptRel , "SoftLeptptRel[nhJets]/F");
  _outTree->Branch("hJet_SoftLeptRelCombIso", hJets.SoftLeptRelCombIso , "SoftLeptRelCombIso[nhJets]/F");
  _outTree->Branch("hJet_genPt",hJets.genPt ,"genPt[nhJets]/F");
  _outTree->Branch("hJet_genEta",hJets.genEta ,"genEta[nhJets]/F");
  _outTree->Branch("hJet_genPhi",hJets.genPhi ,"genPhi[nhJets]/F");
  _outTree->Branch("hJet_JECUnc",hJets.JECUnc ,"JECUnc[nhJets]/F");
  _outTree->Branch("hJet_vtxMass",hJets.vtxMass ,"vtxMass[nhJets]/F");
  _outTree->Branch("hJet_vtxPt",hJets.vtxPt ,"vtxPt[nhJets]/F");
  _outTree->Branch("hJet_vtxEta",hJets.vtxEta ,"vtxEta[nhJets]/F");
  _outTree->Branch("hJet_vtxPhi",hJets.vtxPhi ,"vtxPhi[nhJets]/F");
  _outTree->Branch("hJet_vtxE",hJets.vtxE ,"vtxE[nhJets]/F");
  _outTree->Branch("hJet_vtx3dL",hJets.vtx3dL ,"vtx3dL[nhJets]/F");
  _outTree->Branch("hJet_vtx3deL",hJets.vtx3deL ,"vtx3deL[nhJets]/F");
  _outTree->Branch("hJet_id",hJets.id ,"id[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVL",hJets.SF_CSVL ,"SF_CSVL[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVM",hJets.SF_CSVM ,"SF_CSVM[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVT",hJets.SF_CSVT ,"SF_CSVT[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVLerr",hJets.SF_CSVLerr ,"SF_CSVLerr[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVMerr",hJets.SF_CSVMerr ,"SF_CSVMerr[nhJets]/b");
  _outTree->Branch("hJet_SF_CSVTerr",hJets.SF_CSVTerr ,"SF_CSVTerr[nhJets]/b");
  _outTree->Branch("hJet_ptRaw",hJets.ptRaw ,"ptRaw[nhJets]/F");
  _outTree->Branch("hJet_ptLeadTrack",hJets.ptLeadTrack ,"ptLeadTrack[nhJets]/F");
  _outTree->Branch("hJet_puJetIdL",hJets.puJetIdL ,"puJetIdL[nhJets]/F");
  _outTree->Branch("hJet_puJetIdM",hJets.puJetIdM ,"puJetIdM[nhJets]/F");
  _outTree->Branch("hJet_puJetIdT",hJets.puJetIdT ,"puJetIdT[nhJets]/F");
  _outTree->Branch("hJet_puJetIdMva",hJets.puJetIdMva ,"puJetIdMva[nhJets]/F");
  _outTree->Branch("hJet_charge",hJets.charge ,"charge[nhJets]/F");

  _outTree->Branch("fathFilterJets_pt",fathFilterJets.pt ,"pt[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_eta",fathFilterJets.eta ,"eta[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_phi",fathFilterJets.phi ,"phi[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_e",fathFilterJets.e ,"e[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_csv",fathFilterJets.csv ,"csv[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_chf",fathFilterJets.chf ,"chf[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_ptRaw",fathFilterJets.ptRaw ,"ptRaw[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_ptLeadTrack",fathFilterJets.ptLeadTrack ,"ptLeadTrack[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_flavour",fathFilterJets.flavour ,"flavour[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_isSemiLept",fathFilterJets.isSemiLept ,"isSemiLept[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_isSemiLeptMCtruth",fathFilterJets.isSemiLeptMCtruth ,"isSemiLeptMCtruth[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_genPt",fathFilterJets.genPt ,"genPt[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_genEta",fathFilterJets.genEta ,"genEta[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_genPhi",fathFilterJets.genPhi ,"genPhi[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtxMass",fathFilterJets.vtxMass ,"vtxMass[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtx3dL",fathFilterJets.vtx3dL ,"vtx3dL[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtx3deL",fathFilterJets.vtx3deL ,"vtx3deL[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtxPt",fathFilterJets.vtxPt ,"vtxPt[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtxEta",fathFilterJets.vtxEta ,"vtxEta[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtxPhi",fathFilterJets.vtxPhi ,"vtxPhi[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_vtxE",fathFilterJets.vtxE ,"vtxE[nfathFilterJets]/F");

  //_outTree->Branch("fathFilterJets_pt",fathFilterJets.pt ,"pt[nfathFilterJets]/F");
  //_outTree->Branch("fathFilterJets_eta",fathFilterJets.eta ,"eta[nfathFilterJets]/F");
  //_outTree->Branch("fathFilterJets_phi",fathFilterJets.phi ,"phi[nfathFilterJets]/F");
  //_outTree->Branch("fathFilterJets_e",fathFilterJets.e ,"e[nfathFilterJets]/F");
  //_outTree->Branch("fathFilterJets_csv",fathFilterJets.csv ,"csv[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_csvivf",fathFilterJets.csvivf ,"csvivf[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_cmva",fathFilterJets.cmva ,"cmva[nfathFilterJets]/F");
  //_outTree->Branch("fathFilterJets_flavour",fathFilterJets.flavour ,"flavour[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_cosTheta",fathFilterJets.cosTheta ,"cosTheta[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_jetArea",fathFilterJets.jetArea ,"jetArea[nfathFilterJets]/F");

  //ccla 25Jan13 
  _outTree->Branch("fathFilterJets_nconstituents",fathFilterJets.nconstituents ,"nconstituents[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_JECUnc",fathFilterJets.JECUnc ,"JECUnc[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_SoftLeptPt", fathFilterJets.SoftLeptPt , "SoftLeptPt[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_SoftLeptdR", fathFilterJets.SoftLeptdR , "SoftLeptdR[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_SoftLeptptRel", fathFilterJets.SoftLeptptRel , "SoftLeptptRel[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_csv_nominal",fathFilterJets.csv_nominal ,"csv_nominal[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_nhf",fathFilterJets.nhf ,"nhf[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_cef",fathFilterJets.cef ,"cef[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_nef",fathFilterJets.nef ,"nef[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_nch",fathFilterJets.nch ,"nch[nfathFilterJets]/F");
  _outTree->Branch("fathFilterJets_id",fathFilterJets.id ,"id[nfathFilterJets]/b");
*/
  
  // ------------------ ADDED BRANCHES FOR ALL JET, ALL LEPTON SETS ----------------------- //
  _outTree->Branch("nallJets",  	       	&nallJets,  				"nallJets/I");  
  _outTree->Branch("allJet_pt",				allJets.pt,					"pt[nallJets]/F");
  _outTree->Branch("allJet_eta",			allJets.eta,				"eta[nallJets]/F");
  _outTree->Branch("allJet_phi",			allJets.phi,				"phi[nallJets]/F");
  _outTree->Branch("allJet_e",				allJets.e,           "e[nallJets]/F");
  _outTree->Branch("allJet_pt_withJEC",			allJets.pt_withJEC,  "pt_withJEC[nallJets]/F");
  _outTree->Branch("allJet_eta_withJEC",		allJets.eta_withJEC, "eta_withJEC[nallJets]/F");
  _outTree->Branch("allJet_phi_withJEC",		allJets.phi_withJEC, "phi_withJEC[nallJets]/F");
  _outTree->Branch("allJet_e_withJEC",			allJets.e_withJEC,   "e_withJEC[nallJets]/F");
  _outTree->Branch("allJet_pt_bestMC",			allJets.pt_bestMC,   "pt_bestMC[nallJets]/F");
  _outTree->Branch("allJet_eta_bestMC",			allJets.eta_bestMC,  "eta_bestMC[nallJets]/F");
  _outTree->Branch("allJet_phi_bestMC",			allJets.phi_bestMC,  "phi_bestMC[nallJets]/F");
  _outTree->Branch("allJet_e_bestMC",			allJets.e_bestMC,    "e_bestMC[nallJets]/F");
  _outTree->Branch("allJet_bestMCid",			allJets.bestMCid,    "bestMCid[nallJets]/I");

  _outTree->Branch("allJet_jp",			    allJets.jp, 				"jp[nallJets]/F");
  _outTree->Branch("allJet_jpb",			allJets.jpb, 				"jpb[nallJets]/F");
  _outTree->Branch("allJet_csv",			allJets.csv, 				"csv[nallJets]/F");
  _outTree->Branch("allJet_csv_nominal",	allJets.csv_nominal, 		"csv_nominal[nallJets]/F");
  _outTree->Branch("allJet_csv_upBC",		allJets.csv_upBC, 			"csv_upBC[nallJets]/F");
  _outTree->Branch("allJet_csv_downBC",		allJets.csv_downBC, 		"csv_downBC[nallJets]/F");
  _outTree->Branch("allJet_csv_upL",		allJets.csv_upL, 			"csv_upL[nallJets]/F");
  _outTree->Branch("allJet_csv_downL",		allJets.csv_downL, 			"csv_downL[nallJets]/F");
  _outTree->Branch("allJet_csv_nominal4p",	allJets.csv_nominal4p, 		"csv_nominal4p[nallJets]/F");
  _outTree->Branch("allJet_csv_upBC4p",		allJets.csv_upBC4p, 		"csv_upBC4p[nallJets]/F");
  _outTree->Branch("allJet_csv_downBC4p",	allJets.csv_downBC4p, 		"csv_downBC4p[nallJets]/F");
  _outTree->Branch("allJet_csv_nominal1Bin",allJets.csv_nominal1Bin, 	"csv_nominal1Bin[nallJets]/F");
  _outTree->Branch("allJet_csvivf",			allJets.csvivf, 			"csvivf[nallJets]/F");
  _outTree->Branch("allJet_cmva",			allJets.cmva, 				"cmva[nallJets]/F");
  _outTree->Branch("allJet_cosTheta",		allJets.cosTheta, 			"cosTheta[nallJets]/F");
  _outTree->Branch("allJet_numTracksSV",	allJets.numTracksSV, 		"numTracksSV[nallJets]/I");
  _outTree->Branch("allJet_chf",			allJets.chf, 				"chf[nallJets]/F");
  _outTree->Branch("allJet_nhf",			allJets.nhf, 				"nhf[nallJets]/F");
  _outTree->Branch("allJet_cef",			allJets.cef, 				"cef[nallJets]/F");
  _outTree->Branch("allJet_nef",			allJets.nef, 				"nef[nallJets]/F");
  _outTree->Branch("allJet_nch",			allJets.nch, 				"nch[nallJets]/F");
  _outTree->Branch("allJet_nconstituents",	allJets.nconstituents,		"nconstituents[nallJets]");
  _outTree->Branch("allJet_flavour",		allJets.flavour, 			"flavour[nallJets]/F");
  _outTree->Branch("allJet_isSemiLept",		allJets.isSemiLept, 		"isSemiLept[nallJets]/I");
  _outTree->Branch("allJet_isSemiLeptMCtruth",	allJets.isSemiLeptMCtruth, 	"isSemiLeptMCtruth[nallJets]/I");
  _outTree->Branch("allJet_SoftLeptpdgId",	allJets.SoftLeptpdgId,  	"SoftLeptpdgId[nallJets]/I");
  _outTree->Branch("allJet_SoftLeptIdlooseMu", 	allJets.SoftLeptIdlooseMu,  	"SoftLeptIdlooseMu[nallJets]/I");
  _outTree->Branch("allJet_SoftLeptId95", 	allJets.SoftLeptId95, 		"SoftLeptId95[nallJets]/I");
  _outTree->Branch("allJet_SoftLeptPt", 	allJets.SoftLeptPt, 		"SoftLeptPt[nallJets]/F");
  _outTree->Branch("allJet_SoftLeptdR", 	allJets.SoftLeptdR, 		"SoftLeptdR[nallJets]/F");
  _outTree->Branch("allJet_SoftLeptptRel", 	allJets.SoftLeptptRel, 		"SoftLeptptRel[nallJets]/F");
  _outTree->Branch("allJet_SoftLeptRelCombIso", allJets.SoftLeptRelCombIso,  	"SoftLeptRelCombIso[nallJets]/F");
  _outTree->Branch("allJet_puJetIdL", 		allJets.puJetIdL,  			"puJetIdL[nallJets]/F");
  _outTree->Branch("allJet_puJetIdM", 		allJets.puJetIdM,  			"puJetIdM[nallJets]/F");
  _outTree->Branch("allJet_puJetIdT", 		allJets.puJetIdT,  			"puJetIdT[nallJets]/F");
  _outTree->Branch("allJet_puJetIdMva", 	allJets.puJetIdMva,  		"puJetIdMva[nallJets]/F");
  _outTree->Branch("allJet_charge",			allJets.charge, 			"charge[nallJets]/F");
  _outTree->Branch("allJet_genPt",			allJets.genPt, 				"genPt[nallJets]/F");
  _outTree->Branch("allJet_genEta",			allJets.genEta, 			"genEta[nallJets]/F");
  _outTree->Branch("allJet_genPhi",			allJets.genPhi, 			"genPhi[nallJets]/F");
  _outTree->Branch("allJet_JECUnc",			allJets.JECUnc, 			"JECUnc[nallJets]/F");
  _outTree->Branch("allJet_vtxMass",		allJets.vtxMass, 			"vtxMass[nallJets]/F");
  _outTree->Branch("allJet_vtx3dL",			allJets.vtx3dL, 			"vtx3dL[nallJets]/F");
  _outTree->Branch("allJet_vtx3deL",		allJets.vtx3deL,	 		"vtx3deL[nallJets]/F");
  _outTree->Branch("allJet_vtxProb",        allJets.vtxProb,            "vtxProb[nallJets]/F");
  _outTree->Branch("allJet_ssvhe",          allJets.ssvhe,              "ssvhe[nallJets]/F");
  _outTree->Branch("allJet_id",				allJets.id, 				"id[nallJets]/b");
  _outTree->Branch("allJet_SF_CSVL",		allJets.SF_CSVL, 			"SF_CSVL[nallJets]/F");
  _outTree->Branch("allJet_SF_CSVM",		allJets.SF_CSVM, 			"SF_CSVM[nallJets]/F");
  _outTree->Branch("allJet_SF_CSVT",		allJets.SF_CSVT, 			"SF_CSVT[nallJets]/F");
  _outTree->Branch("allJet_SF_CSVLerr",		allJets.SF_CSVLerr, 		"SF_CSVLerr[nallJets]/F");
  _outTree->Branch("allJet_SF_CSVMerr",		allJets.SF_CSVMerr, 		"SF_CSVMerr[nallJets]/F");
  _outTree->Branch("allJet_SF_CSVTerr",     allJets.SF_CSVTerr,         "SF_CSVTerr[nallJets]/F");
  _outTree->Branch("allJet_tche",           allJets.tche,               "tche[nallJets]/F");
  _outTree->Branch("allJet_tchp",           allJets.tchp,               "tchp[nallJets]/F");
  _outTree->Branch("allJet_vtxPosition_x",  allJets.vtxPosition_x,      "vtxPosition_x[nallJets]/F");
  _outTree->Branch("allJet_vtxPosition_y",  allJets.vtxPosition_y,      "vtxPosition_y[nallJets]/F");
  _outTree->Branch("allJet_vtxPosition_z",  allJets.vtxPosition_z,      "vtxPosition_z[nallJets]/F");


  _outTree->Branch("nallMuons",    				&nallMuons,				"nallMuons/I"    );
  _outTree->Branch("allMuon_mass",				allMuons.mass,			"mass[nallMuons]/F");
  _outTree->Branch("allMuon_pt",				allMuons.pt,			"pt[nallMuons]/F");
  _outTree->Branch("allMuon_eta",				allMuons.eta,			"eta[nallMuons]");
  _outTree->Branch("allMuon_phi",				allMuons.phi,			"phi[nallMuons]/F");
  _outTree->Branch("allMuon_aodCombRelIso",		allMuons.aodCombRelIso,	"aodCombRelIso[nallMuons]/F");
  _outTree->Branch("allMuon_pfCombRelIso",		allMuons.pfCombRelIso,	"pfCombRelIso[nallMuons]/F");
  _outTree->Branch("allMuon_photonIso",			allMuons.photonIso,		"photonIso[nallMuons]/F");
  _outTree->Branch("allMuon_neutralHadIso",		allMuons.neutralHadIso,	"neutralHadIso[nallMuons]/F");
  _outTree->Branch("allMuon_chargedHadIso",		allMuons.chargedHadIso,	"chargedHadIso[nallMuons]/F");
  _outTree->Branch("allMuon_chargedPUIso",		allMuons.chargedPUIso,	"chargedPUIso[nallMuons]/F");
  _outTree->Branch("allMuon_particleIso",		allMuons.particleIso,	"particleIso[nallMuons]/F");
  _outTree->Branch("allMuon_dxy",				allMuons.dxy,			"dxy[nallMuons]/F");
  _outTree->Branch("allMuon_dz",				allMuons.dz,			"dz[nallMuons]/F");
  _outTree->Branch("allMuon_type",				allMuons.type,			"type[nallMuons]/I");
  _outTree->Branch("allMuon_id80",				allMuons.id80,			"id80[nallMuons]/F");
  _outTree->Branch("allMuon_id95",				allMuons.id95,			"id95[nallMuons]/F");
  _outTree->Branch("allMuon_vbtf",				allMuons.vbtf,			"vbtf[nallMuons]/F");
  _outTree->Branch("allMuon_id80NoIso",			allMuons.id80NoIso,		"id80NoIso[nallMuons]/F");
  _outTree->Branch("allMuon_genPt",				allMuons.genPt,			"genPt[nallMuons]/F");
  _outTree->Branch("allMuon_genEta",			allMuons.genEta,		"genEta[nallMuons]");
  _outTree->Branch("allMuon_genPhi",			allMuons.genPhi,		"genPhi[nallMuons]/F");
  _outTree->Branch("allMuon_charge",			allMuons.charge,		"charge[nallMuons]/F");
  _outTree->Branch("allMuon_pfCorrIso",			allMuons.pfCorrIso,		"pfCorrIso[nallMuons]/F");
  _outTree->Branch("allMuon_pfCorrIsoHCP",		allMuons.pfCorrIsoHCP,	"pfCorrIsoHCP[nallMuons]/F");
  _outTree->Branch("allMuon_id2012tight",		allMuons.id2012tight,	"id2012tight[nallMuons]/F");
  _outTree->Branch("allMuon_idMVAnotrig",		allMuons.idMVAnotrig,	"idMVAnotrig[nallMuons]/F");
  _outTree->Branch("allMuon_idMVAtrig",			allMuons.idMVAtrig,		"idMVAtrig[nallMuons]/F");
  _outTree->Branch("allMuon_idMVApresel",		allMuons.idMVApresel,	"idMVApresel[nallMuons]/F");
  _outTree->Branch("allMuon_innerHits",			allMuons.innerHits,		"innerHits[nallMuons]/F");
  _outTree->Branch("allMuon_photonIsoDoubleCount",allMuons.photonIsoDoubleCount,"photonIsoDoubleCount[nallMuons]/F");
  _outTree->Branch("allMuon_wpHWW",				allMuons.wpHWW,			"wpHWW[nallMuons]/F");
  _outTree->Branch("allMuon_wp95",				allMuons.wp95,			"wp95[nallMuons]/F");
  _outTree->Branch("allMuon_wp90",				allMuons.wp90,			"wp90[nallMuons]/F");
  _outTree->Branch("allMuon_wp85",				allMuons.wp85,			"wp85[nallMuons]/F");
  _outTree->Branch("allMuon_wp80",				allMuons.wp80,			"wp80[nallMuons]/F");
  _outTree->Branch("allMuon_wp70",				allMuons.wp70,			"wp70[nallMuons]/F");
  _outTree->Branch("allMuon_tIso",              allMuons.tIso      ,    "tIso[nallMuons]/F"       );
  _outTree->Branch("allMuon_eIso",              allMuons.eIso      ,    "eIso[nallMuons]/F"       );
  _outTree->Branch("allMuon_hIso",              allMuons.hIso      ,    "hIso[nallMuons]/F"       );
  _outTree->Branch("allMuon_pfChaIso",          allMuons.pfChaIso  ,    "pfChaIso[nallMuons]/F"   );
  _outTree->Branch("allMuon_pfChaPUIso",        allMuons.pfChaPUIso,    "pfChaPUIso[nallMuons]/F" );
  _outTree->Branch("allMuon_pfPhoIso",          allMuons.pfPhoIso  ,    "pfPhoIso[nallMuons]/F"   );
  _outTree->Branch("allMuon_pfNeuIso",          allMuons.pfNeuIso  ,    "pfNeuIso[nallMuons]/F"   );

  _outTree->Branch("nallElectrons",    			&nallElectrons,				"nallElectrons/I"    );
  _outTree->Branch("allElectron_mass",			allElectrons.mass,			"mass[nallElectrons]/F");
  _outTree->Branch("allElectron_pt",			allElectrons.pt,			"pt[nallElectrons]/F");
  _outTree->Branch("allElectron_eta",			allElectrons.eta,			"eta[nallElectrons]");
  _outTree->Branch("allElectron_phi",			allElectrons.phi,			"phi[nallElectrons]/F");
  _outTree->Branch("allElectron_aodCombRelIso",	allElectrons.aodCombRelIso,	"aodCombRelIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfCombRelIso",	allElectrons.pfCombRelIso,	"pfCombRelIso[nallElectrons]/F");
  _outTree->Branch("allElectron_photonIso",		allElectrons.photonIso,		"photonIso[nallElectrons]/F");
  _outTree->Branch("allElectron_neutralHadIso",	allElectrons.neutralHadIso,	"neutralHadIso[nallElectrons]/F");
  _outTree->Branch("allElectron_chargedHadIso",	allElectrons.chargedHadIso,	"chargedHadIso[nallElectrons]/F");
  _outTree->Branch("allElectron_chargedPUIso",	allElectrons.chargedPUIso,	"chargedPUIso[nallElectrons]/F");
  _outTree->Branch("allElectron_particleIso",	allElectrons.particleIso,	"particleIso[nallElectrons]/F");
  _outTree->Branch("allElectron_dxy",			allElectrons.dxy,			"dxy[nallElectrons]/F");
  _outTree->Branch("allElectron_dz",			allElectrons.dz,			"dz[nallElectrons]/F");
  _outTree->Branch("allElectron_type",			allElectrons.type,			"type[nallElectrons]/I");
  _outTree->Branch("allElectron_id80",			allElectrons.id80,			"id80[nallElectrons]/F");
  _outTree->Branch("allElectron_id95",			allElectrons.id95,			"id95[nallElectrons]/F");
  _outTree->Branch("allElectron_vbtf",			allElectrons.vbtf,			"vbtf[nallElectrons]/F");
  _outTree->Branch("allElectron_id80NoIso",		allElectrons.id80NoIso,		"id80NoIso[nallElectrons]/F");
  _outTree->Branch("allElectron_genPt",			allElectrons.genPt,			"genPt[nallElectrons]/F");
  _outTree->Branch("allElectron_genEta",		allElectrons.genEta,		"genEta[nallElectrons]");
  _outTree->Branch("allElectron_genPhi",		allElectrons.genPhi,		"genPhi[nallElectrons]/F");
  _outTree->Branch("allElectron_charge",		allElectrons.charge,		"charge[nallElectrons]/F");
  _outTree->Branch("allElectron_pfCorrIso",		allElectrons.pfCorrIso,		"pfCorrIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfCorrIsoHCP",	allElectrons.pfCorrIsoHCP,	"pfCorrIsoHCP[nallElectrons]/F");
  _outTree->Branch("allElectron_id2012tight",	allElectrons.id2012tight,	"id2012tight[nallElectrons]/F");
  _outTree->Branch("allElectron_idMVAnotrig",	allElectrons.idMVAnotrig,	"idMVAnotrig[nallElectrons]/F");
  _outTree->Branch("allElectron_idMVAtrig",		allElectrons.idMVAtrig,		"idMVAtrig[nallElectrons]/F");
  _outTree->Branch("allElectron_idMVApresel",	allElectrons.idMVApresel,	"idMVApresel[nallElectrons]/F");
  _outTree->Branch("allElectron_innerHits",		allElectrons.innerHits,		"innerHits[nallElectrons]/F");
  _outTree->Branch("allElectron_photonIsoDoubleCount",allElectrons.photonIsoDoubleCount,"photonIsoDoubleCount[nallElectrons]/F");
  _outTree->Branch("allElectron_wpHWW",			allElectrons.wpHWW,			"wpHWW[nallElectrons]/F");
  _outTree->Branch("allElectron_wp95",			allElectrons.wp95,			"wp95[nallElectrons]/F");
  _outTree->Branch("allElectron_wp90",			allElectrons.wp90,			"wp90[nallElectrons]/F");
  _outTree->Branch("allElectron_wp85",			allElectrons.wp85,			"wp85[nallElectrons]/F");
  _outTree->Branch("allElectron_wp80",			allElectrons.wp80,			"wp80[nallElectrons]/F");
  _outTree->Branch("allElectron_wp70",          allElectrons.wp70,          "wp70[nallElectrons]/F");
  _outTree->Branch("allElectron_tIso",          allElectrons.tIso,          "tIso[nallElectrons]/F");
  _outTree->Branch("allElectron_eIso",          allElectrons.eIso,          "eIso[nallElectrons]/F");
  _outTree->Branch("allElectron_hIso",          allElectrons.hIso,          "hIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfChaIso",      allElectrons.pfChaIso,      "pfChaIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfChaPUIso",    allElectrons.pfChaPUIso,    "pfChaPUIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfPhoIso",      allElectrons.pfPhoIso,      "pfPhoIso[nallElectrons]/F");
  _outTree->Branch("allElectron_pfNeuIso",      allElectrons.pfNeuIso,      "pfNeuIso[nallElectrons]/F");

  // ------------------------------------------------------------------//

//===NOTE : for saving allJet outside of VH candidate (Duong 10-08-2015)=======
 
  _outTree->Branch("naJets",  	       	&naJets,  				"naJets/I");  
  _outTree->Branch("aJet_pt",				aJets.pt,					"pt[naJets]/F");
  _outTree->Branch("aJet_eta",			aJets.eta,				"eta[naJets]/F");
  _outTree->Branch("aJet_phi",			aJets.phi,				"phi[naJets]/F");
  _outTree->Branch("aJet_e",				aJets.e,					"e[naJets]/F");
  _outTree->Branch("aJet_pt_withJEC",	  aJets.pt_withJEC,			"pt_withJEC[naJets]/F");
  _outTree->Branch("aJet_eta_withJEC",  aJets.eta_withJEC,			"eta_withJEC[naJets]/F");
  _outTree->Branch("aJet_phi_withJEC",  aJets.phi_withJEC,			"phi_withJEC[naJets]/F");
  _outTree->Branch("aJet_e_withJEC",    aJets.e_withJEC,			"e_withJEC[naJets]/F");
  _outTree->Branch("aJet_pt_bestMC",	  aJets.pt_bestMC,			"pt_bestMC[naJets]/F");
  _outTree->Branch("aJet_eta_bestMC",	  aJets.eta_bestMC,			"eta_bestMC[naJets]/F");
  _outTree->Branch("aJet_phi_bestMC",	  aJets.phi_bestMC,			"phi_bestMC[naJets]/F");
  _outTree->Branch("aJet_e_bestMC",	  aJets.e_bestMC,			"e_bestMC[naJets]/F");
  _outTree->Branch("aJet_bestMCid",	  aJets.bestMCid,                       "bestMCid[naJets]/I");

  _outTree->Branch("aJet_jp",			aJets.jp, 			"jp[naJets]/F");
  _outTree->Branch("aJet_jpb",			aJets.jpb, 			"jpb[naJets]/F");
  _outTree->Branch("aJet_csv",			aJets.csv, 			"csv[naJets]/F");
  _outTree->Branch("aJet_csv_nominal",	        aJets.csv_nominal, 		"csv_nominal[naJets]/F");
  _outTree->Branch("aJet_csv_upBC",		aJets.csv_upBC, 		"csv_upBC[naJets]/F");
  _outTree->Branch("aJet_csv_downBC",		aJets.csv_downBC, 		"csv_downBC[naJets]/F");
  _outTree->Branch("aJet_csv_upL",		aJets.csv_upL, 			"csv_upL[naJets]/F");
  _outTree->Branch("aJet_csv_downL",		aJets.csv_downL, 		"csv_downL[naJets]/F");
  _outTree->Branch("aJet_csv_nominal4p",	aJets.csv_nominal4p, 		"csv_nominal4p[naJets]/F");
  _outTree->Branch("aJet_csv_upBC4p",		aJets.csv_upBC4p, 		"csv_upBC4p[naJets]/F");
  _outTree->Branch("aJet_csv_downBC4p",	        aJets.csv_downBC4p, 		"csv_downBC4p[naJets]/F");
  _outTree->Branch("aJet_csv_nominal1Bin",      aJets.csv_nominal1Bin, 	        "csv_nominal1Bin[naJets]/F");
  _outTree->Branch("aJet_csvivf",		aJets.csvivf, 			"csvivf[naJets]/F");
  _outTree->Branch("aJet_cmva",			aJets.cmva, 			"cmva[naJets]/F");
  _outTree->Branch("aJet_cosTheta",		aJets.cosTheta, 		"cosTheta[naJets]/F");
  _outTree->Branch("aJet_numTracksSV",	        aJets.numTracksSV, 		"numTracksSV[naJets]/I");
  _outTree->Branch("aJet_chf",			aJets.chf, 				"chf[naJets]/F");
  _outTree->Branch("aJet_nhf",			aJets.nhf, 				"nhf[naJets]/F");
  _outTree->Branch("aJet_cef",			aJets.cef, 				"cef[naJets]/F");
  _outTree->Branch("aJet_nef",			aJets.nef, 				"nef[naJets]/F");
  _outTree->Branch("aJet_nch",			aJets.nch, 				"nch[naJets]/F");
  _outTree->Branch("aJet_nconstituents",	aJets.nconstituents,		"nconstituents[naJets]");
  _outTree->Branch("aJet_flavour",		aJets.flavour, 			"flavour[naJets]/F");
  _outTree->Branch("aJet_isSemiLept",		aJets.isSemiLept, 		"isSemiLept[naJets]/I");
  _outTree->Branch("aJet_isSemiLeptMCtruth",	aJets.isSemiLeptMCtruth, 	"isSemiLeptMCtruth[naJets]/I");
  _outTree->Branch("aJet_SoftLeptpdgId",	aJets.SoftLeptpdgId,  	"SoftLeptpdgId[naJets]/I");
  _outTree->Branch("aJet_SoftLeptIdlooseMu", 	aJets.SoftLeptIdlooseMu,  	"SoftLeptIdlooseMu[naJets]/I");
  _outTree->Branch("aJet_SoftLeptId95", 	aJets.SoftLeptId95, 		"SoftLeptId95[naJets]/I");
  _outTree->Branch("aJet_SoftLeptPt", 	aJets.SoftLeptPt, 		"SoftLeptPt[naJets]/F");
  _outTree->Branch("aJet_SoftLeptdR", 	aJets.SoftLeptdR, 		"SoftLeptdR[naJets]/F");
  _outTree->Branch("aJet_SoftLeptptRel", 	aJets.SoftLeptptRel, 		"SoftLeptptRel[naJets]/F");
  _outTree->Branch("aJet_SoftLeptRelCombIso", aJets.SoftLeptRelCombIso,  	"SoftLeptRelCombIso[naJets]/F");
  _outTree->Branch("aJet_puJetIdL", 		aJets.puJetIdL,  			"puJetIdL[naJets]/F");
  _outTree->Branch("aJet_puJetIdM", 		aJets.puJetIdM,  			"puJetIdM[naJets]/F");
  _outTree->Branch("aJet_puJetIdT", 		aJets.puJetIdT,  			"puJetIdT[naJets]/F");
  _outTree->Branch("aJet_puJetIdMva", 	aJets.puJetIdMva,  		"puJetIdMva[naJets]/F");
  _outTree->Branch("aJet_charge",			aJets.charge, 			"charge[naJets]/F");
  _outTree->Branch("aJet_genPt",			aJets.genPt, 				"genPt[naJets]/F");
  _outTree->Branch("aJet_genEta",			aJets.genEta, 			"genEta[naJets]/F");
  _outTree->Branch("aJet_genPhi",			aJets.genPhi, 			"genPhi[naJets]/F");
  _outTree->Branch("aJet_JECUnc",			aJets.JECUnc, 			"JECUnc[naJets]/F");
  _outTree->Branch("aJet_vtxMass",		aJets.vtxMass, 			"vtxMass[naJets]/F");
  _outTree->Branch("aJet_vtx3dL",			aJets.vtx3dL, 			"vtx3dL[naJets]/F");
  _outTree->Branch("aJet_vtx3deL",		aJets.vtx3deL,	 		"vtx3deL[naJets]/F");
  _outTree->Branch("aJet_vtxProb",        aJets.vtxProb,            "vtxProb[naJets]/F");
  _outTree->Branch("aJet_ssvhe",          aJets.ssvhe,              "ssvhe[naJets]/F");
  _outTree->Branch("aJet_id",				aJets.id, 				"id[naJets]/b");
  _outTree->Branch("aJet_SF_CSVL",		aJets.SF_CSVL, 			"SF_CSVL[naJets]/F");
  _outTree->Branch("aJet_SF_CSVM",		aJets.SF_CSVM, 			"SF_CSVM[naJets]/F");
  _outTree->Branch("aJet_SF_CSVT",		aJets.SF_CSVT, 			"SF_CSVT[naJets]/F");
  _outTree->Branch("aJet_SF_CSVLerr",		aJets.SF_CSVLerr, 		"SF_CSVLerr[naJets]/F");
  _outTree->Branch("aJet_SF_CSVMerr",		aJets.SF_CSVMerr, 		"SF_CSVMerr[naJets]/F");
  _outTree->Branch("aJet_SF_CSVTerr",     aJets.SF_CSVTerr,         "SF_CSVTerr[naJets]/F");
  _outTree->Branch("aJet_tche",           aJets.tche,               "tche[naJets]/F");
  _outTree->Branch("aJet_tchp",           aJets.tchp,               "tchp[naJets]/F");
  _outTree->Branch("aJet_vtxPosition_x",  aJets.vtxPosition_x,      "vtxPosition_x[naJets]/F");
  _outTree->Branch("aJet_vtxPosition_y",  aJets.vtxPosition_y,      "vtxPosition_y[naJets]/F");
  _outTree->Branch("aJet_vtxPosition_z",  aJets.vtxPosition_z,      "vtxPosition_z[naJets]/F");
  _outTree->Branch("aJet_selectedTauDR",  aJets.selectedTauDR ,"selectedTauDR[naJets]/F");
 
//===================
//===NOTE : remove jetFat (Duong 10-11-2015)===
  _outTree->Branch("naJetsFat",  	&naJetsFat,		"naJetsFat/I"		);
/*
  _outTree->Branch("aJetFat_pt",	aJetsFat.pt,	"pt[naJetsFat]/F"	);
  _outTree->Branch("aJetFat_eta",	aJetsFat.eta,	"eta[naJetsFat]/F"	);
  _outTree->Branch("aJetFat_phi",	aJetsFat.phi,	"phi[naJetsFat]/F"	);
  _outTree->Branch("aJetFat_e",		aJetsFat.e,	"e[naJetsFat]/F"	);
  _outTree->Branch("aJetFat_csv",	aJetsFat.csv,	"csv[naJetsFat]/F"	);
*/
  _outTree->Branch("numJets",  		  &numJets,  		"numJets/I"       	);                
  _outTree->Branch("numBJets",		  &numBJets,  		"numBJets/I"       	);                
  _outTree->Branch("deltaPullAngle",  &deltaPullAngle,  "deltaPullAngle/F"	);
  _outTree->Branch("deltaPullAngle2", &deltaPullAngle2, "deltaPullAngle2/F"	);
  _outTree->Branch("gendrcc", 		  &gendrcc,  		"gendrcc/F"			);
  _outTree->Branch("gendrbb", 		  &gendrbb,  		"gendrbb/F"			);
  _outTree->Branch("genZpt", 		  &genZpt,  		"genZpt/F"			);
  _outTree->Branch("genWpt", 		  &genWpt,  		"genWpt/F"			);
  _outTree->Branch("genHpt", 		  &genHpt,  		"genHpt/F"			);
  _outTree->Branch("weightSignalEWK", &weightSignalEWK, "weightSignalEWK/F" );
  _outTree->Branch("weightSignalQCD", &weightSignalQCD, "weightSignalQCD/F" );

  _outTree->Branch("weightTrig"        , &weightTrig          ,  "weightTrig/F");
  _outTree->Branch("weightTrigMay"        , &weightTrigMay          ,  "weightTrigMay/F");
  _outTree->Branch("weightTrigV4"        , &weightTrigV4          ,  "weightTrigV4/F");
  _outTree->Branch("weightTrigMET"        , &weightTrigMET          ,  "weightTrigMET/F");
  _outTree->Branch("weightTrigOrMu30"  , &weightTrigOrMu30    , "weightTrigOrMu30/F");
  _outTree->Branch("weightEleRecoAndId"        , &weightEleRecoAndId     ,  "weightEleRecoAndId/F");
  _outTree->Branch("weightEleTrigJetMETPart"        , &weightEleTrigJetMETPart          ,  "weightEleTrigJetMETPart/F");
  _outTree->Branch("weightEleTrigElePart"        , &weightEleTrigElePart          ,  "weightEleTrigElePart/F");
  _outTree->Branch("weightEleTrigEleAugPart"        , &weightEleTrigEleAugPart          ,  "weightEleTrigEleAugPart/F");

  _outTree->Branch("weightTrigMET80"        , &weightTrigMET80          , "weightTrigMET80/F");
  _outTree->Branch("weightTrigMET100"        , &weightTrigMET100          , "weightTrigMET100/F");
  _outTree->Branch("weightTrig2CJet20"        , &weightTrig2CJet20          , "weightTrig2CJet20/F");
  _outTree->Branch("weightTrigMET150"        , &weightTrigMET150          , "weightTrigMET150/F");
  _outTree->Branch("weightTrigMET802CJet"        , &weightTrigMET802CJet          , "weightTrigMET802CJet/F");
  _outTree->Branch("weightTrigMET1002CJet"        , &weightTrigMET1002CJet          , "weightTrigMET1002CJet/F");
  _outTree->Branch("weightTrigMETLP"        , &weightTrigMETLP          , "weightTrigMETLP/F");

  _outTree->Branch("weightTrig2012A", &weightTrig2012A,"weightTrig2012A/F");
  _outTree->Branch("weightTrig2012ADiMuon", &weightTrig2012ADiMuon,"weightTrig2012ADiMuon/F");
  _outTree->Branch("weightTrig2012ADiEle", &weightTrig2012ADiEle,"weightTrig2012ADiEle/F");
  _outTree->Branch("weightTrig2012ASingleMuon", &weightTrig2012ASingleMuon,"weightTrig2012ASingleMuon/F");
  _outTree->Branch("weightTrig2012ASingleEle", &weightTrig2012ASingleEle,"weightTrig2012ASingleEle/F");
  _outTree->Branch("weightTrig2012AMuonPlusWCandPt", &weightTrig2012AMuonPlusWCandPt,"weightTrig2012AMuonPlusWCandPt/F");

  _outTree->Branch("weightTrig2012", &weightTrig2012,"weightTrig2012/F");
  _outTree->Branch("weightTrig2012DiMuon", &weightTrig2012DiMuon,"weightTrig2012DiMuon/F");
  _outTree->Branch("weightTrig2012DiEle", &weightTrig2012DiEle,"weightTrig2012DiEle/F");
  _outTree->Branch("weightTrig2012SingleMuon", &weightTrig2012SingleMuon,"weightTrig2012SingleMuon/F");
  _outTree->Branch("weightTrig2012SingleEle", &weightTrig2012SingleEle,"weightTrig2012SingleEle/F");
  _outTree->Branch("weightTrig2012MuonPlusWCandPt", &weightTrig2012MuonPlusWCandPt,"weightTrig2012MuonPlusWCandPt/F");

  _outTree->Branch("weightTrig2012AB", &weightTrig2012AB,"weightTrig2012AB/F");
  _outTree->Branch("weightTrig2012ABDiMuon", &weightTrig2012ABDiMuon,"weightTrig2012ABDiMuon/F");
  _outTree->Branch("weightTrig2012ABDiEle", &weightTrig2012ABDiEle,"weightTrig2012ABDiEle/F");
  _outTree->Branch("weightTrig2012ABSingleMuon", &weightTrig2012ABSingleMuon,"weightTrig2012ABSingleMuon/F");
  _outTree->Branch("weightTrig2012ABSingleEle", &weightTrig2012ABSingleEle,"weightTrig2012ABSingleEle/F");
  _outTree->Branch("weightTrig2012ABMuonPlusWCandPt", &weightTrig2012ABMuonPlusWCandPt,"weightTrig2012ABMuonPlusWCandPt/F");
  
  _outTree->Branch("weightTrig2012DiJet30MHT80", &weightTrig2012DiJet30MHT80,"weightTrig2012DiJet30MHT80/F");
  _outTree->Branch("weightTrig2012PFMET150", &weightTrig2012PFMET150,"weightTrig2012PFMET150/F");
  _outTree->Branch("weightTrig2012SumpT100MET100", &weightTrig2012SumpT100MET100,"weightTrig2012SumpT100MET100/F");
  _outTree->Branch("weightTrig2012APFMET150orDiJetMET", &weightTrig2012APFMET150orDiJetMET,"weightTrig2012APFMET150orDiJetMET/F");
  _outTree->Branch("weightTrig2012BPFMET150orDiJetMET", &weightTrig2012BPFMET150orDiJetMET,"weightTrig2012BPFMET150orDiJetMET/F");
  _outTree->Branch("weightTrig2012CPFMET150orDiJetMET", &weightTrig2012CPFMET150orDiJetMET,"weightTrig2012CPFMET150orDiJetMET/F");

  _outTree->Branch("deltaPullAngleAK7", &deltaPullAngleAK7  ,  "deltaPullAngleAK7/F");
  _outTree->Branch("deltaPullAngle2AK7", &deltaPullAngle2AK7  ,  "deltaPullAngle2AK7/F");
  _outTree->Branch("PU0",       	&PU0,  			"PU0/F"				);
  _outTree->Branch("PUm1",      	&PUm1,  		"PUm1/F"			);
  _outTree->Branch("PUp1",      	&PUp1,  		"PUp1/F"			);
  _outTree->Branch("PUweight",      &PUweight,  	"PUweight/F"		);
  _outTree->Branch("PUweightP",     &PUweightP,  	"PUweightP/F"	  	);
  _outTree->Branch("PUweightM",     &PUweightM,  	"PUweightM/F" 	  	);
  _outTree->Branch("PUweightAB",    &PUweightAB,  	"PUweightAB/F"	  	);
  _outTree->Branch("PUweight2011B",	&PUweight2011B, "PUweight2011B/F" 	);
  _outTree->Branch("PUweight1DObs", &PUweight1DObs, "PUweight1DObs/F" 	);
  _outTree->Branch("eventFlav",     &eventFlav,     "eventFlav/I"	  	);
 
  _outTree->Branch("Vtype",		    &Vtype, 		"Vtype/I" 		);                
  _outTree->Branch("VtypeWithTau",  &VtypeWithTau, 	"VtypeWithTau/I");                
  _outTree->Branch("HVdPhi",  		&HVdPhi,    	"HVdPhi/F" 		);                
  _outTree->Branch("HVMass",  		&HVMass,   		"HVMass/F" 		);                
  _outTree->Branch("HMETdPhi",  	&HMETdPhi,   	"HMETdPhi/F" 	);                
  _outTree->Branch("VMt",  			&VMt,   		"VMt/F"    		);             	

  _outTree->Branch("nvlep",    &nvlep,    "nvlep/I"    );

//===NOTE : don't save vlep (Duong 10-09-2015)
/*
  _outTree->Branch("vLepton_mass",vLeptons.mass ,"mass[nvlep]/F");
  _outTree->Branch("vLepton_pt",vLeptons.pt ,"pt[nvlep]/F");
  _outTree->Branch("vLepton_eta",vLeptons.eta ,"eta[nvlep]");
  _outTree->Branch("vLepton_phi",vLeptons.phi ,"phi[nvlep]/F");
  _outTree->Branch("vLepton_aodCombRelIso",vLeptons.aodCombRelIso ,"aodCombRelIso[nvlep]/F");
  _outTree->Branch("vLepton_pfCombRelIso",vLeptons.pfCombRelIso ,"pfCombRelIso[nvlep]/F");
  _outTree->Branch("vLepton_photonIso",vLeptons.photonIso ,"photonIso[nvlep]/F");
  _outTree->Branch("vLepton_neutralHadIso",vLeptons.neutralHadIso ,"neutralHadIso[nvlep]/F");
  _outTree->Branch("vLepton_chargedHadIso",vLeptons.chargedHadIso ,"chargedHadIso[nvlep]/F");
  _outTree->Branch("vLepton_chargedPUIso",vLeptons.chargedPUIso ,"chargedPUIso[nvlep]/F");
  _outTree->Branch("vLepton_particleIso",vLeptons.particleIso ,"particleIso[nvlep]/F");
  _outTree->Branch("vLepton_dxy",vLeptons.dxy ,"dxy[nvlep]/F");
  _outTree->Branch("vLepton_dz",vLeptons.dz ,"dz[nvlep]/F");
  _outTree->Branch("vLepton_type",vLeptons.type ,"type[nvlep]/I");
  _outTree->Branch("vLepton_id80",vLeptons.id80 ,"id80[nvlep]/F");
  _outTree->Branch("vLepton_id95",vLeptons.id95 ,"id95[nvlep]/F");
  _outTree->Branch("vLepton_vbtf",vLeptons.vbtf ,"vbtf[nvlep]/F");
  _outTree->Branch("vLepton_id80NoIso",vLeptons.id80NoIso ,"id80NoIso[nvlep]/F");
  _outTree->Branch("vLepton_genPt",vLeptons.genPt ,"genPt[nvlep]/F");
  _outTree->Branch("vLepton_genEta",vLeptons.genEta ,"genEta[nvlep]");
  _outTree->Branch("vLepton_genPhi",vLeptons.genPhi ,"genPhi[nvlep]/F");
  _outTree->Branch("vLepton_charge",vLeptons.charge ,"charge[nvlep]/F");
  _outTree->Branch("vLepton_pfCorrIso",vLeptons.pfCorrIso,"pfCorrIso[nvlep]/F");
  _outTree->Branch("vLepton_pfCorrIsoHCP",vLeptons.pfCorrIsoHCP,"pfCorrIsoHCP[nvlep]/F");
  _outTree->Branch("vLepton_id2012tight",vLeptons.id2012tight,"id2012tight[nvlep]/F");
  _outTree->Branch("vLepton_idMVAnotrig",vLeptons.idMVAnotrig,"idMVAnotrig[nvlep]/F");
  _outTree->Branch("vLepton_idMVAtrig",vLeptons.idMVAtrig,"idMVAtrig[nvlep]/F");
  _outTree->Branch("vLepton_idMVApresel",vLeptons.idMVApresel,"idMVApresel[nvlep]/F");
  _outTree->Branch("vLepton_innerHits",vLeptons.innerHits,"innerHits[nvlep]/F");
  _outTree->Branch("vLepton_photonIsoDoubleCount",vLeptons.photonIsoDoubleCount,"photonIsoDoubleCount[nvlep]/F");
  _outTree->Branch("vLepton_wpHWW",vLeptons.wpHWW,"wpHWW[nvlep]/F");
  _outTree->Branch("vLepton_wp95",vLeptons.wp95,"wp95[nvlep]/F");
  _outTree->Branch("vLepton_wp90",vLeptons.wp90,"wp90[nvlep]/F");
  _outTree->Branch("vLepton_wp85",vLeptons.wp85,"wp85[nvlep]/F");
  _outTree->Branch("vLepton_wp80",vLeptons.wp80,"wp80[nvlep]/F");
  _outTree->Branch("vLepton_wp70",vLeptons.wp70,"wp70[nvlep]/F");
*/

  // Adding variables for tau leptons
  _outTree->Branch("nvlepTau", &nvlepTau, "nvlepTau/I" );

//==NOTE : don't save vLepton== Duong (10-09-2015)
/*  
  _outTree->Branch("vLeptonTaus_mass",vLeptonsTaus.mass ,"mass[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_pt",vLeptonsTaus.pt ,"pt[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_eta",vLeptonsTaus.eta ,"eta[nvlepTau]");
  _outTree->Branch("vLeptonTaus_phi",vLeptonsTaus.phi ,"phi[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_pfCombRelIso",vLeptonsTaus.pfCombRelIso ,"pfCombRelIso[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_genPt",vLeptonsTaus.genPt ,"genPt[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_genEta",vLeptonsTaus.genEta ,"genEta[nvlepTau]");
  _outTree->Branch("vLeptonTaus_genPhi",vLeptonsTaus.genPhi ,"genPhi[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_charge",vLeptonsTaus.charge ,"charge[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_decayModeFinding",vLeptonsTaus.decayModeFinding,"decayModeFinding[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byLooseCombinedIsolationDeltaBetaCorr",vLeptonsTaus.byLooseCombinedIsolationDeltaBetaCorr,"byLooseCombinedIsolationDeltaBetaCorr[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstMuonTight",vLeptonsTaus.againstMuonTight,"againstMuonTight[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronLoose",vLeptonsTaus.againstElectronLoose,"againstElectronLoose[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronMedium",vLeptonsTaus.againstElectronMedium,"againstElectronMedium[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronMVA",vLeptonsTaus.againstElectronMVA,"againstElectronMVA[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_NsignalPFChargedHadrCands",vLeptonsTaus.NsignalPFChargedHadrCands,"NsignalPFChargedHadrCands[nvlepTau]/I");
  _outTree->Branch("vLeptonTaus_NsignalPFGammaCands",vLeptonsTaus.NsignalPFGammaCands,"NsignalPFGammaCands[nvlepTau]/I");
  _outTree->Branch("vLeptonTaus_leadPFChargedHadrCandPt",vLeptonsTaus.leadPFChargedHadrCandPt,"leadPFChargedHadrCandPt[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byLooseIsolation",vLeptonsTaus.byLooseIsolation,"byLooseIsolation[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byMediumIsolation",vLeptonsTaus.byMediumIsolation,"byMediumIsolation[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byTightIsolation",vLeptonsTaus.byTightIsolation,"byTightIsolation[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byLooseCombinedIsolationDeltaBetaCorr3Hits",vLeptonsTaus.byLooseCombinedIsolationDeltaBetaCorr3Hits,"byLooseCombinedIsolationDeltaBetaCorr3Hits[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byMediumCombinedIsolationDeltaBetaCorr3Hits",vLeptonsTaus.byMediumCombinedIsolationDeltaBetaCorr3Hits,"byMediumCombinedIsolationDeltaBetaCorr3Hits[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byTightCombinedIsolationDeltaBetaCorr3Hits",vLeptonsTaus.byTightCombinedIsolationDeltaBetaCorr3Hits,"byTightCombinedIsolationDeltaBetaCorr3Hits[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronMVA3raw",vLeptonsTaus.againstElectronMVA3raw,"againstElectronMVA3raw[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronMVA3category",vLeptonsTaus.againstElectronMVA3category,"againstElectronMVA3category[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronLooseMVA3",vLeptonsTaus.againstElectronLooseMVA3,"againstElectronLooseMVA3[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronMediumMVA3",vLeptonsTaus.againstElectronMediumMVA3,"againstElectronMediumMVA3[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronTightMVA3",vLeptonsTaus.againstElectronTightMVA3,"againstElectronTightMVA3[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronVTightMVA3",vLeptonsTaus.againstElectronVTightMVA3,"againstElectronVTightMVA3[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstElectronDeadECAL",vLeptonsTaus.againstElectronDeadECAL,"againstElectronDeadECAL[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byLooseIsolationMVA",vLeptonsTaus.byLooseIsolationMVA,"byLooseIsolationMVA[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byMediumIsolationMVA",vLeptonsTaus.byMediumIsolationMVA,"byMediumIsolationMVA[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byTightIsolationMVA",vLeptonsTaus.byTightIsolationMVA,"byTightIsolationMVA[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byLooseIsolationMVA2",vLeptonsTaus.byLooseIsolationMVA2,"byLooseIsolationMVA2[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byMediumIsolationMVA2",vLeptonsTaus.byMediumIsolationMVA2,"byMediumIsolationMVA2[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_byTightIsolationMVA2",vLeptonsTaus.byTightIsolationMVA2,"byTightIsolationMVA2[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstMuonLoose2",vLeptonsTaus.againstMuonLoose2,"againstMuonLoose2[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstMuonMedium2",vLeptonsTaus.againstMuonMedium2,"againstMuonMedium2[nvlepTau]/F");
  _outTree->Branch("vLeptonTaus_againstMuonTight2",vLeptonsTaus.againstMuonTight2,"againstMuonTight2[nvlepTau]/F");
*/
  _outTree->Branch("tauPlusMode"          ,  &tauPlusMode    ,   "tauPlusMode/I");
  _outTree->Branch("tauMinusMode"         ,  &tauMinusMode   ,   "tauMinusMode/I");

  

  _outTree->Branch("nalep",    &nalep,    "nalep/I"    );
  _outTree->Branch("aLepton_mass",aLeptons.mass ,"mass[nalep]/F");
  _outTree->Branch("aLepton_pt",aLeptons.pt ,"pt[nalep]/F");
  _outTree->Branch("aLepton_eta",aLeptons.eta ,"eta[nalep]");
  _outTree->Branch("aLepton_phi",aLeptons.phi ,"phi[nalep]/F");
  _outTree->Branch("aLepton_aodCombRelIso",aLeptons.aodCombRelIso ,"aodCombRelIso[nalep]/F");
  _outTree->Branch("aLepton_pfCombRelIso",aLeptons.pfCombRelIso ,"pfCombRelIso[nalep]/F");
  _outTree->Branch("aLepton_photonIso",aLeptons.photonIso ,"photonIso[nalep]/F");
  _outTree->Branch("aLepton_neutralHadIso",aLeptons.neutralHadIso ,"neutralHadIso[nalep]/F");
  _outTree->Branch("aLepton_chargedHadIso",aLeptons.chargedHadIso ,"chargedHadIso[nalep]/F");
  _outTree->Branch("aLepton_chargedPUIso",aLeptons.chargedPUIso ,"chargedPUIso[nalep]/F");
  _outTree->Branch("aLepton_particleIso",aLeptons.particleIso ,"particleIso[nalep]/F");
  _outTree->Branch("aLepton_dxy",aLeptons.dxy ,"dxy[nalep]/F");
  _outTree->Branch("aLepton_dz",aLeptons.dz ,"dz[nalep]/F");
  _outTree->Branch("aLepton_type",aLeptons.type ,"type[nalep]/I");
  _outTree->Branch("aLepton_id80",aLeptons.id80 ,"id80[nalep]/F");
  _outTree->Branch("aLepton_id95",aLeptons.id95 ,"id95[nalep]/F");
  _outTree->Branch("aLepton_vbtf",aLeptons.vbtf ,"vbtf[nalep]/F");
  _outTree->Branch("aLepton_id80NoIso",aLeptons.id80NoIso ,"id80NoIso[nalep]/F");
  _outTree->Branch("aLepton_genPt",aLeptons.genPt ,"genPt[nalep]/F");
  _outTree->Branch("aLepton_genEta",aLeptons.genEta ,"genEta[nalep]");
  _outTree->Branch("aLepton_genPhi",aLeptons.genPhi ,"genPhi[nalep]/F");
  _outTree->Branch("aLepton_charge",aLeptons.charge ,"charge[nalep]/F");
  _outTree->Branch("aLepton_pfCorrIso",aLeptons.pfCorrIso,"pfCorrIso[nalep]/F");
  _outTree->Branch("aLepton_pfCorrIsoHCP",aLeptons.pfCorrIsoHCP,"pfCorrIsoHCP[nalep]/F");
  _outTree->Branch("aLepton_id2012tight",aLeptons.id2012tight,"id2012tight[nalep]/F");
  _outTree->Branch("aLepton_idMVAnotrig",aLeptons.idMVAnotrig,"idMVAnotrig[nalep]/F");
  _outTree->Branch("aLepton_idMVAtrig",aLeptons.idMVAtrig,"idMVAtrig[nalep]/F");
  _outTree->Branch("aLepton_idMVApresel",aLeptons.idMVApresel,"idMVApresel[nalep]/F");
  _outTree->Branch("aLepton_innerHits",aLeptons.innerHits,"innerHits[nalep]/F");
  _outTree->Branch("aLepton_photonIsoDoubleCount",aLeptons.photonIsoDoubleCount,"photonIsoDoubleCount[nalep]/F");
  _outTree->Branch("aLepton_wpHWW",aLeptons.wpHWW,"wpHWW[nalep]/F");
  _outTree->Branch("aLepton_wp95",aLeptons.wp95,"wp95[nalep]/F");
  _outTree->Branch("aLepton_wp90",aLeptons.wp90,"wp90[nalep]/F");
  _outTree->Branch("aLepton_wp85",aLeptons.wp85,"wp85[nalep]/F");
  _outTree->Branch("aLepton_wp80",aLeptons.wp80,"wp80[nalep]/F");
  _outTree->Branch("aLepton_wp70",aLeptons.wp70,"wp70[nalep]/F");

  _outTree->Branch("top"		,  &top	         ,   "mass/F:pt/F:wMass/F");
  _outTree->Branch("WplusMode"		,  &WplusMode	 ,   "WplusMode/I");
  _outTree->Branch("WminusMode"		,  &WminusMode	 ,   "WminusMode/I");

  //IVF
  _outTree->Branch("nSvs",&nSvs ,"nSvs/I");
  _outTree->Branch("Sv_massBCand", &IVF.massBcand,"massBcand[nSvs]/F");
  _outTree->Branch("Sv_massSv", &IVF.massSv,"massSv[nSvs]/F");
  _outTree->Branch("Sv_pt", &IVF.pt,"pt[nSvs]/F");
  _outTree->Branch("Sv_eta", &IVF.eta,"eta[nSvs]/F");
  _outTree->Branch("Sv_phi", &IVF.phi,"phi[nSvs]/F");
  _outTree->Branch("Sv_dist3D", &IVF.dist3D,"dist3D[nSvs]/F");
  _outTree->Branch("Sv_dist2D", &IVF.dist2D,"dist2D[nSvs]/F");
  _outTree->Branch("Sv_distSim2D", &IVF.distSig2D,"distSig2D[nSvs]/F");
  _outTree->Branch("Sv_distSig3D", &IVF.distSig3D,"distSig3D[nSvs]/F");
  _outTree->Branch("Sv_dist3D_norm", &IVF.dist3D_norm,"dist3D_norm[nSvs]/F");
  //IVF higgs candidate
  _outTree->Branch("SVH"          ,  &SVH               ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");

  //   //SimBHadron
  _outTree->Branch("nSimBs",&nSimBs ,"nSimBs/I");
  _outTree->Branch("SimBs_mass", &SimBs.mass,"mass[nSimBs]/F");
  _outTree->Branch("SimBs_pt", &SimBs.pt,"pt[nSimBs]/F");
  _outTree->Branch("SimBs_eta", &SimBs.eta,"eta[nSimBs]/F");
  _outTree->Branch("SimBs_phi", &SimBs.phi,"phi[nSimBs]/F");
  _outTree->Branch("SimBs_vtx_x", &SimBs.vtx_x,"vtx_x[nSimBs]/F");
  _outTree->Branch("SimBs_vtx_y", &SimBs.vtx_y,"vtx_y[nSimBs]/F");
  _outTree->Branch("SimBs_vtx_z", &SimBs.vtx_z,"vtx_z[nSimBs]/F");
  _outTree->Branch("SimBs_pdgId", &SimBs.pdgId,"pdgId[nSimBs]/F");
  _outTree->Branch("SimBs_status", &SimBs.status,"status[nSimBs]/F");
  //SimBHadron Higgs Candidate
  _outTree->Branch("SimBsH"          ,  &SimBsH               ,  "mass/F:pt/F:eta:phi/F:dR/F:dPhi/F:dEta/F");

  _outTree->Branch("rho"		,  &rho	         ,   "rho/F");
  _outTree->Branch("rho25"		,  &rho25	         ,   "rho25/F");
  _outTree->Branch("rhoN"		,  &rhoN	         ,   "rhoN/F");
  _outTree->Branch("nPVs"		,  &nPVs	         ,   "nPVs/I");
  _outTree->Branch("nPV1s"		,  &nPV1s	         ,   "nPV1s/I");
  _outTree->Branch("nTruePUs"		,  &nTruePUs	         ,   "nTruePUs/I");
  _outTree->Branch("METnoPU"		,  &METnoPU	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METnoPUCh"		,  &METnoPUCh	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MET"        ,  &MET          ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("pfMET"        ,  &pfMET          ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METtype1corr"		,  &METtype1corr	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METtype1p2corr"		,  &METtype1p2corr	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METnoPUtype1corr"		,  &METnoPUtype1corr	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METnoPUtype1p2corr"		,  &METnoPUtype1p2corr	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("METtype1diff"		,  &METtype1diff	         ,   "et/F:sumet:sig/F:phi/F");
  
//==NOTE : MET outside of VH candidate==
  _outTree->Branch("aMET"        ,  &aMET          ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("apfMET"        ,  &apfMET          ,   "et/F:sumet:sig/F:phi/F");

  _outTree->Branch("metUnc_et",&metUnc.et ,"et[24]/F");
  _outTree->Branch("metUnc_phi",&metUnc.phi ,"phi[24]/F");
  _outTree->Branch("metUnc_sumet",&metUnc.sumet ,"sumet[24]/F");

  _outTree->Branch("fakeMET"		,  &fakeMET	         ,   "et/F:sumet:sig/F:phi/F");
  _outTree->Branch("MHT"		,  &MHT	         ,   "mht/F:ht:sig/F:phi/F");
  _outTree->Branch("minDeltaPhijetMET"		,  &minDeltaPhijetMET	         ,   "minDeltaPhijetMET/F");
  _outTree->Branch("jetPt_minDeltaPhijetMET"		,  &jetPt_minDeltaPhijetMET	         ,   "jetPt_minDeltaPhijetMET/F");

  std::stringstream s;
  s << "triggerFlags[" << triggers.size() << "]/b";
  _outTree->Branch("triggerFlags", triggerFlags, s.str().c_str()); 
 
  _outTree->Branch("EVENT"		,  &EVENT	         ,   "run/I:lumi/I:event/I:json/I");
  _outTree->Branch("hbhe"		,  &hbhe	         ,   "hbhe/b");
  _outTree->Branch("totalKinematics"		,  &totalKinematics	         ,   "totalKinematics/b");
  _outTree->Branch("ecalFlag"		,  &ecalFlag	         ,   "ecalFlag/b");
  _outTree->Branch("cschaloFlag"		,  &cschaloFlag	         ,   "cschaloFlag/b");
  _outTree->Branch("hcallaserFlag"		,  &hcallaserFlag	         ,   "hcallaserFlag/b");
  _outTree->Branch("trackingfailureFlag"		,  &trackingfailureFlag	         ,   "trackingfailureFlag/b");
  _outTree->Branch("eebadscFlag"		,  &eebadscFlag	         ,   "eebadscFlag/b");
  _outTree->Branch("btag1TSF"		,  &btag1TSF	         ,   "btag1TSF/F");
  _outTree->Branch("btag2TSF"		,  &btag2TSF	         ,   "btag2TSF/F");
  _outTree->Branch("btag1T2CSF"	,  &btag1T2CSF	         ,   "btag1T2CSF/F");
  _outTree->Branch("btag2CSF"	,  &btag2CSF	         ,   "btag2CSF/F");
  _outTree->Branch("btagA0CSF"	,  &btagA0CSF	         ,   "btagA0CSF/F");
  _outTree->Branch("btagA0TSF"	,  &btagA0TSF	         ,   "btagA0TSF/F");
  _outTree->Branch("btag1TA1C"	,  &btag1TA1C	         ,   "btag1TA1C/F");

//==NOTE : some gen info stuffs (Duong 10-21-2015)==
  _outTree->Branch("zdecayMode",    			&zdecayMode,				"zdecayMode/I"    );
  _outTree->Branch("nGenLep",    			&nGenLep,				"nGenLep/I"    );
  _outTree->Branch("nGenPho",    			&nGenPho,				"nGenPho/I"    );
  _outTree->Branch("nGenSta3obj",    			&nGenSta3obj,				"nGenSta3obj/I"    );
  _outTree->Branch("nGenAk5Jet",    			&nGenAk5Jet,				"nGenAk5Jet/I"    );
  _outTree->Branch("nGenPatPFjet",    			&nGenPatPFjet,				"nGenPatPFjet/I"    );
//==TODO: add mass or energy so that we can build the 4 vector==
  _outTree->Branch("genLep_pt", genLeps.pt, "genLep_pt[nGenLep]/F") ;  
  _outTree->Branch("genLep_eta", genLeps.eta, "genLep_eta[nGenLep]/F") ;  
  _outTree->Branch("genLep_phi", genLeps.phi, "genLep_phi[nGenLep]/F") ;  
  _outTree->Branch("genLep_charge", genLeps.charge, "genLep_charge[nGenLep]/I") ;  
  _outTree->Branch("genLep_flavor", genLeps.flavor, "genLep_flavor[nGenLep]/I") ;  
  _outTree->Branch("genLep_status", genLeps.status, "genLep_status[nGenLep]/I") ;  
  _outTree->Branch("genPho_pt", genPhos.pt, "genPho_pt[nGenPho]/F") ;  
  _outTree->Branch("genPho_eta", genPhos.eta, "genPho_eta[nGenPho]/F") ;  
  _outTree->Branch("genPho_phi", genPhos.phi, "genPho_phi[nGenPho]/F") ;  
  _outTree->Branch("genPho_charge", genPhos.charge, "genPho_charge[nGenPho]/I") ;  
  _outTree->Branch("genPho_flavor", genPhos.flavor, "genPho_flavor[nGenPho]/I") ;  
  _outTree->Branch("genPho_status", genPhos.status, "genPho_status[nGenPho]/I") ;  
  _outTree->Branch("genSta3obj_pt", genSta3objs.pt, "genSta3obj_pt[nGenSta3obj]/F") ;  
  _outTree->Branch("genSta3obj_eta", genSta3objs.eta, "genSta3obj_eta[nGenSta3obj]/F") ;  
  _outTree->Branch("genSta3obj_phi", genSta3objs.phi, "genSta3obj_phi[nGenSta3obj]/F") ;  
  _outTree->Branch("genSta3obj_charge", genSta3objs.charge, "genSta3obj_charge[nGenSta3obj]/I") ;  
  _outTree->Branch("genSta3obj_flavor", genSta3objs.flavor, "genSta3obj_flavor[nGenSta3obj]/I") ;  
  _outTree->Branch("genSta3obj_status", genSta3objs.status, "genSta3obj_status[nGenSta3obj]/I") ;  
  _outTree->Branch("genAk5Jet_pt", genAk5Jets.pt, "genAk5Jet_pt[nGenAk5Jet]/F") ;  
  _outTree->Branch("genAk5Jet_eta", genAk5Jets.eta, "genAk5Jet_eta[nGenAk5Jet]/F") ;  
  _outTree->Branch("genAk5Jet_phi", genAk5Jets.phi, "genAk5Jet_phi[nGenAk5Jet]/F") ;  
  _outTree->Branch("genAk5Jet_charge", genAk5Jets.charge, "genAk5Jet_charge[nGenAk5Jet]/I") ;  
  _outTree->Branch("genAk5Jet_flavor", genAk5Jets.flavor, "genAk5Jet_flavor[nGenAk5Jet]/I") ;  
  _outTree->Branch("genAk5Jet_status", genAk5Jets.status, "genAk5Jet_status[nGenAk5Jet]/I") ;  
  _outTree->Branch("genPatPFjet_pt", genPatPFjets.pt, "genPatPFjet_pt[nGenPatPFjet]/F") ;  
  _outTree->Branch("genPatPFjet_eta", genPatPFjets.eta, "genPatPFjet_eta[nGenPatPFjet]/F") ;  
  _outTree->Branch("genPatPFjet_phi", genPatPFjets.phi, "genPatPFjet_phi[nGenPatPFjet]/F") ;  
  _outTree->Branch("genPatPFjet_charge", genPatPFjets.charge, "genPatPFjet_charge[nGenPatPFjet]/I") ;  
  _outTree->Branch("genPatPFjet_flavor", genPatPFjets.flavor, "genPatPFjet_flavor[nGenPatPFjet]/I") ;  
  _outTree->Branch("genPatPFjet_status", genPatPFjets.status, "genPatPFjet_status[nGenPatPFjet]/I") ; 

  int ievt=0;  
  int totalcount=0;


  //
  //  LOOP OVER ALL FILES, EVENTS
  //
  
  //  TFile* inFile = new TFile(inputFile.c_str(), "read");
  for(unsigned int iFile=0; iFile<inputFiles_.size(); ++iFile)
  {
	// Open input file
    std::cout << iFile << std::endl;
    TFile* inFile = TFile::Open(inputFiles_[iFile].c_str());
    if(inFile==0) {std::cout << "FAILED " << inputFiles_[iFile] << std::endl; continue; }
    else std::cout << "File " << iFile << " opened successfully." << std::endl;

    // loop the events      
    fwlite::Event ev(inFile);
    fwlite::Handle< VHbbEventAuxInfo > vhbbAuxHandle; 
    for(ev.toBegin(); !ev.atEnd() ; ++ev, ++ievt)
    {
      
      hEvent->Fill(0.5) ;
      
      if (ievt < skipEvents_ && skipEvents_ > 0) continue;
      if (maxEvents_ >= 0)
        if (ievt > maxEvents_ + skipEvents_) break;
	  // Get HbbAnalyzer Object from Step1 PATTuble
      const char * lab = "HbbAnalyzerNew";
      vhbbAuxHandle.getByLabel(ev,lab,0,0);
      const VHbbEventAuxInfo & aux = *vhbbAuxHandle.product();

	  // Extract Event Info (EVENT = EventInfo object)
      EVENT.run = ev.id().run();
      EVENT.lumi = ev.id().luminosityBlock();
      EVENT.event = ev.id().event();
      EVENT.json = jsonContainsEvent (jsonVector, ev);
      isBadHcalEvent = std::binary_search (badEvents.begin(), badEvents.end(), RLE( EVENT.run,EVENT.lumi,EVENT.event),evcomp);
      weightMCProd = aux.weightMCProd;
      weightSignalQCD=1.;
      weightSignalEWK=1.;

      if(EVENT.run < runMin_ && runMin_ > 0) continue;
      if(EVENT.run > runMax_ && runMax_ > 0) continue;

      hEvent->Fill(1.5) ; //after skip
      
      count->Fill(1.);



//===NOTE : reset all objects which is going to be filled into tree (Duong 10-16-2015)===
      
      naJets=0 ;
      nhJets=0 ;
      nfathFilterJets=0 ;
      naJetsFat=0 ;
      nallJets = 0 ;
      nallMuons = 0 ;
      nallElectrons = 0 ;
      MET.reset() ;
      pfMET.reset() ;
      fakeMET.reset() ;
      METnoPU.reset() ;
      METnoPUCh.reset() ;
      METtype1corr.reset() ;
      METtype1p2corr.reset() ;
      METnoPUtype1corr.reset() ;
      METnoPUtype1p2corr.reset() ;
      METtype1diff.reset() ;

      nvlep=0 ;
      nalep=0 ;
      nvlepTau=0 ;

      //Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      //event.getByLabel(edm::InputTag("addPileupInfo"), PupInfo);
      //std::vector<PileupSummaryInfo>::const_iterator PVI;
      //float Tnpv = -1;
      //for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI)
	  //{
      //  int BX = PVI->getBunchCrossing();
      //  if(BX == 0)
	  //  { 
      //    Tnpv = PVI->getTrueNumInteractions();
      //    continue;
      //  }
      //}
      //double MyWeight = LumiWeights_.weight( Tnpv );
      //  double MyWeight = LumiWeights_.weight( (*iEventB) );
      nPV1s=aux.pvInfo.nVertices ;
      nTruePUs=aux.puInfo.truePU ;
      PUweight=1.;
      PUweightP=1.;
      PUweightM=1.;
      PUweightAB=1.;
      PUweight2011B=1.;
      PUweight1DObs=1.;
      weightSignalQCD=1;
      if(isMC_)
      {
		// PU weights // Run2011A
		std::map<int, unsigned int>::const_iterator puit = aux.puInfo.pus.find(0);
		int npu =puit->second;
		PUweight =   lumiWeights.weight( (int) aux.puInfo.truePU ); //use new method with "true PU"        
		PUweightP =  lumiWeightsPl.weight( (int) aux.puInfo.truePU ); //use new method with "true PU"        
		PUweightM =  lumiWeightsMi.weight( (int) aux.puInfo.truePU ); //use new method with "true PU"        
		PUweightAB = lumiWeightsAB.weight( (int) aux.puInfo.truePU ); //use new method with "true PU"        
		pu->Fill(puit->second);
		// PU weight Run2011B
		// PU weight Run2011B
		std::map<int, unsigned int>::const_iterator puit0 =  aux.puInfo.pus.find(0);
		std::map<int, unsigned int>::const_iterator puitm1 = aux.puInfo.pus.find(-1);
		std::map<int, unsigned int>::const_iterator puitp1 = aux.puInfo.pus.find(+1);
		PU0=puit0->second;
		PUp1=puitp1->second;
		PUm1=puitm1->second;
		input3DPU->Fill(PUm1,PU0,PUp1);	  
		//PUweight2011B = lumiWeights2011B.weight3D( puitm1->second, puit0->second,puitp1->second); 
		PUweight1DObs = lumiWeights1DObs.weight( npu); 
      }
      countWithPU->Fill(1,PUweight);
      countWithPUP->Fill(1,PUweightP);
      countWithPUM->Fill(1,PUweightM);
      countWithPUAB->Fill(1,PUweightAB);
      countWithPU2011B->Fill(1,PUweight2011B);
      coutnWithMCProd->Fill(1,weightMCProd);
      countWithPUMCProd->Fill(1,weightMCProd*PUweight);
      //LHE Infos
      fwlite::Handle<LHEEventProduct> evt;

      //	std::cout << "Label for lhe = " << evt.getBranchNameFor(ev,"source") << std::endl;
      if( !((evt.getBranchNameFor(ev,"source")).empty()) )
      {
		evt.getByLabel(ev,"source");
		//std::cout << "LHEEventProduct found!" << std::endl;
		bool lCheck=false;
		bool lbarCheck=false;
		bool vlCheck=false;
		bool vlbarCheck=false;
		//int idl, idlbar;
		lheHT=0.;
		lheNj=0;
		TLorentzVector l,lbar,vl,vlbar,V_tlv;
		const lhef::HEPEUP hepeup_ = evt->hepeup();
		const std::vector<lhef::HEPEUP::FiveVector> pup_ = hepeup_.PUP; // px, py, pz, E, M

		//################# PDF CODE ##################

		//############# PDF reweighting ##################
		//change it if you want to rescale the energy
		double _newECMS = 8000; // energy in GeV
		double _origECMS = 8000; // energy in GeV
		float Q = hepeup_.SCALUP;
		int id1        = hepeup_.IDUP[0];
		double x1      = fabs(hepeup_.PUP[0][2]/(_origECMS/2));
		double x1prime = fabs(hepeup_.PUP[0][2]/(_newECMS/2));
		int id2        = hepeup_.IDUP[1];
		double x2      = fabs(hepeup_.PUP[1][2]/(_origECMS/2));
		double x2prime = fabs(hepeup_.PUP[1][2]/(_newECMS/2));
		//      gluon is 0 in the LHAPDF numberin
		if (id1 == 21) id1 = 0;
		if (id2 == 21) id2 = 0;
		//int _pdfmember = 0; //member = 0 is the pdf central value.
		//madgrpah default pdf
		//LHAPDF::usePDFMember(1,_pdfmember);
		double oldpdf1 = LHAPDF::xfx(1, x1, Q, id1)/x1;
		double oldpdf2 = LHAPDF::xfx(1, x2, Q, id2)/x2;
		double newpdf1=0;
		double newpdf2=0;

		newpdf1 = LHAPDF::xfx(1, x1prime, Q, id1)/x1prime;
		newpdf2 = LHAPDF::xfx(1, x2prime, Q, id2)/x2prime;
		PDFweight[0] = (newpdf1/oldpdf1)*(newpdf2/oldpdf2);


          //new pdf set
        for (unsigned int setpdf=2; setpdf <= pdfNames.size()+1; ++setpdf)
		{
          for(int pdfmember = 0; pdfmember < LHAPDF::numberPDF(setpdf); pdfmember++)
		  {
            LHAPDF::usePDFMember(setpdf, pdfmember); //load the right memeber once again to be sure...
            newpdf1 = LHAPDF::xfx(setpdf, x1prime, Q, id1)/x1prime;
            newpdf2 = LHAPDF::xfx(setpdf, x2prime, Q, id2)/x2prime;
            int norm = (setpdf-2)*LHAPDF::numberPDF(setpdf); // normalisation
            PDFweight[norm+setpdf-1+pdfmember] = (newpdf1/oldpdf1)*(newpdf2/oldpdf2);
          }
        }
        //############### LHE stitching information ##################
	//################# END OF PDF CODE ##################

        for(unsigned int i=0; i<pup_.size(); ++i)
		{
		  int id=hepeup_.IDUP[i]; //pdgId
		  int status = hepeup_.ISTUP[i];
		  int idabs=TMath::Abs(id); 
		  if( status == 1 && ( ( idabs == 21 ) || (idabs > 0 && idabs < 7) ) ){ // gluons and quarks
			lheHT += TMath::Sqrt( TMath::Power(hepeup_.PUP[i][0],2) + TMath::Power(hepeup_.PUP[i][1],2) ); // first entry is px, second py
			lheNj++;
		  }

		  if(id==11){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
		  if(id==-11){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
		  if(id==12){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
		  if(id==-12){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
		  
		  if(id==13){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
		  if(id==-13){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
		  if(id==14){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
		  if(id==-14){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
		  
		  if(id==15){ l.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lCheck=true;}
		  if(id==-15){ lbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); lbarCheck=true;}
		  if(id==16){ vl.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlCheck=true;}
		  if(id==-16){ vlbar.SetPxPyPzE(hepeup_.PUP[i][0],hepeup_.PUP[i][1],hepeup_.PUP[i][2],hepeup_.PUP[i][3]); vlbarCheck=true;}
		  
		}
		if( lCheck && lbarCheck ) V_tlv = l + lbar; // ZtoLL
		if( vlCheck && vlbarCheck ) V_tlv = vl + vlbar; // ZtoNuNu
		if( lCheck && vlbarCheck ) V_tlv = l + vlbar; // WToLNu
		if( lbarCheck && vlCheck ) V_tlv = lbar + vl; // WToLNu       
		lheV_pt = V_tlv.Pt();
      
      } //end LHE information

      //std::cout << "lhe V pt = " << lheV_pt << std::endl;

      //Write event info
      
      // simBHadrons
      const SimBHadronCollection *sbhc = NULL;
      if(isMC_)
      {
        fwlite::Handle<SimBHadronCollection> SBHC;
        SBHC.getByLabel(ev, "bhadrons");
        sbhc = SBHC.product();
      }
	
      // adapted from PhysicsTools/PatUtils/plugins/PATPFJetMETcorrInputProducer.cc
      // To propagate the effect of JEC at Step 2 to METtype1corr
      double METtype1diff_mex=0.0, METtype1diff_mey=0.0, METtype1diff_sumet=0.0, METtype1diff_type1JetPtThreshold=10.0, METtype1diff_skipEMfractionThreshold=0.9;

      const std::vector<VHbbCandidate> * candZ ;
      const std::vector<VHbbCandidate> * candW ;
      VHbbEvent modifiedEvent;
      const VHbbEvent *  iEvent =0;
      
      std::vector<TLorentzVector> jetWithJEC_p4_vec ;
      std::vector<TLorentzVector> jetBestMC_p4_vec ;
      std::vector<int> bestMCid_vec ;
      if(fromCandidate) //fromCandiate = False
      {
		fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandleZ;
		vhbbCandHandleZ.getByLabel(ev,"hbbBestCSVPt20Candidates");
		candZ = vhbbCandHandleZ.product();

		fwlite::Handle< std::vector<VHbbCandidate> > vhbbCandHandle;
		vhbbCandHandle.getByLabel(ev,"hbbHighestPtHiggsPt30Candidates");
		candW = vhbbCandHandle.product();
      }
      else
      {
		candZlocal->clear();
		candWlocal->clear();
		fwlite::Handle< VHbbEvent > vhbbHandle;
		vhbbHandle.getByLabel(ev,"HbbAnalyzerNew");
		modifiedEvent = *vhbbHandle.product();

		if(aux.mcH.size()>0)	// If VHbbEventAuxInfo has H MC Info
		{
		  float pt=aux.mcH[0].p4.Pt();
		  int nMCJets=0;
		  if(aux.mcW.size()>1) pt=aux.mcW[1].p4.Pt();
		  if(aux.mcZ.size()>1) pt=aux.mcZ[1].p4.Pt();
		  for(size_t j=0; j< modifiedEvent.simpleJets2.size() ; j++)
		  {
			//um$(VHbbEvent_HbbAnalyzerNew__VH.obj.simpleJets2.p4.Pt()>20 && abs(VHbbEvent_HbbAnalyzerNew__VH.obj.simpleJets2.p4.Eta() )< 2.5)
			if( modifiedEvent.simpleJets2[j].bestMCp4.Pt()>20 && fabs(modifiedEvent.simpleJets2[j].bestMCp4.Eta())<2.5) nMCJets++;
		  }
		  weightSignalQCD=weightNNLOQCDsignal(pt,nMCJets);
		  weightSignalEWK=weightNLOEWKsignal(pt);
		}
		countWithSignalQCDcorrections->Fill(1,weightSignalQCD);

		// 
		// CODE FOR JEC CORRECTIONS
		// 
	
	
		if(isMC_)
		{
		  iEvent= &modifiedEvent;
	
		  for(size_t j=0; j< modifiedEvent.subJets.size() ; j++)   modifiedEvent.subJets[j] = jec.correct( modifiedEvent.subJets[j],aux.puInfo.rho,true);
		  for(size_t j=0; j< modifiedEvent.filterJets.size() ; j++)   modifiedEvent.filterJets[j] = jec.correct( modifiedEvent.filterJets[j],aux.puInfo.rho,true);
	
			
		  for(size_t j=0; j< modifiedEvent.simpleJets2.size() ; j++)
		  {
			//VHbbEvent::SimpleJet orig=modifiedEvent.simpleJets2[j];
			//VHbbEvent::SimpleJet origRemade = jec.correct( modifiedEvent.simpleJets2[j],aux.puInfo.rho,true,true); // do ref check, can be commented out 
			//VHbbEvent::SimpleJet corr2011 = jec.correctRight( modifiedEvent.simpleJets2[j],aux.puInfo.rho,true,true); // do ref check, can be commented out
			if (modifiedEvent.simpleJets2[j].p4.Pt() > METtype1diff_type1JetPtThreshold && (modifiedEvent.simpleJets2[j].chargedEmEFraction+modifiedEvent.simpleJets2[j].neutralEmEFraction) < METtype1diff_skipEMfractionThreshold)
			{
			  METtype1diff_mex += modifiedEvent.simpleJets2[j].p4.Px();
			  METtype1diff_mey += modifiedEvent.simpleJets2[j].p4.Py();
			  METtype1diff_sumet -= modifiedEvent.simpleJets2[j].p4.Et();
			}
			modifiedEvent.simpleJets2[j] = jec.correct( modifiedEvent.simpleJets2[j],aux.puInfo.rho,true); 
                        jetWithJEC_p4_vec.push_back(modifiedEvent.simpleJets2[j].p4) ;
                        jetBestMC_p4_vec.push_back(modifiedEvent.simpleJets2[j].bestMCp4) ;
                        bestMCid_vec.push_back(modifiedEvent.simpleJets2[j].bestMCid) ;
			//std::cout << "Original " << orig.p4.Pt() << " == " << origRemade.p4.Pt() << " using CHS2011 " << corr2011.p4.Pt() << " final: " << modifiedEvent.simpleJets2[j].p4.Pt() << std::endl;
			if (modifiedEvent.simpleJets2[j].p4.Pt() > METtype1diff_type1JetPtThreshold && (modifiedEvent.simpleJets2[j].chargedEmEFraction+modifiedEvent.simpleJets2[j].neutralEmEFraction) < METtype1diff_skipEMfractionThreshold)
			{
			  METtype1diff_mex -= modifiedEvent.simpleJets2[j].p4.Px();
			  METtype1diff_mey -= modifiedEvent.simpleJets2[j].p4.Py();
			  METtype1diff_sumet += modifiedEvent.simpleJets2[j].p4.Et();
			}
//===resolution smearing MC jet==
			TLorentzVector & p4 = modifiedEvent.simpleJets2[j].p4; 
			TLorentzVector & mcp4 = modifiedEvent.simpleJets2[j].bestMCp4;
                        if ((fabs(p4.Pt() - mcp4.Pt())/ p4.Pt())<0.5)
			{ //Limit the effect to the core 
			  float cor = (p4.Pt()+resolutionBias(fabs(p4.Eta()))*(p4.Pt()-mcp4.Pt()))/p4.Pt();
			  p4.SetPtEtaPhiE(p4.Pt()*cor,p4.Eta(), p4.Phi(), p4.E()*cor);
			}
		  }
		}
		
                else 
		{  // for data
		  //iEvent = vhbbHandle.product();
		  // modify also the real data now to apply JEC 2012
		  iEvent= &modifiedEvent;
		  for(size_t j=0; j< modifiedEvent.subJets.size() ; j++)    modifiedEvent.subJets[j]    = jec.correct( modifiedEvent.subJets[j],   aux.puInfo.rho,false);
		  for(size_t j=0; j< modifiedEvent.filterJets.size() ; j++) modifiedEvent.filterJets[j] = jec.correct( modifiedEvent.filterJets[j],aux.puInfo.rho,false);
		  for(size_t j=0; j< modifiedEvent.simpleJets2.size() ; j++)
		  { 
			//jec.correct( modifiedEvent.simpleJets2[j],aux.puInfo.rho,false,true); // do ref check, can be commented out 
			if(modifiedEvent.simpleJets2[j].p4.Pt() > METtype1diff_type1JetPtThreshold && (modifiedEvent.simpleJets2[j].chargedEmEFraction+modifiedEvent.simpleJets2[j].neutralEmEFraction) < METtype1diff_skipEMfractionThreshold)
			{
			  METtype1diff_mex += modifiedEvent.simpleJets2[j].p4.Px();
			  METtype1diff_mey += modifiedEvent.simpleJets2[j].p4.Py();
			  METtype1diff_sumet -= modifiedEvent.simpleJets2[j].p4.Et();
			}
			modifiedEvent.simpleJets2[j] = jec.correct( modifiedEvent.simpleJets2[j],aux.puInfo.rho,false);                        jetWithJEC_p4_vec.push_back(modifiedEvent.simpleJets2[j].p4) ;
                        jetBestMC_p4_vec.push_back(modifiedEvent.simpleJets2[j].bestMCp4) ;
                        bestMCid_vec.push_back(modifiedEvent.simpleJets2[j].bestMCid) ;

			if (modifiedEvent.simpleJets2[j].p4.Pt() > METtype1diff_type1JetPtThreshold && (modifiedEvent.simpleJets2[j].chargedEmEFraction+modifiedEvent.simpleJets2[j].neutralEmEFraction) < METtype1diff_skipEMfractionThreshold)
			{
			  METtype1diff_mex -= modifiedEvent.simpleJets2[j].p4.Px();
			  METtype1diff_mey -= modifiedEvent.simpleJets2[j].p4.Py();
			  METtype1diff_sumet += modifiedEvent.simpleJets2[j].p4.Et();
			}
		  }
		}  

		// USE HbbCandidateFinderAlgo TO FIND Z,W Candidates
		// note: Modified to return only Z-ll candidates
		algoZ->run(iEvent,*candZlocal,aux);
		algoW->run(iEvent,*candWlocal,aux);

		//if(candZlocal->size() == 0 or candZlocal->at(0).H.jets.size() < 2)  //recover low pt 
		//if(candZlocal->size() == 0 or candZlocal->at(0).H.jets.size() < 1)  //recover low pt 
		//{
		//  candZlocal->clear();
		//  candWlocal->clear();
		//  algoRecoverLowPt->run(iEvent,*candZlocal,aux);
		//  algoRecoverLowPt->run(iEvent,*candWlocal,aux);
		//}
	
		candZ= candZlocal; 
		candW= candWlocal; 
		//for(size_t m=0;m<iEvent->muInfo.size();m++)
		//  if( fabs(iEvent->muInfo[m].p4.Pt()-28.118684) < 0.0001 || fabs(iEvent->muInfo[m].p4.Pt()-34.853199) < 0.0001 )  
		//    std::cout << "FOUND " << iEvent->muInfo[m].p4.Pt() <<  " " << EVENT.event << " " << candW->size() << " " << candZ->size() << std::endl;
      } //end getting Z candiate

      const std::vector<VHbbCandidate> * cand = candZ;
      //std::clog << "Filling tree "<< std::endl;
      bool isW=false;

      // to check how much we gain with jets subtraction 
      tauMinusMode = -99;
      tauPlusMode = -99;
      for(unsigned int j=0; j< aux.mcTau.size();j++)
      {
		for(unsigned int k=0;k< aux.mcTau[j].dauid.size();k++)
		{
		  if (ana.getParameter<bool>("verbose"))
		  {
			std::cout << "(j,k,m_Tau,dauid,pt_dau,momid)=(" << j << "," << k << "," << aux.mcTau[j].p4.M() << "," << aux.mcTau[j].dauid[k] << ","
				  << aux.mcTau[j].dauFourMomentum[k].Pt() << "," << aux.mcTau[j].momid  <<")" << std::endl;
		  }
		  int idd=abs(aux.mcTau[j].dauid[k]);
		  if( idd != 15 && idd!=16  //(idd==11 || idd==13 )
			)
		  {
			if(tauMinusMode==-99 && aux.mcTau[j].charge ==-1) tauMinusMode = idd;
			if(tauPlusMode==-99  && aux.mcTau[j].charge ==+1) tauPlusMode  = idd;
		  }
		}
      } 
	
      genHpt=aux.mcH.size() > 0 ? aux.mcH[0].p4.Pt():-99;

      trigger.setEvent(&ev);
      for(size_t j=0;j < triggers.size();j++) triggerFlags[j]=trigger.accept(triggers[j]);

//==NOTE : bypass the filter basing on the existence of Z candidate (Duong 10-08-2015)
      //if(cand->size() == 0 ) continue;
      if (cand->size() > 0) {
      //if(cand->size() == 0 or cand->at(0).H.jets.size() < 2) continue;
      //std::cout << "cand->size() " << cand->size() << std::endl;
      //std::cout << "cand->at(0).H.jets.size() " << cand->at(0).H.jets.size() << std::endl;
      if(cand->size() > 1 ) std::cout << "MULTIPLE CANDIDATES: " << cand->size() << std::endl;
      if(cand->at(0).candidateType == VHbbCandidate::Wmun || cand->at(0).candidateType == VHbbCandidate::Wen ) { cand=candW; isW=true; }

//==NOTE : bypass the filter basing on the existence of Z candiate (Duong 10-08-2015)
//      if(cand->size() == 0)
//      {
		//std::cout << "W event loss due to tigther cuts" << std::endl;
//        continue;
//      }
      
      // secondary vtxs
      fwlite::Handle<std::vector<reco::Vertex> > SVC;
      SVC.getByLabel(ev,"bcandidates");
      const std::vector<reco::Vertex> svc = *(SVC.product());

	  // Extract leading Candidate
      const VHbbCandidate & vhCand =  cand->at(0);

	  // Log triggers
      patFilters.setEvent(&ev,"VH");
      hbhe                      = patFilters.accept("hbhe");
      ecalFlag 	      		= patFilters.accept("ecalFilter");
      totalKinematics 		= patFilters.accept("totalKinematics");
      cschaloFlag   		= patFilters.accept("cschaloFilter");   
      hcallaserFlag   		= patFilters.accept("hcallaserFilter");   
      trackingfailureFlag 	= patFilters.accept("trackingfailureFilter");   
      eebadscFlag     		= patFilters.accept("eebadscFilter");
   
//===NOTE move trigger flag out of Z candidate finding (Duong Nov. 3, 2015) ===
//      trigger.setEvent(&ev);
//      for(size_t j=0;j < triggers.size();j++) triggerFlags[j]=trigger.accept(triggers[j]);

      eventFlav=0;

      if(aux.mcBbar.size() > 0 || aux.mcB.size() > 0) eventFlav=5;
      else if(aux.mcC.size() > 0) eventFlav=4;

	  // Log VB Information
      Vtype  = vhCand.candidateType;
      V.pt   = vhCand.V.p4.Pt();
      V.eta  = vhCand.V.p4.Eta();
      V.phi  = vhCand.V.p4.Phi();
      VMt    = vhCand.Mt();
      V.mass = vhCand.V.p4.M();
      if(isW) V.mass = vhCand.Mt();
      
      // CHANGED: Added loop to store allJets, allMuons, allElectrons info
      allJets.reset();
      nallJets = vhCand.allJets.size();
      int indTmp(-1) ;
      if (vhCand.allJets.size() != vhCand.allJets_oriInd.size()) cout << "\n allJets != allJets_oriInd " ;
      for(int j=0; j<nallJets && j < MAXJ; j++) {
        indTmp = vhCand.allJets_oriInd[j] ;
        allJets.set(vhCand.allJets[j],j,jetWithJEC_p4_vec[indTmp], jetBestMC_p4_vec[indTmp], iEvent->simpleJets2[indTmp].bestMCid);
      }
 
      allMuons.reset();
      nallMuons = vhCand.V.muons.size();
      for(int j=0; j<nallMuons && j < MAXL; j++) allMuons.set(vhCand.V.muons[j],j,13,aux, j);

      allElectrons.reset();
      nallElectrons = vhCand.V.electrons.size(); //electrons = all electron inside HbbAnalyzerNew of step 1 
      for(int j=0; j<nallElectrons && j < MAXL; j++) allElectrons.set(vhCand.V.electrons[j],j,11,aux, j);

	  // Log Higgs Information
      if(vhCand.H.HiggsFlag) H.HiggsFlag=1;
	  else H.HiggsFlag=0;
      if(vhCand.H.HiggsFlag)
	  {
		H.mass = vhCand.H.p4.M();
		H.pt   = vhCand.H.p4.Pt();
		H.eta  = vhCand.H.p4.Eta();
		H.phi  = vhCand.H.p4.Phi();
      }

	  // Log FatHiggs Information
      nfathFilterJets=0;
      if(vhCand.FatH.FatHiggsFlag) FatH.FatHiggsFlag =1;
	  else FatH.FatHiggsFlag=0;
      fathFilterJets.reset();
      
//==fill additional fat jets===
      aJetsFat.reset();
      if(vhCand.FatH.FatHiggsFlag)
	  { 
		FatH.mass= vhCand.FatH.p4.M(); 
		FatH.pt = vhCand.FatH.p4.Pt();
		if(FatH.pt!=0) FatH.eta = vhCand.FatH.p4.Eta();
		else		   FatH.eta = -99.;
		FatH.phi = vhCand.FatH.p4.Phi();

		//if(vhCand.FatH.FatHiggsFlag)  vhCand.FatH.subjetsSize; 
        nfathFilterJets=vhCand.FatH.subjetsSize;  
		for( int j=0; j < nfathFilterJets; j++ ) fathFilterJets.set(vhCand.FatH.jets[j],j, vhCand.FatH.jets[j].p4, vhCand.FatH.jets[j].p4, -99); //==TODO : need to use correct jetWithJEC and jetWithMC and bestMCid

		if(nfathFilterJets==2)
		{
		  FatH.filteredmass=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).M();
		  FatH.filteredpt=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).Pt();
		  FatH.filteredeta=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).Eta();
		  FatH.filteredphi=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4).Phi();
		  fathFilterJets.cosTheta[0]=  vhCand.FatH.helicities[0];
		  fathFilterJets.cosTheta[1]=  vhCand.FatH.helicities[1];
		}
		else if(nfathFilterJets==3)
		{
		  FatH.filteredmass=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).M();
		  FatH.filteredpt=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).Pt();
		  FatH.filteredeta=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).Eta();
		  FatH.filteredphi=(vhCand.FatH.jets[0].p4+vhCand.FatH.jets[1].p4+vhCand.FatH.jets[2].p4).Phi();
		  fathFilterJets.cosTheta[0]=  vhCand.FatH.helicities[0];
		  fathFilterJets.cosTheta[1]=  vhCand.FatH.helicities[1];
		  fathFilterJets.cosTheta[2]=  vhCand.FatH.helicities[2];
		}

		naJetsFat=vhCand.additionalJetsFat.size();
		for( int j=0; j < naJetsFat && j < MAXJ; j++ ) aJetsFat.set(vhCand.additionalJetsFat[j],j, vhCand.additionalJetsFat[j].p4, vhCand.additionalJetsFat[j].p4, -99); //==TODO : need use correct jetWithJEC and jetBestMC and bestMCid

      } // FatHiggsFlag

// Log Higgs Candidate Info

//==fill higgs jet and additional jets==
      hJets.reset();
//==NOTE : move aJets out of VH candidate==
//      aJets.reset();
      if(vhCand.H.HiggsFlag)
	  {
		nhJets=2;
		hJets.set(vhCand.H.jets[0],0,vhCand.H.jets[0].p4,vhCand.H.jets[0].p4,-99) ; //==TODO : need to use correct jetWithJEC and jetBestMC
		hJets.set(vhCand.H.jets[1],1,vhCand.H.jets[1].p4,vhCand.H.jets[1].p4,-99);

		VtypeWithTau=vhCand.candidateTypeWithTau;
//		aJets.reset();
  
		naJets=vhCand.additionalJets.size();
		numBJets=0;
		if(vhCand.H.jets[0].csv> btagThr) numBJets++;
		if(vhCand.H.jets[1].csv> btagThr) numBJets++;
		for( int j=0; j < naJets && j < MAXJ; j++ ) 
	    {
//		  aJets.set(vhCand.additionalJets[j],j);
		  if(vhCand.additionalJets[j].csv> btagThr) numBJets++;
//		  if (VtypeWithTau==VHbbCandidate::Wtaun) aJets.selectedTauDR[j] = vhCand.VTau.taus[0].p4.DeltaR(vhCand.additionalJets[j].p4);  
		}

		numJets = vhCand.additionalJets.size()+2;
		H.dR = deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi());
		H.dPhi = deltaPhi(vhCand.H.jets[0].p4.Phi(),vhCand.H.jets[1].p4.Phi());
		H.dEta= TMath::Abs( vhCand.H.jets[0].p4.Eta() - vhCand.H.jets[1].p4.Eta() );
		HVdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.p4.Phi()) ) ;
		HVMass = (vhCand.H.p4 + vhCand.V.p4).M() ;
		HMETdPhi = fabs( deltaPhi(vhCand.H.p4.Phi(),vhCand.V.mets.at(0).p4.Phi()) ) ;
		//eltaPullAngle = vhCand.H.deltaTheta;
  
		deltaPullAngle  =  VHbbCandidateTools::getDeltaTheta(vhCand.H.jets[0],vhCand.H.jets[1]);
		deltaPullAngle2 =  VHbbCandidateTools::getDeltaTheta(vhCand.H.jets[1],vhCand.H.jets[0]);
		hJets.cosTheta[0]=  vhCand.H.helicities[0];
		hJets.cosTheta[1]=  vhCand.H.helicities[1];
      } // Higgs Flag

      if(VtypeWithTau==VHbbCandidate::Wtaun)
	  {
		//VTau.mass = vhCand.VTau.p4.M();
		VTau.mass = vhCand.MtTau();
		VTau.pt = vhCand.VTau.p4.Pt();
		VTau.eta = vhCand.VTau.p4.Eta();
		VTau.phi = vhCand.VTau.p4.Phi();
	  }

	  // METInfo calomet;  METInfo tcmet;  METInfo pfmet;  METInfo mht;  METInfo metNoPU
	  //
	  //====pfMETtype1corr with some filter but currently don't apply see HbbCandidateFinderAlgo.cc findMET()
      MET.et = vhCand.V.mets.at(0).p4.Pt();
      MET.phi = vhCand.V.mets.at(0).p4.Phi();
      MET.sumet = vhCand.V.mets.at(0).sumEt;
      MET.sig = vhCand.V.mets.at(0).metSig;
      
      //====raw MET===
      pfMET.et    = iEvent->pfmet.p4.Pt();
      pfMET.phi   = iEvent->pfmet.p4.Phi();
      pfMET.sumet = iEvent->pfmet.sumEt;
      pfMET.sig   = iEvent->pfmet.metSig;

      fakeMET.sumet = 0;
      fakeMET.sig = 0;
      fakeMET.et = 0;
      fakeMET.phi = 0;

      if( Vtype == VHbbCandidate::Zmumu)
	  {
		TVector3 mu1 = vhCand.V.muons[0].p4.Vect();
		TVector3 mu2 = vhCand.V.muons[vhCand.V.secondLepton].p4.Vect();
		// Not needed with PFMET
		//mu1.SetMag( mu1.Mag() - vhCand.V.muons[0].emEnergy - vhCand.V.muons[0].hadEnergy);
		//mu2.SetMag( mu2.Mag() - vhCand.V.muons[1].emEnergy - vhCand.V.muons[1].hadEnergy);
		TVector3 sum = vhCand.V.mets.at(0).p4.Vect() + mu1 + mu2;
		fakeMET.et = sum.Pt();
		fakeMET.phi = sum.Phi();
		fakeMET.sumet = vhCand.V.mets.at(0).sumEt - mu1.Pt() - mu2.Pt();
	  }

      METnoPU.et = iEvent->metNoPU.p4.Pt();
      METnoPU.phi = iEvent->metNoPU.p4.Phi();
      METnoPU.sumet = iEvent->metNoPU.sumEt;
      METnoPU.sig = iEvent->metNoPU.metSig;
      METnoPUCh.et = iEvent->metCh.p4.Pt();
      METnoPUCh.phi = iEvent->metCh.p4.Phi();
      METnoPUCh.sumet = iEvent->metCh.sumEt;
      METnoPUCh.sig = iEvent->metCh.metSig;

      METnoPUCh.et = iEvent->metCh.p4.Pt();
      METnoPUCh.phi = iEvent->metCh.p4.Phi();
      METnoPUCh.sumet = iEvent->metCh.sumEt;
      METnoPUCh.sig = iEvent->metCh.metSig;

      METtype1corr.et = iEvent->pfmetType1corr.p4.Pt();
      METtype1corr.phi = iEvent->pfmetType1corr.p4.Phi();
      METtype1corr.sumet = iEvent->pfmetType1corr.sumEt;
      METtype1corr.sig = iEvent->pfmetType1corr.metSig;

      METtype1p2corr.et = iEvent->pfmetType1p2corr.p4.Pt();
      METtype1p2corr.phi = iEvent->pfmetType1p2corr.p4.Phi();
      METtype1p2corr.sumet = iEvent->pfmetType1p2corr.sumEt;
      METtype1p2corr.sig = iEvent->pfmetType1p2corr.metSig;

      METnoPUtype1corr.et = iEvent->pfmetNoPUType1corr.p4.Pt();
      METnoPUtype1corr.phi = iEvent->pfmetNoPUType1corr.p4.Phi();
      METnoPUtype1corr.sumet = iEvent->pfmetNoPUType1corr.sumEt;
      METnoPUtype1corr.sig = iEvent->pfmetNoPUType1corr.metSig;

      METnoPUtype1p2corr.et = iEvent->pfmetNoPUType1p2corr.p4.Pt();
      METnoPUtype1p2corr.phi = iEvent->pfmetNoPUType1p2corr.p4.Phi();
      METnoPUtype1p2corr.sumet = iEvent->pfmetNoPUType1p2corr.sumEt;
      METnoPUtype1p2corr.sig = iEvent->pfmetNoPUType1p2corr.metSig;
  
      TVector2 METtype1diff_vec2(METtype1diff_mex, METtype1diff_mey);
      METtype1diff.et = METtype1diff_vec2.Mod();
      METtype1diff.phi = METtype1diff_vec2.Phi();
      METtype1diff.sumet = METtype1diff_sumet;
      METnoPUtype1p2corr.sig = 0;

      //	std::cout << " iEvent->metUncInfo.size() " << iEvent->metUncInfo.size() << std::endl;
      for(size_t m=0;m<iEvent->metUncInfo.size();m++)
	  { 
		metUnc.set(iEvent->metUncInfo[m], m );
		//std::cout << "metUncInfo[" << m <<" ].et = " << metUnc.et[m] << std::endl; 
	  }

      rho = aux.puInfo.rho;
      rho25 = aux.puInfo.rho25;
      rhoN = aux.puInfo.rhoNeutral;
      nPVs=aux.pvInfo.nVertices; 

      if(!fromCandidate)
	  {
		MHT.mht = iEvent->mht.p4.Pt(); 
		MHT.phi = iEvent->mht.p4.Phi(); 
		MHT.ht = iEvent->mht.sumEt; 
		MHT.sig = iEvent->mht.metSig; 
      }

      /////////
      // track sharing flags:
      ////////
      TkSharing.HiggsCSVtkSharing = TkSharing.HiggsIPtkSharing = TkSharing.HiggsSVtkSharing = TkSharing.FatHiggsCSVtkSharing = TkSharing.FatHiggsIPtkSharing = TkSharing.FatHiggsSVtkSharing = false;

      // csv tracks
      if(vhCand.H.HiggsFlag)
	  {
		if (vhCand.H.jets[0].csvNTracks > 0 && vhCand.H.jets[1].csvNTracks > 0){
		  for (int t=0;t!=vhCand.H.jets[0].csvNTracks;t++){
			for (int ti=0;ti!=vhCand.H.jets[1].csvNTracks;ti++){
			  if ((int)vhCand.H.jets[0].csvTrackIds[t] == (int)vhCand.H.jets[1].csvTrackIds[ti]){
				TkSharing.HiggsCSVtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet
  
		// ip tracks
		if (vhCand.H.jets[0].btagNTracks > 0 && vhCand.H.jets[1].btagNTracks > 0){
		  for (int t=0;t!=vhCand.H.jets[0].btagNTracks;t++){
			for (int ti=0;ti!=vhCand.H.jets[1].btagNTracks;ti++){
			  if ((int)vhCand.H.jets[0].btagTrackIds[t] == (int)vhCand.H.jets[1].btagTrackIds[ti]){
				TkSharing.HiggsIPtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet
		
		// sv tracks
		if (vhCand.H.jets[0].vtxNTracks > 0 && vhCand.H.jets[1].vtxNTracks > 0){
		  for (int t=0;t!=vhCand.H.jets[0].vtxNTracks;t++){
			for (int ti=0;ti!=vhCand.H.jets[1].vtxNTracks;ti++){
			  if ((int)vhCand.H.jets[0].vtxTrackIds[t] == (int)vhCand.H.jets[1].vtxTrackIds[ti]){
				TkSharing.HiggsSVtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet

	  } // Di-jet Higgs Flag
   
	  // tracksharing for Filtered jets:
	  if(vhCand.FatH.FatHiggsFlag && nfathFilterJets > 1)
	  {
		// csv tracks
		if (vhCand.FatH.jets[0].csvNTracks > 0 && vhCand.FatH.jets[1].csvNTracks > 0){
		  for (int t=0;t!=vhCand.FatH.jets[0].csvNTracks;t++){
			for (int ti=0;ti!=vhCand.FatH.jets[1].csvNTracks;ti++){
			  if ((int)vhCand.FatH.jets[0].csvTrackIds[t] == (int)vhCand.FatH.jets[1].csvTrackIds[ti]){
				TkSharing.FatHiggsCSVtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet
		
		// ip tracks
		if (vhCand.FatH.jets[0].btagNTracks > 0 && vhCand.FatH.jets[1].btagNTracks > 0){
		  for (int t=0;t!=vhCand.FatH.jets[0].btagNTracks;t++){
			for (int ti=0;ti!=vhCand.FatH.jets[1].btagNTracks;ti++){
			  if ((int)vhCand.FatH.jets[0].btagTrackIds[t] == (int)vhCand.FatH.jets[1].btagTrackIds[ti]){
				TkSharing.FatHiggsIPtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet
		
		// sv tracks
		if (vhCand.FatH.jets[0].vtxNTracks > 0 && vhCand.FatH.jets[1].vtxNTracks > 0){
		  for (int t=0;t!=vhCand.FatH.jets[0].vtxNTracks;t++){
			for (int ti=0;ti!=vhCand.FatH.jets[1].vtxNTracks;ti++){
			  if ((int)vhCand.FatH.jets[0].vtxTrackIds[t] == (int)vhCand.FatH.jets[1].vtxTrackIds[ti]){
				TkSharing.FatHiggsSVtkSharing = true;
			  }// same trackID
			}// loop tracks in second hjet
		  }// loop tracks in first hjet
		}// if tracks in jet

	  }// fatH
      
      //Secondary Vertices
      IVF.reset();
      nSvs = svc.size();
      const TVector3 recoPv = aux.pvInfo.firstPVInPT2;
      const math::XYZPoint myPv(recoPv);

      // FAKE ERROR MATRIX
	  // look here for Matrix filling info http://project-mathlibs.web.cern.ch/project-mathlibs/sw/html/SMatrixDoc.html
	  //std::vector<double> fillMatrix(6);
	  //for (int i = 0; i<6; ++i) fillMatrix[i] = 0.;
	  //fillMatrix[0] = TMath::Power(0.002,2);
	  //fillMatrix[2] = TMath::Power(0.002,2);
	  //fillMatrix[5] = TMath::Power(0.002,2);
	  //const ROOT::Math::SMatrix<double, 3, 3, ROOT::Math::MatRepSym<double, 3> > myFakeMatrixError(fillMatrix.begin(),fillMatrix.end());
	  //const reco::Vertex recoVtxPv(myPv, myFakeMatrixError);

	  // REAL ERROR MATRIX
      const reco::Vertex recoVtxPv(myPv, aux.pvInfo.efirstPVInPT2);
      for( int j=0; j < nSvs && j < MAXB; ++j )
	  {
		const GlobalVector flightDir = flightDirection(recoPv,svc[j]);
		reco::SecondaryVertex recoSv(recoVtxPv, svc[j], flightDir ,true);
		IVF.set( recoSv, recoPv ,j);
      }
	  
      if(nSvs > 1)
	  {
		TLorentzVector BCands_H1, BCands_H2, BCands_H;
		BCands_H1.SetPtEtaPhiM(IVF.pt[0], IVF.eta[0], IVF.phi[0], IVF.massBcand[0]);
		BCands_H2.SetPtEtaPhiM(IVF.pt[1], IVF.eta[1], IVF.phi[1], IVF.massBcand[1]);
		BCands_H = BCands_H1 + BCands_H2;
		SVH.dR = deltaR(IVF.eta[0], IVF.phi[0], IVF.eta[1], IVF.phi[1] );
		SVH.dPhi = deltaPhi(IVF.phi[0], IVF.phi[1] );
		SVH.dEta = TMath::Abs(IVF.eta[0] - IVF.eta[1] );
		SVH.mass = BCands_H.M();
		SVH.pt = BCands_H.Pt();
		SVH.eta = BCands_H.Eta();
		SVH.phi = BCands_H.Phi();
      }

      //SimBhadron
      SimBs.reset();
      if(isMC_)
	  {
		nSimBs = sbhc->size();
		for( int j=0; j < nSimBs && j < MAXB; ++j ) SimBs.set( sbhc->at(j), j);
		if(nSimBs > 1)
		{
		  TLorentzVector SimBs_H1, SimBs_H2, SimBs_H;
		  SimBs_H1.SetPtEtaPhiM(SimBs.pt[0], SimBs.eta[0], SimBs.phi[0], SimBs.mass[0]);
		  SimBs_H2.SetPtEtaPhiM(SimBs.pt[1], SimBs.eta[1], SimBs.phi[1], SimBs.mass[1]);
		  SimBs_H = SimBs_H1 + SimBs_H2;
		  SimBsH.dR = deltaR(SimBs.eta[0], SimBs.phi[0], SimBs.eta[1], SimBs.phi[1] );
		  SimBsH.dPhi = deltaPhi(SimBs.phi[0], SimBs.phi[1] );
		  SimBsH.dEta = TMath::Abs(SimBs.eta[0] - SimBs.eta[1] );
		  SimBsH.mass = SimBs_H.M();
		  SimBsH.pt = SimBs_H.Pt();
		  SimBsH.eta = SimBs_H.Eta();
		  SimBsH.phi = SimBs_H.Phi();
		}
      }

      //Loop on jets
      double maxBtag=-99999;
      minDeltaPhijetMET = 999;
      TLorentzVector bJet;
      std::vector<std::vector<BTagWeight::JetInfo> > btagJetInfos;
      std::vector<float> jet10eta;
      std::vector<float> jet10pt;
      std::vector<float> jet30eta;
      std::vector<float> jet30pt;
      if(fromCandidate)
	  {
		//Loop on Higgs Jets
		for(unsigned int j=0; j < vhCand.H.jets.size(); j++ )
		{
		  if (vhCand.H.jets[j].csv > maxBtag)
		  {
			bJet=vhCand.H.jets[j].p4;
			maxBtag =vhCand.H.jets[j].csv;
		  }
		  if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())) < minDeltaPhijetMET) 
		  {
			minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.H.jets[j].p4.Phi())); 
			jetPt_minDeltaPhijetMET=vhCand.H.jets[j].p4.Pt();
		  }
		  btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
		}
		//Loop on Additional Jets
		for(unsigned int j=0; j < vhCand.additionalJets.size(); j++ )
		{
		  if (vhCand.additionalJets[j].csv > maxBtag)
		  {
			bJet=vhCand.additionalJets[j].p4;
			maxBtag =vhCand.additionalJets[j].csv;
		  }
		  //if (fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi())) < minDeltaPhijetMET) 
		  //{
		  //  minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), vhCand.additionalJets[j].p4.Phi()));
		  //  jetPt_minDeltaPhijetMET=vhCand.additionalJets[j].p4.Pt();
		  //}
		  if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  )  // btag SF computed using only H-jets if best-H made with dijetPt rather than best CSV
			if(vhCand.additionalJets[j].p4.Pt() > 20) btagJetInfos.push_back(btagEff.jetInfo(vhCand.additionalJets[j]));
		}
	  } 
      else
	  {
		//Loop on all jets
		for(unsigned int j=0; j < iEvent->simpleJets2.size(); j++ )
		{
		  if (iEvent->simpleJets2[j].csv > maxBtag)
		  {
			bJet=iEvent->simpleJets2[j].p4;
			maxBtag =iEvent->simpleJets2[j].csv;
		  }
		  if (iEvent->simpleJets2[j].p4.Pt() > 20 &&  fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5 && fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi())) < minDeltaPhijetMET)
		  {
			minDeltaPhijetMET=fabs(deltaPhi( vhCand.V.mets.at(0).p4.Phi(), iEvent->simpleJets2[j].p4.Phi()));
			jetPt_minDeltaPhijetMET=iEvent->simpleJets2[j].p4.Pt();
		  }
		  if(iEvent->simpleJets2[j].p4.Pt() > 10)
		  {
			jet10eta.push_back(iEvent->simpleJets2[j].p4.Eta());
			jet10pt.push_back(iEvent->simpleJets2[j].p4.Pt());
		  }
		  if(iEvent->simpleJets2[j].p4.Pt() > 30)
		  {
			jet30eta.push_back(iEvent->simpleJets2[j].p4.Eta());
			jet30pt.push_back(iEvent->simpleJets2[j].p4.Pt());
		  }
   
		  //For events made with highest CSV, all jets in the event should be taken into account for "tagging" SF (anti tagging is a mess)
		  // because for example a light jet not used for the Higgs can have in reality a higher CSV due to SF > 1 and become a higgs jet
		  if( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ )  ) 
		    if(iEvent->simpleJets2[j].p4.Pt() > 20 && fabs(iEvent->simpleJets2[j].p4.Eta()) < 2.5)
			  btagJetInfos.push_back(btagEff.jetInfo(iEvent->simpleJets2[j]));
		}

		//Loop on Higgs jets
		if(vhCand.H.HiggsFlag)
		{
		  for(unsigned int j=0; j < vhCand.H.jets.size(); j++ )
			//if we use the highest pt pair, only the two higgs jet should be used to compute the SF because the other jets are excluded 
			// by a criteria (pt of the dijet) that is not btag SF dependent 
			if(!( ( isW && ! useHighestPtHiggsW ) ||  ( ! isW && ! useHighestPtHiggsZ ) ))
			  btagJetInfos.push_back(btagEff.jetInfo(vhCand.H.jets[j]));
		} // HiggsFlag
	  }

      vLeptons.reset();
      weightTrig = 1.; // better to default to 1 
      weightTrigMay = -1.;
      weightTrigV4 = -1.; 
      weightTrigOrMu30 = 1.;
      TLorentzVector leptonForTop;
//      size_t firstAddMu = 0 ;
//      size_t firstAddEle = 0 ;
//      size_t firstAddTau=0;

      if(Vtype == VHbbCandidate::Zmumu )
	  {
		vLeptons.set(vhCand.V.muons[0],0,13,aux,1); 
		vLeptons.set(vhCand.V.muons[vhCand.V.secondLepton],1,13,aux,2);
		float cweightID   = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleMuID(vLeptons.pt[1],vLeptons.eta[1]) ;
		float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
		float weightTrig2 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[1],vLeptons.eta[1]);
		float cweightTrig = weightTrig1 + weightTrig2 - weightTrig1*weightTrig2;
		//2011
		weightTrig = cweightID * cweightTrig;
	
	
		weightTrig2012DiMuon = triggerWeight.doubleMuon2012A(vLeptons.pt[0],vLeptons.eta[0],vLeptons.pt[1],vLeptons.eta[1]);
		float weightTrig2012SingleMuonMu1 = triggerWeight.singleMuon2012A(vLeptons.pt[0],vLeptons.eta[0]);         
		float weightTrig2012SingleMuonMu2 = triggerWeight.singleMuon2012A(vLeptons.pt[1],vLeptons.eta[1]);         
		weightTrig2012SingleMuon = weightTrig2012SingleMuonMu1+weightTrig2012SingleMuonMu2-weightTrig2012SingleMuonMu1*weightTrig2012SingleMuonMu2;
		weightTrig2012= weightTrig2012SingleMuon * triggerWeight.muId2012A(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.muId2012A(vLeptons.pt[1],vLeptons.eta[1]) ;
	
		weightTrig2012ABDiMuon = triggerWeightAB.doubleMuon2012A(vLeptons.pt[0],vLeptons.eta[0],vLeptons.pt[1],vLeptons.eta[1]);
		float weightTrig2012ABSingleMuonMu1 = triggerWeightAB.singleMuon2012A(vLeptons.pt[0],vLeptons.eta[0]);         
		float weightTrig2012ABSingleMuonMu2 = triggerWeightAB.singleMuon2012A(vLeptons.pt[1],vLeptons.eta[1]);         
		weightTrig2012ABSingleMuon = weightTrig2012ABSingleMuonMu1+weightTrig2012ABSingleMuonMu2-weightTrig2012ABSingleMuonMu1*weightTrig2012ABSingleMuonMu2;
		weightTrig2012AB = weightTrig2012ABSingleMuon * triggerWeightAB.muId2012A(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeightAB.muId2012A(vLeptons.pt[1],vLeptons.eta[1]) ;
	
		weightTrig2012A=weightTrig2012;
		weightTrig2012ADiMuon=weightTrig2012DiMuon;
		weightTrig2012ASingleMuon=weightTrig2012SingleMuon;
	
		nvlep=2;
//		firstAddMu=2;
      }

      if( Vtype == VHbbCandidate::Zee )
	  {
		vLeptons.set(vhCand.V.electrons[0],0,11,aux, 1);
		//std::cout << vhCand.V.secondLepton << std::endl;
		vLeptons.set(vhCand.V.electrons[vhCand.V.secondLepton],1,11,aux, 2);
		nvlep=2;
//		firstAddEle=2;
		std::vector<float> pt,eta;
		pt.push_back(vLeptons.pt[0]); eta.push_back(vLeptons.eta[0]);
		pt.push_back(vLeptons.pt[1]); eta.push_back(vLeptons.eta[1]);
		weightEleRecoAndId =  triggerWeight.scaleID95Ele(vLeptons.pt[0],vLeptons.eta[0])
							* triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0])
							* triggerWeight.scaleID95Ele(vLeptons.pt[1],vLeptons.eta[1])
							* triggerWeight.scaleRecoEle(vLeptons.pt[1],vLeptons.eta[1]);
		weightEleTrigElePart 	= triggerWeight.scaleDoubleEle17Ele8(pt,eta); 
		weightEleTrigEleAugPart = triggerWeight.scaleDoubleEle17Ele8Aug(pt,eta); 
		weightTrig = (weightEleTrigElePart*1.14+weightEleTrigEleAugPart*0.98 )/2.12 * weightEleRecoAndId;
	
		weightTrig2012ADiEle = triggerWeight.doubleEle2012A(vLeptons.pt[0],vLeptons.eta[0],vLeptons.pt[1],vLeptons.eta[1]);
		float weightTrig2012ASingleEle1 = triggerWeight.singleEle2012Awp95(vLeptons.pt[0],vLeptons.eta[0]);
		float weightTrig2012ASingleEle2 = triggerWeight.singleEle2012Awp95(vLeptons.pt[1],vLeptons.eta[1]);
		weightTrig2012ASingleEle = weightTrig2012ASingleEle1+weightTrig2012ASingleEle2-weightTrig2012ASingleEle1*weightTrig2012ASingleEle2;
		weightTrig2012A = weightTrig2012ADiEle * triggerWeight.eleId2012A(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.eleId2012A(vLeptons.pt[1],vLeptons.eta[1]) ; 
			   
		weightTrig2012DiEle = triggerWeight.doubleEle2012A(vLeptons.pt[0],vLeptons.eta[0],vLeptons.pt[1],vLeptons.eta[1]);
		float weightTrig2012SingleEle1 = triggerWeight.singleEle2012Awp95(vLeptons.pt[0],vLeptons.eta[0]);
		float weightTrig2012SingleEle2 = triggerWeight.singleEle2012Awp95(vLeptons.pt[1],vLeptons.eta[1]);
		weightTrig2012SingleEle = weightTrig2012SingleEle1+weightTrig2012SingleEle2-weightTrig2012SingleEle1*weightTrig2012SingleEle2;
		weightTrig2012 = weightTrig2012DiEle * triggerWeight.eleId2012A(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.eleId2012A(vLeptons.pt[1],vLeptons.eta[1]) ; 
	
		weightTrig2012ABDiEle = triggerWeightAB.doubleEle2012A(vLeptons.pt[0],vLeptons.eta[0],vLeptons.pt[1],vLeptons.eta[1]);
		float weightTrig2012ABSingleEle1 = triggerWeightAB.singleEle2012Awp95(vLeptons.pt[0],vLeptons.eta[0]);
		float weightTrig2012ABSingleEle2 = triggerWeightAB.singleEle2012Awp95(vLeptons.pt[1],vLeptons.eta[1]);
		weightTrig2012ABSingleEle = weightTrig2012ABSingleEle1+weightTrig2012ABSingleEle2-weightTrig2012ABSingleEle1*weightTrig2012ABSingleEle2;
		weightTrig2012AB = weightTrig2012ABDiEle * triggerWeightAB.eleId2012A(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeightAB.eleId2012A(vLeptons.pt[1],vLeptons.eta[1]) ; 
      }

      if(Vtype == VHbbCandidate::Wmun )
	  {
		leptonForTop=vhCand.V.muons[0].p4;
		vLeptons.set(vhCand.V.muons[0],0,13,aux, 1); 
		float cweightID = triggerWeight.scaleMuID(vLeptons.pt[0],vLeptons.eta[0]);
		float weightTrig1 = triggerWeight.scaleMuIsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
		float cweightTrig = weightTrig1;
		weightTrig = cweightID * cweightTrig;
		float weightTrig1OrMu30 = triggerWeight.scaleMuOr30IsoHLT(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrigOrMu30 = cweightID*weightTrig1OrMu30;
		weightTrig2012ASingleMuon = triggerWeight.singleMuon2012A(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrig2012AMuonPlusWCandPt = weightTrig2012ASingleMuon + 
						 triggerWeight.muPlusWCandPt2012A_legMu(vLeptons.pt[0],vLeptons.eta[0])
						 *triggerWeight.muPlusWCandPt2012A_legW(vhCand.V.p4.Pt(),0);
		weightTrig2012A =  weightTrig2012ASingleMuon *  triggerWeight.muId2012A(vLeptons.pt[0],vLeptons.eta[0]) ;          
		weightTrig2012SingleMuon = weightTrig2012ASingleMuon;
		weightTrig2012MuonPlusWCandPt = weightTrig2012AMuonPlusWCandPt;
		weightTrig2012 = weightTrig2012A;

		weightTrig2012ABSingleMuon = triggerWeightAB.singleMuon2012A(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrig2012ABMuonPlusWCandPt = weightTrig2012ABSingleMuon +
						  triggerWeightAB.muPlusWCandPt2012A_legMu(vLeptons.pt[0],vLeptons.eta[0])
						 *triggerWeightAB.muPlusWCandPt2012A_legW(vhCand.V.p4.Pt(),0);
		weightTrig2012AB =  weightTrig2012ABSingleMuon *  triggerWeightAB.muId2012A(vLeptons.pt[0],vLeptons.eta[0]) ;

		nvlep=1;
//		firstAddMu=1;
      }

      if( Vtype == VHbbCandidate::Wen )
	  {
		leptonForTop=vhCand.V.electrons[0].p4;
		vLeptons.set(vhCand.V.electrons[0],0,11,aux, 1);
		nvlep=1;
//		firstAddEle=1;
		weightTrigMay = triggerWeight.scaleSingleEleMay(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrigV4 = triggerWeight.scaleSingleEleV4(vLeptons.pt[0],vLeptons.eta[0]);
		weightEleRecoAndId=triggerWeight.scaleID80Ele(vLeptons.pt[0],vLeptons.eta[0]) * triggerWeight.scaleRecoEle(vLeptons.pt[0],vLeptons.eta[0]);
		weightEleTrigJetMETPart=triggerWeight.scaleJet30Jet25(jet30pt,jet30eta)*triggerWeight.scalePFMHTEle(MET.et);
		weightEleTrigElePart= weightTrigV4; //this is for debugging only, checking only the V4 part
	
		weightTrigMay*=weightEleRecoAndId;
		weightTrigV4*=weightEleRecoAndId;
		weightTrigV4*=weightEleTrigJetMETPart;
		//weightTrig = weightTrigMay * 0.187 + weightTrigV4 * (1.-0.187); //FIXME: use proper lumi if we reload 2.fb
		weightTrig = (weightTrigMay * 0.215 + weightTrigV4 * 1.915)/ 2.13; //FIXME: use proper lumi if we reload 2.fb

		weightTrig2012ASingleEle = triggerWeight.singleEle2012Awp80(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrig2012A =  weightTrig2012ASingleEle * triggerWeight.eleId2012Awp80(vLeptons.pt[0],vLeptons.eta[0]) ;
	
		weightTrig2012SingleEle = weightTrig2012ASingleEle;
		weightTrig2012 = weightTrig2012A;
		  
		weightTrig2012ABSingleEle = triggerWeightAB.singleEle2012Awp80(vLeptons.pt[0],vLeptons.eta[0]);
		weightTrig2012AB =  weightTrig2012ABSingleEle * triggerWeightAB.eleId2012Awp80(vLeptons.pt[0],vLeptons.eta[0]) ;
      }

      vLeptonsTaus.reset();

      if(VtypeWithTau==VHbbCandidate::Wtaun)
	  {

		if ( vhCand.VTau.taus.size() > 0 )
		{
		  vLeptonsTaus.set(vhCand.VTau.taus[0],0,15,aux,1); 
		  //cout << vLeptonsTaus.pt[0] << ", " << vLeptonsTaus.decayModeFinding[0] << endl;
		  nvlepTau=1;
		  //firstAddTau=1;
		}
      }
      
      if(isMC_)
	  {
        weightTrigMET80 =  triggerWeight.scaleMET80(MET.et);
        weightTrigMET100 =  triggerWeight.scaleMET80(MET.et);
        weightTrig2CJet20 = triggerWeight.scale2CentralJet( jet10pt, jet10eta);
        weightTrigMET150 = triggerWeight.scaleMET150(MET.et);
        weightTrigMET802CJet= weightTrigMET80 * weightTrig2CJet20;
        weightTrigMET1002CJet= weightTrigMET100 * weightTrig2CJet20;
        weightTrig2012DiJet30MHT80=triggerWeight.scaleDiJet30MHT80_2012A(vhCand.V.mets.at(0).p4.Pt()); 
        weightTrig2012PFMET150=triggerWeight.scalePFMET150_2012AB(vhCand.V.mets.at(0).p4.Pt()); // demonstrated to hold also for RunC (used averaged efficiency)
        weightTrig2012SumpT100MET100=triggerWeight.scaleSumpT100MET100_2012B(vhCand.V.mets.at(0).p4.Pt()); // demonstrated to hold also for RunC (used averaged efficiency)
        weightTrig2012APFMET150orDiJetMET=triggerWeight.scalePFMET150orDiJetMET_2012A(vhCand.V.mets.at(0).p4.Pt()); 
        weightTrig2012BPFMET150orDiJetMET=triggerWeight.scalePFMET150orDiJetMET_2012B(vhCand.V.mets.at(0).p4.Pt()); 
        weightTrig2012CPFMET150orDiJetMET=triggerWeight.scalePFMET150orDiJetMET_2012C(vhCand.V.mets.at(0).p4.Pt()); 
	  }

	  if( Vtype == VHbbCandidate::Znn )
	  {
		nvlep=0;
		float weightTrig1 = triggerWeight.scaleMetHLT(vhCand.V.mets.at(0).p4.Pt());
		weightTrigMETLP = weightTrig1;
		weightTrig = weightTrigMET150 + weightTrigMET802CJet  - weightTrigMET802CJet*weightTrigMET150; 
		//weightTrig2012A = 
	  }
      
      if(weightTrigMay < 0) weightTrigMay=weightTrig;
      if(weightTrigV4 < 0) weightTrigV4=weightTrig;
      if(!isMC_)
      {
		weightTrig = 1.; 
		weightTrigMay = 1.;
		weightTrigV4 = 1.;
		weightEleRecoAndId= 1.;
		weightEleTrigJetMETPart= 1.;
		weightEleTrigElePart= 1.;
		weightEleTrigEleAugPart=1.;
		weightTrigMET80= 1.;
		weightTrigMET100= 1.;
		weightTrig2CJet20= 1.;
		weightTrigMET150= 1.;
		weightTrigMET802CJet= 1.;
		weightTrigMET1002CJet= 1.;
		weightTrigMETLP = 1.;
      }

//==fill additional lepton==
//NOTE : move outside of candidate check (Duong, 10-08-2015)
//      aLeptons.reset();
//      nalep=0;
//      if(fromCandidate)
//      {
//		for(size_t j=firstAddMu; j< vhCand.V.muons.size();    j++) aLeptons.set(vhCand.V.muons[j],     nalep++, 13, aux);
//		for(size_t j=firstAddEle;j< vhCand.V.electrons.size();j++) aLeptons.set(vhCand.V.electrons[j], nalep++, 11, aux);
//      }
//      else
//      {
//          for(size_t j=0;j< iEvent->muInfo.size();j++)
//            if((j!= vhCand.V.firstLeptonOrig && j!= vhCand.V.secondLeptonOrig) || ((Vtype != VHbbCandidate::Wmun ) && (Vtype != VHbbCandidate::Zmumu )) )
//			aLeptons.set(iEvent->muInfo[j],nalep++,13,aux);
//          for(size_t j=0;j< iEvent->eleInfo.size();j++)
//            if((j!= vhCand.V.firstLeptonOrig && j!= vhCand.V.secondLeptonOrig) || ((Vtype != VHbbCandidate::Wen ) && (Vtype != VHbbCandidate::Zee )))
//			aLeptons.set(iEvent->eleInfo[j],nalep++,11,aux);
//      }
	  //std::cout << "Leptons done" << std::endl;

      if(isMC_)
      {
	    //std::cout << "BTAGSF " <<  btagJetInfos.size() << " " << btag.weight<BTag1Tight2CustomFilter>(btagJetInfos) << std::endl;
	    if ( btagJetInfos.size()< 10) 
	    {
		  btag1T2CSF = btag.weight<BTag1Tight2CustomFilter>(btagJetInfos);
		  btag2TSF = btag.weight<BTag2TightFilter>(btagJetInfos);
		  btag1TSF = btag.weight<BTag1TightFilter>(btagJetInfos);
		  btagA0CSF = btag.weight<BTagAntiMax0CustomFilter>(btagJetInfos);
		  btagA0TSF = btag.weight<BTagAntiMax0TightFilter>(btagJetInfos);
		  btag2CSF = btag.weight<BTag2CustomFilter>(btagJetInfos);
		  btag1TA1C = btag.weight<BTag1TightAndMax1CustomFilter>(btagJetInfos);
	    }
	    else
	    {
		  std::cout << "WARNING:  combinatorics for " << btagJetInfos.size() << " jets is too high (>=10). use SF=1 " << std::endl;
		  //TODO: revert to random throw  for this cases
		  btag1T2CSF = 1.;
		  btag2TSF =  1.;
		  btag1TSF =  1.;
		  btagA0CSF = 1.;
		  btagA0TSF = 1.;
		  btag2CSF = 1.;
		  btag1TA1C = 1.;
	    }
	  }
            
      if(maxBtag > -99999)
      { 
	    TopHypo topQuark = TopMassReco::topMass(leptonForTop,bJet,vhCand.V.mets.at(0).p4);
	    top.mass = topQuark.p4.M();
	    top.pt = topQuark.p4.Pt();
	    top.wMass = topQuark.p4W.M();
      }
	  else
	  {
		top.mass = -99;
		top.pt = -99;
		top.wMass = -99;
	  }
  
	   //FIXME: too much  warnings... figure out why 
	   //gendrcc=aux.genCCDeltaR(); 
	   //gendrbb=aux.genBBDeltaR(); 

	  genZpt=aux.mcZ.size() > 0 ? aux.mcZ[0].p4.Pt():-99;
	  genWpt=aux.mcW.size() > 0 ? aux.mcW[0].p4.Pt():-99;

	  // Z* is status=3 and Nmother=2 (q and qbar)  
	  // Z is status=2 and Ndau=2 (mup and mum)  
      for (unsigned int i=0; i<aux.mcZ.size(); i++)
	  {
		if (aux.mcZ[i].status==3)
		{
		  genZstar.mass = aux.mcZ[i].p4.M();
		  genZstar.pt = aux.mcZ[i].p4.Pt();
		  //genZstar.eta = aux.mcZ[i].p4.Eta();
		  if(genZstar.pt>0.1) genZstar.eta = aux.mcZ[i].p4.Eta();
		  else genZstar.eta=-99;
		  genZstar.phi = aux.mcZ[i].p4.Phi();
		  genZstar.status = aux.mcZ[i].status;
		  genZstar.charge = aux.mcZ[i].charge;
		  if(aux.mcZ[i].momid!=-99) genZstar.momid  =  aux.mcZ[i].momid;
		}

        if (aux.mcZ[i].dauid.size()>1)
		{
		  if( abs(aux.mcZ[i].dauid[0])==11 || abs(aux.mcZ[i].dauid[0])==12 ||
			  abs(aux.mcZ[i].dauid[0])==13 || abs(aux.mcZ[i].dauid[0])==14 ||
			  abs(aux.mcZ[i].dauid[0])==15 || abs(aux.mcZ[i].dauid[0])==16    )
		  {
			genZ.mass = aux.mcZ[i].p4.M();
			genZ.pt = aux.mcZ[i].p4.Pt();
			//genZ.eta = aux.mcZ[i].p4.Eta();
			if(genZ.pt>0.1) genZ.eta = aux.mcZ[i].p4.Eta();
			else genZ.eta=-99;
			genZ.phi = aux.mcZ[i].p4.Phi();
			genZ.status = aux.mcZ[i].status;
			genZ.charge = aux.mcZ[i].charge;
			if(aux.mcZ[i].momid!=-99) genZ.momid = aux.mcZ[i].momid;
		  }
		}
	  }

      for (unsigned int i=0; i<aux.mcW.size(); i++)
	  {
        if(aux.mcW[i].momid==6 && aux.mcW[i].dauid.size()>1)
		{
		  genTop.wdau1mass = aux.mcW[i].dauFourMomentum[0].M();
		  genTop.wdau1pt   = aux.mcW[i].dauFourMomentum[0].Pt();
		  if(genTop.wdau1pt > 0.1) genTop.wdau1eta= aux.mcW[i].dauFourMomentum[0].Eta();
          genTop.wdau1phi  = aux.mcW[i].dauFourMomentum[0].Phi();
          genTop.wdau1id   = aux.mcW[i].dauid[0];
		  genTop.wdau2mass = aux.mcW[i].dauFourMomentum[1].M();
		  genTop.wdau2pt   = aux.mcW[i].dauFourMomentum[1].Pt();
		  if(genTop.wdau2pt > 0.1) genTop.wdau2eta= aux.mcW[i].dauFourMomentum[1].Eta();
		  genTop.wdau2phi  = aux.mcW[i].dauFourMomentum[1].Phi();
		  genTop.wdau2id   = aux.mcW[i].dauid[1];
		}

        if(aux.mcW[i].momid==-6 && aux.mcW[i].dauid.size()>1)
		{
          genTbar.wdau1mass = aux.mcW[i].dauFourMomentum[0].M();
          genTbar.wdau1pt   = aux.mcW[i].dauFourMomentum[0].Pt();
		  if(genTbar.wdau1pt>0) genTbar.wdau1eta= aux.mcW[i].dauFourMomentum[0].Eta();
		  genTbar.wdau1phi  = aux.mcW[i].dauFourMomentum[0].Phi();
		  genTbar.wdau1id   = aux.mcW[i].dauid[0];
		  genTbar.wdau2mass = aux.mcW[i].dauFourMomentum[1].M();
		  genTbar.wdau2pt   = aux.mcW[i].dauFourMomentum[1].Pt();
		  if(genTbar.wdau2pt>0) genTbar.wdau2eta= aux.mcW[i].dauFourMomentum[1].Eta();
		  genTbar.wdau2phi  = aux.mcW[i].dauFourMomentum[1].Phi();
		  genTbar.wdau2id   = aux.mcW[i].dauid[1];
	    }

	    if (aux.mcW[i].status==3 )
		{
	      genWstar.mass   = aux.mcW[i].p4.M();
	      genWstar.pt     = aux.mcW[i].p4.Pt();
		  if(genWstar.pt>0.1) genWstar.eta = aux.mcW[i].p4.Eta();
		  genWstar.phi    = aux.mcW[i].p4.Phi();
		  genWstar.status = aux.mcW[i].status;
		  genWstar.charge = aux.mcW[i].charge;
		  if(aux.mcW[i].momid!=-99) genWstar.momid = aux.mcW[i].momid;
		}

		if (aux.mcW[i].dauid.size()>1 && (abs(aux.mcW[i].dauid[0])==13 || abs(aux.mcW[i].dauid[0])==11 ))
		{
		  genW.mass   = aux.mcW[i].p4.M();
		  genW.pt     = aux.mcW[i].p4.Pt();
		  if(genW.pt>0) genW.eta = aux.mcW[i].p4.Eta();
		  genW.phi    = aux.mcW[i].p4.Phi();
		  genW.status = aux.mcW[i].status;
		  genW.charge = aux.mcW[i].charge;
		  if(aux.mcW[i].momid!=-99) genW.momid = aux.mcW[i].momid;
		}
	  }
	  
	  // b coming from Higgs
	  for (unsigned int i=0; i<aux.mcB.size(); i++)
	  {
		if(abs(aux.mcB[i].momid)!=5)
		{
		  genB.mass   = aux.mcB[i].p4.M();
		  genB.pt     = aux.mcB[i].p4.Pt();
		  if(genB.pt>0) genB.eta = aux.mcB[i].p4.Eta();
		  genB.phi    = aux.mcB[i].p4.Phi();
		  genB.status = aux.mcB[i].status;
		  genB.charge = aux.mcB[i].charge;
		  if(aux.mcB[i].momid!=-99) genB.momid = aux.mcB[i].momid;
		}

        if(aux.mcB[i].momid==6)
		{
		  genTop.bmass   = aux.mcB[i].p4.M();
		  genTop.bpt     = aux.mcB[i].p4.Pt();
		  if(genTop.bpt >0) genTop.beta = aux.mcB[i].p4.Eta();
		  genTop.bphi    = aux.mcB[i].p4.Phi();
		  genTop.bstatus = aux.mcB[i].status;
	    }
	  }
	
	  for (unsigned int i=0; i<aux.mcBbar.size(); i++)
	  {
		if (abs(aux.mcBbar[i].momid)!=5 )
		{
		  genBbar.mass   = aux.mcBbar[i].p4.M();
		  genBbar.pt     = aux.mcBbar[i].p4.Pt();
		  if(genBbar.pt >0) genBbar.eta = aux.mcBbar[i].p4.Eta();
		  genBbar.phi    = aux.mcBbar[i].p4.Phi();
		  genBbar.status = aux.mcBbar[i].status;
		  if(aux.mcBbar[i].momid!=-99) genBbar.momid = aux.mcBbar[i].momid;
	    }
		
        if  ( aux.mcBbar[i].momid==-6 )
		{
		  genTbar.bmass = aux.mcBbar[i].p4.M();
		  genTbar.bpt = aux.mcBbar[i].p4.Pt();
		  if(genTbar.bpt>0) genTbar.beta = aux.mcBbar[i].p4.Eta();
		  genTbar.bphi = aux.mcBbar[i].p4.Phi();
		  genTbar.bstatus = aux.mcBbar[i].status;
	    }

	  }

	  if (aux.mcH.size()>0)
	  {
		genH.mass   = aux.mcH[0].p4.M();
		genH.pt     = aux.mcH[0].p4.Pt();
		if(genH.pt>0) genH.eta = aux.mcH[0].p4.Eta();
		genH.phi    = aux.mcH[0].p4.Phi();
		genH.status = aux.mcH[0].status;
		genH.charge = aux.mcH[0].charge;
		if (aux.mcH[0].momid!=-99) genH.momid =  aux.mcH[0].momid;
	  }

	  WminusMode=-99;
	  WplusMode=-99;
	  for(unsigned int j=0; j< aux.mcW.size();j++)
      {
		for(unsigned int k=0;k< aux.mcW[j].dauid.size();k++)
	    {
	      int idd=abs(aux.mcW[j].dauid[k]);
	      if(idd==11 || idd==13 || idd==15|| (idd<=5 && idd >=1)) 
		  {
			if(WminusMode==-99 && aux.mcW[j].charge ==-1) WminusMode = idd;
			if(WplusMode ==-99 && aux.mcW[j].charge ==+1) WplusMode  = idd;
		  }
	    }

	    // now check if a semileptonic W is also in a bjets....      
	    //if ( ( (WminusMode==11 || WminusMode==13 || WminusMode==15  ) || (WplusMode==11 || WplusMode==13 || WplusMode==15  ))  && deltaR(vhCand.H.jets[0].p4.Eta(),vhCand.H.jets[0].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi())<0.3 ) hJets.isSemiLeptMCtruth[0]=1;
	    //if ( ( (WminusMode==11 || WminusMode==13 || WminusMode==15  ) || (WplusMode==11 || WplusMode==13 || WplusMode==15  ))  && deltaR(vhCand.H.jets[1].p4.Eta(),vhCand.H.jets[1].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi())<0.3 ) hJets.isSemiLeptMCtruth[1]=1;
	    //for( int j=0; j < naJets && j < MAXJ; j++ ) 
	    //  if ((idd==11 || idd==13 || idd==15) && deltaR(vhCand.additionalJets[j].p4.Eta(),vhCand.additionalJets[j].p4.Phi(), aux.mcW[j].p4.Eta(), aux.mcW[j].p4.Phi()) <0.3)
		//	  aJets.isSemiLept[j]=1;
	  }

      /// Compute pull angle from AK7
      if(vhCand.H.HiggsFlag)
	  {
        if(!fromCandidate)
		{
          std::vector<VHbbEvent::SimpleJet> ak7wrt1(iEvent->simpleJets3);
          std::vector<VHbbEvent::SimpleJet> ak7wrt2(iEvent->simpleJets3);
          if(ak7wrt1.size() > 1)
		  {
            CompareDeltaR deltaRComparatorJ1(vhCand.H.jets[0].p4);
            CompareDeltaR deltaRComparatorJ2(vhCand.H.jets[1].p4);
            std::sort( ak7wrt1.begin(),ak7wrt1.end(),deltaRComparatorJ1 );
            std::sort( ak7wrt2.begin(),ak7wrt2.end(),deltaRComparatorJ2 );
            std::vector<VHbbEvent::SimpleJet> ak7_matched;
            // if the matched are different save them
            if(ak7wrt1[0].p4.DeltaR(ak7wrt2[0].p4) > 0.1)
			{
              ak7_matched.push_back(ak7wrt1[0]);
              ak7_matched.push_back(ak7wrt2[0]);
            }
            // else look at the second best
            else
			{ // ak7wrt1 is best
              if( ak7wrt1[1].p4.DeltaR(vhCand.H.jets[0].p4) < ak7wrt2[1].p4.DeltaR(vhCand.H.jets[1].p4))
              {
                ak7_matched.push_back(ak7wrt1[1]);
                ak7_matched.push_back(ak7wrt2[0]);
              }
              else
              {
                ak7_matched.push_back(ak7wrt1[0]);
                ak7_matched.push_back(ak7wrt2[1]);
              }
            }

            CompareJetPt ptComparator;
            std::sort( ak7_matched.begin(),ak7_matched.end(),ptComparator );
            if(ak7_matched[0].p4.DeltaR(vhCand.H.jets[0].p4) < 0.5 and ak7_matched[1].p4.DeltaR(vhCand.H.jets[1].p4) < 0.5)
            {
              deltaPullAngleAK7  = VHbbCandidateTools::getDeltaTheta(ak7_matched[0],ak7_matched[1]);
              deltaPullAngle2AK7 = VHbbCandidateTools::getDeltaTheta(ak7_matched[1],ak7_matched[0]);
            }
          }
        }

      }//HiggsFlag
      } //end if has a Z candidate

//==NOTE : iEvent is modified with jet correction, resolution seamering applied== 
      nalep=0;
      aLeptons.reset();
      int doBelongToVhCand(0) ;
      unsigned firstLeptInd(1000) ;
      unsigned secondLeptInd(1000) ;
      unsigned nGoodLept = nallElectrons + nallMuons ;
      if (cand->size() > 0) { 
        const VHbbCandidate & vhCand =  cand->at(0);
        firstLeptInd = vhCand.V.firstLeptonOrig ;
        secondLeptInd = vhCand.V.secondLeptonOrig ;
      }

      for(size_t j=0;j< iEvent->muInfo.size();j++) {
        if (j == firstLeptInd) doBelongToVhCand = 1 ;
        if (j == secondLeptInd) doBelongToVhCand = 2 ;
        aLeptons.set(iEvent->muInfo[j],nalep++,13,aux, doBelongToVhCand);
        if (iEvent->muInfo[j].p4.Pt() > 15 && fabs(iEvent->muInfo[j].p4.Eta()) < 2.6) nGoodLept += 1 ;  
      }
      for(size_t j=0;j< iEvent->eleInfo.size();j++) { 
        if (j == firstLeptInd) doBelongToVhCand = 1 ;
        if (j == secondLeptInd) doBelongToVhCand = 2 ;
        aLeptons.set(iEvent->eleInfo[j],nalep++,11,aux, doBelongToVhCand);
        if (iEvent->eleInfo[j].p4.Pt() > 15 && fabs(iEvent->eleInfo[j].p4.Eta()) < 2.6) nGoodLept += 1 ;  
      }
      
      aJets.reset();
      naJets = 0 ;
      for(unsigned j=0 ; j < iEvent->simpleJets2.size() ; ++j) {
        if (iEvent->simpleJets2[j].p4.Pt() > 20 && fabs(iEvent->simpleJets2[j].p4.Eta()) < 3) {
          aJets.set(iEvent->simpleJets2[j],naJets, jetWithJEC_p4_vec[j], jetBestMC_p4_vec[j], iEvent->simpleJets2[j].bestMCid) ;
          naJets += 1 ;
        }
      } 
      
      aMET.reset() ;
      apfMET.reset() ; 
      aMET.et = iEvent->pfmetType1corr.p4.Pt();
      aMET.phi = iEvent->pfmetType1corr.p4.Phi();
      aMET.sumet = iEvent->pfmetType1corr.sumEt;
      aMET.sig = iEvent->pfmetType1corr.metSig;
      apfMET.et    = iEvent->pfmet.p4.Pt();
      apfMET.phi   = iEvent->pfmet.p4.Phi();
      apfMET.sumet = iEvent->pfmet.sumEt;
      apfMET.sig   = iEvent->pfmet.metSig;
      
      if(isMC_ && doFillMoreMCtruth) {  
//==fill gen particle from Z+ll==
//Z(3)->ll(3)->ll(2)->X(1): with photon radiation, Z(3)->ll(3)->X(1) no photon radiation
      fwlite::Handle<std::vector<reco::GenParticle> > genPartH ;
      //genPartH.getByLabel(ev,"savedGenParticles");
      genPartH.getByLabel(ev,"genParticles");
      std::vector<reco::GenParticle> genPartV = *(genPartH.product());
//==find Z particle==
      zdecayMode = 0 ;
      int pdgId(0) ;
      int status(-10) ;
      vector<reco::GenParticleRefVector::const_iterator> lepRefItrs ; //status 3 lepton from Z decays
      vector<int> zPartInds ;
      for (unsigned i = 0; i < genPartV.size(); i++) {
        pdgId = genPartV[i].pdgId() ;
        status = genPartV[i].status() ;
        if (pdgId == 23 && status == 3) {
          const reco::GenParticleRefVector& daughterRefs = genPartV[i].daughterRefVector();
          zPartInds.push_back(i) ;
          //cout << "\n ===================" ;
          for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
            //cout<<"    - Daughter "<<(*idr).key()<<" "<<(*idr)->pdgId()<<endl;
            //if (abs((*idr)->pdgId()) == 11) { zdecayMode = 1 ; break ;} 
            //if (abs((*idr)->pdgId()) == 13) { zdecayMode = 2 ; break ;} 
            //if (abs((*idr)->pdgId()) == 15) { zdecayMode = 3 ; break ;} 
            int pdgIdAbs = abs((*idr)->pdgId()) ;
            if ( pdgIdAbs == 11 || pdgIdAbs == 13 || pdgIdAbs == 15 ) lepRefItrs.push_back(idr) ; 
            if (abs((*idr)->pdgId()) == 11) { zdecayMode = 1 ; } 
            if (abs((*idr)->pdgId()) == 13) { zdecayMode = 2 ; } 
            if (abs((*idr)->pdgId()) == 15) { zdecayMode = 3 ; } 
          }

          break ;

        }
      
      } //end loop over genParticle
      
      //===get daughters of lepton meaning lepton with status 2 particle===
      vector<reco::GenParticleRefVector::const_iterator> lepDauItrs ;
      if (lepRefItrs.size() != 2) cout << "\n Warning: not a Z->ll" ;
      if (lepRefItrs.size() == 2) {
        //lepDaughter with status == 2
         //cout << "\n Two leptons from Z are: " << (*lepRefItrs[0])->pdgId() << "  " << (*lepRefItrs[0])->status() << "  " << (*lepRefItrs[1])->pdgId() << "  " << (*lepRefItrs[1])->status() ;
         if (((*lepRefItrs[0])->daughterRefVector()).size() == 1) lepDauItrs.push_back(((*lepRefItrs[0])->daughterRefVector()).begin()) ;
         else cout << "\n Warning: there are " << ((*lepRefItrs[0])->daughterRefVector()).size() << " daughters for first lepton  " << (*lepRefItrs[0])->pdgId() << "  " << (*lepRefItrs[0])->status() ;  
         if (((*lepRefItrs[1])->daughterRefVector()).size() != 0) lepDauItrs.push_back(((*lepRefItrs[1])->daughterRefVector()).begin()) ;
         else cout << "\n Warning: there are " << ((*lepRefItrs[1])->daughterRefVector()).size() << " daughters for second lepton  " << (*lepRefItrs[1])->pdgId() << "  " << (*lepRefItrs[1])->status() ;  
      }
      
      //===fill leptons and photons with status == 1
      nGenLep = 0 ;
      genLeps.reset() ; 
      nGenPho = 0 ;
      genPhos.reset() ;
      if (lepDauItrs.size() == 2) {
        reco::GenParticleRefVector::const_iterator lepOut ;
        vector<reco::GenParticleRefVector::const_iterator> phoItrs ;
        if (FindStatus1Lepton(lepDauItrs[0], lepOut, phoItrs)) genLeps.set(lepOut,nGenLep++) ;
        if (FindStatus1Lepton(lepDauItrs[1], lepOut, phoItrs)) genLeps.set(lepOut,nGenLep++) ;
        for (unsigned iPho = 0 ; iPho < phoItrs.size() ; ++iPho) genPhos.set(phoItrs[iPho],nGenPho++) ;
      }

/*      
      if (lepDauItrs.size() == 2) {
        int status = (*lepDauItrs[0])->status() ;
        int pdgId = (*lepDauItrs[0])->pdgId() ;
        //===no photon radiation, status 3 lep -> status 1 lep==
        if (status == 1) {
          genLeps.set(lepDauItrs[0],nGenLep++) ;
        }
        //===with photon radiation or tau decay, status 3 lep -> status 2 lep====
        else if (status == 2) {
          const reco::GenParticleRefVector& particleRefs = (*lepDauItrs[0])->daughterRefVector() ;
          for(reco::GenParticleRefVector::const_iterator idr = particleRefs.begin(); idr!= particleRefs.end(); ++idr) {
            int pdgIdAbsTmp = abs((*idr)->pdgId()) ;
            int statusTmp = (*idr)->status() ;
            if ((pdgIdAbsTmp == 11 || pdgIdAbsTmp == 13) && statusTmp == 1) { //only store electron and muons
              genLeps.set(idr,nGenLep++) ;
            } 
            else if ((pdgIdAbsTmp == 22) && statusTmp == 1) {
              genPhos.set(idr,nGenPho++) ;
            }
            else {
              cout << "\n no gen ele, muon, or photon found 0: " << pdgId << "  " << status << "  " << (*idr)->pdgId() << "  " << (*idr)->status() ; 
            } 
          }
        } //end else if (status == 2
        else {
          cout << "\n Warning: Unidentified daughter of lepton 0 status" ;
        }

        status = (*lepDauItrs[1])->status() ;
        pdgId = (*lepDauItrs[1])->pdgId() ;
        //===no photon radiation, status 3 lep -> status 1 lep==
        if (status == 1) {
          genLeps.set(lepDauItrs[1],nGenLep++) ;
        }
        //===with photon radiation or tau decay, status 3 lep -> status 2 lep====
        else if (status == 2) {
          const reco::GenParticleRefVector& particleRefs = (*lepDauItrs[1])->daughterRefVector() ;
          for(reco::GenParticleRefVector::const_iterator idr = particleRefs.begin(); idr!= particleRefs.end(); ++idr) {
            int pdgIdAbsTmp = abs((*idr)->pdgId()) ;
            int statusTmp = (*idr)->status() ;
            if ((pdgIdAbsTmp == 11 || pdgIdAbsTmp == 13) && statusTmp == 1) { //only store electron and muons
              genLeps.set(idr,nGenLep++) ;
            } 
            else if ((pdgIdAbsTmp == 22) && statusTmp == 1) {
              genPhos.set(idr,nGenPho++) ;
            }
            else {
              cout << "\n no gen ele, muon, or photon found 1: " << pdgId << "  " << status << "  " << (*idr)->pdgId() << "  " << (*idr)->status() ; 
            } 
          }
        } //end else if (status == 2
        else {
          cout << "\n Warning: Unidentified daughter of lepton 1 status" ;
        }
 
      } //end lepDauItrs.size() == 2      
*/       
      if (lepDauItrs.size() != 2) cout << "\n Warning: found " << lepDauItrs.size() << " daughter lepton (should be always 2)" ;

//===fill status 3 leptons + partons===
      nGenSta3obj = 0 ;
      genSta3objs.reset() ;
      if (zPartInds.size() != 1) cout << "\n Warning: there are " << zPartInds.size() << " status 3 Z found, weird, will not fill" ;
      else {
      const reco::GenParticleRefVector& zMomRefs = genPartV[zPartInds[0]].motherRefVector() ; //mother of Z, neeeded to find Z sibling
      if (zMomRefs.size() != 2) cout << "\n Warning: there are " << zMomRefs.size() << " mothers of Z, will not fill sibling" ;
      else {
        //===fill two leptons status 3===
        if (lepRefItrs.size() != 2) cout << "\n Warning: odd Z decay, can't find 2 leptons, will not fill sibling" ;
        else {
          genSta3objs.set(lepRefItrs[0], nGenSta3obj++) ;
          genSta3objs.set(lepRefItrs[1], nGenSta3obj++) ;
        } 
        //===fill parton siblin with Z===
        const reco::GenParticleRefVector& daughterRefs = zMomRefs[0]->daughterRefVector() ;
        const reco::GenParticleRefVector& daughterRef1s = zMomRefs[1]->daughterRefVector() ;
        bool doPassSameDaughter(true) ;
        for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
          //===check that both moms give the same daughters==
          int pdgId = (*idr)->pdgId() ;
          int status = (*idr)->status() ;
          bool foundIt = false ;
          for(reco::GenParticleRefVector::const_iterator idr1 = daughterRef1s.begin(); idr1!= daughterRef1s.end(); ++idr1) {
            
            int pdgId1 = (*idr1)->pdgId() ;
            int status1 = (*idr1)->status() ;
            if (pdgId == pdgId1 && status == status1) {
              foundIt = true ;
              break ;
            }
          }
          if (!foundIt) {
            cout << "\n Warning: two mothers (0 compares to 1) not giving the same daughters, will not fill" ; 
            doPassSameDaughter = false ;
            break ;
          }
        }
        
        for(reco::GenParticleRefVector::const_iterator idr = daughterRef1s.begin(); idr!= daughterRef1s.end(); ++idr) {
          //===check that both moms give the same daughters==
          int pdgId = (*idr)->pdgId() ;
          int status = (*idr)->status() ;
          bool foundIt = false ;
          for(reco::GenParticleRefVector::const_iterator idr1 = daughterRefs.begin(); idr1!= daughterRefs.end(); ++idr1) {
            
            int pdgId1 = (*idr1)->pdgId() ;
            int status1 = (*idr1)->status() ;
            if (pdgId == pdgId1 && status == status1) {
              foundIt = true ;
              break ;
            }
          }
          if (!foundIt) {
            cout << "\n Warning: two mothers (1 compares to 0) not giving the same daughters, will not fill" ; 
            doPassSameDaughter = false ;
            break ;
          }
        }

        if (doPassSameDaughter) {
          for(reco::GenParticleRefVector::const_iterator idr = daughterRefs.begin(); idr!= daughterRefs.end(); ++idr) {
            if ((*idr)->pdgId() != 23) genSta3objs.set(idr, nGenSta3obj++) ;
          }

        }
   
        } //end if zMom size == 2
           
 
      } //end elseif of zPartInds.size() != 2
      
      //===Now fill gen jets====
      //fwlite::Handle<std::vector<reco::JetFlavourMatching> > genJetH ;
      fwlite::Handle<reco::JetFlavourMatchingCollection> genJetH ;
      genJetH.getByLabel(ev,"flavourByVal"); //for ak5GenJet
      nGenAk5Jet = 0 ;
      genAk5Jets.reset() ;
      for ( reco::JetFlavourMatchingCollection::const_iterator j  = genJetH->begin(); j != genJetH->end() && nGenAk5Jet < MAXGENOBJ ; j ++ ) {
        genAk5Jets.set(j, nGenAk5Jet++);
      }
      
      genJetH.getByLabel(ev,"flavourByVal1"); //for patJetsPFlow
      nGenPatPFjet = 0 ;
      genPatPFjets.reset() ;
      for ( reco::JetFlavourMatchingCollection::const_iterator j  = genJetH->begin(); j != genJetH->end() && nGenPatPFjet < MAXGENOBJ ; j ++ ) {
        genPatPFjets.set(j, nGenPatPFjet++);
      }

       
      } //end if (isMC_)
      
      //==NOTE : for data require at least one good lepton = nallElectrons + nallMuons + naLepton (pt > 15 and fabs(eta) < 2.6)== Duong (10-09-2015)
      if (!isMC_ && (nGoodLept == 0)) continue ;  

      hEvent->Fill(2.5) ; //number of events filled

      _outTree->Fill();

    
    }// closed event loop

    std::cout << "closing the file: " << inputFiles_[iFile] << std::endl;
    inFile->Close();  // close input file

  } // loop on files

  
  std::cout << "Events: " << ievt <<std::endl;
  std::cout << "TotalCount: " << totalcount <<std::endl;

  _outFile->cd();
  //_outTree->Write();
  _outFile->Write();
  _outFile->Close();

  return 0;
}


