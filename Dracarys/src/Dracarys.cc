// -*- C++ -*-
//
// Package:    VLF_ANA/Dracarys
// Class:      Dracarys
// 
/**\class Dracarys Dracarys.cc VLF_ANA/Dracarys/plugins/Dracarys.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
//Original Authors:  Jose Ruiz
//                   Camilo Salazar
//         Created:  Tue, 08 Aug 2017 09:38:40 GMT
//



// system include files
#include <memory>
#include <cmath>
#include <string> /*Using of strigs*/

// user include files

#include "VLF_ANA/Dracarys/interface/Dracarys.h"

// Counters
int Dracarys::NoCuts;
int Dracarys::TriggerPathCut;
int Dracarys::GoodVertex;
int Dracarys::aJetatLessCut;
int Dracarys::LeadingMuPtM3;
int Dracarys::MissingETCut;
int Dracarys::BasicJetsCut;
int Dracarys::bJetsCut;
int Dracarys::MuonMetMTCut;

Dracarys::Dracarys(const edm::ParameterSet& iConfig):
  triggerBits_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("bits"))),
  triggerObjects_(consumes<pat::TriggerObjectStandAlone>(iConfig.getParameter<edm::InputTag>("objects"))),
  triggerPrescales_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("prescales"))),
  tok_vertex_(consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("vertices"))),
  tok_pileup_(consumes<std::vector<PileupSummaryInfo> >(iConfig.getParameter<edm::InputTag>("pileupInfo"))),
  tok_jets_(consumes<edm::View<pat::Jet> >(iConfig.getParameter<edm::InputTag>("objet"))),
  tok_met_(consumes<pat::METCollection>(iConfig.getParameter<edm::InputTag>("obmet"))),
  tok_muons_(consumes<edm::View<pat::Muon> >(iConfig.getParameter<edm::InputTag>("obmuon")))
  //,histContainer_()
{
  //now do whatever initialization is needed
  usesResource("TFileService"); 
  // register to the TFileService
  edm::Service<TFileService> fs;
  //Is Data Boolean
  is_data_ = (iConfig.getParameter<bool>("is_data"));
  //trigger variables
  isTrigger_ = (iConfig.getParameter<bool>("isTrigger"));
  isTriggerToo_ = (iConfig.getParameter<bool>("isTriggerToo"));
  TriggerPath1_ = (iConfig.getParameter<string>("TriggerPath1"));
  TriggerPath2_ = (iConfig.getParameter<string>("TriggerPath2"));
  //Debugging option
  debug_ = (iConfig.getParameter<bool>("debug"));
  //Cuts
  Pvtx_ndof_min_ = (iConfig.getParameter<int>("Pvtx_ndof_min"));
  Pvtx_vtx_max_ = (iConfig.getParameter<double>("Pvtx_vtx_max"));
  Pvtx_vtxdxy_max_ = (iConfig.getParameter<double>("Pvtx_vtxdxy_max"));
  MinMuonPt_ = (iConfig.getParameter<double>("MinMuonPt"));
  MaxMuonPt_ = (iConfig.getParameter<double>("MaxMuonPt"));
  MuonIso_ = (iConfig.getParameter<double>("MuonIso"));
  MuonIsoCut_ = (iConfig.getParameter<bool>("MuonIsoCut"));
  MuonID_ = (iConfig.getParameter<int>("MuonID"));
  MuonIDCut_ = (iConfig.getParameter<bool>("MuonIDCut"));
  MinNMuons_ = (iConfig.getParameter<int>("MinNMuons"));
  MaxNMuons_ = (iConfig.getParameter<int>("MaxNMuons"));
  MinMET_ = (iConfig.getParameter<double>("MinMET"));
  MinJetPt_ = (iConfig.getParameter<double>("MinJetPt"));
  MaxJetEta_ = (iConfig.getParameter<double>("MaxJetEta"));
  MinNJets_ = (iConfig.getParameter<int>("MinNJets"));
  MaxNJets_ = (iConfig.getParameter<int>("MaxNJets"));
  bJetTag_ = (iConfig.getParameter<double>("bJetTag"));
  MinbJetPt_ = (iConfig.getParameter<double>("MinbJetPt"));
  MaxbJetEta_ = (iConfig.getParameter<double>("MaxbJetEta"));
  MinNbJets_ = (iConfig.getParameter<int>("MinNbJets"));
  MaxNbJets_ = (iConfig.getParameter<int>("MaxNbJets"));
  MinMTMuonMet_ = (iConfig.getParameter<double>("MinMTMuonMet"));
  MaxMTMuonMet_ = (iConfig.getParameter<double>("MaxMTMuonMet"));
  //Create a TTree
  tree_ = fs->make<TTree>("VLFTree","VLFTree");
}


Dracarys::~Dracarys()
{
  //delete tree_;
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
  
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Dracarys::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace edm;
  
  
  //Trigger
  edm::Handle<edm::TriggerResults> triggerBits;
  edm::Handle<pat::TriggerObjectStandAlone> triggerObjects;
  edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
  
  iEvent.getByToken(triggerBits_, triggerBits);
  iEvent.getByToken(triggerObjects_, triggerObjects);
  iEvent.getByToken(triggerPrescales_, triggerPrescales);
  
  //Vertex info
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_vertex_, vertices);
  
  /*Handling Muons*/
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(tok_muons_, muons);
  /*Handling Jets*/
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(tok_jets_,jets);
  /*Handling MET*/   
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(tok_met_,mets);
  const pat::MET &met = mets->front();
  
  ////////////////////////////

  //Counting events
  Dracarys::NoCuts++;
  
  //Cleaning all variables
  Clean();
  
  bool flagIsTrigger=false;
  bool flagIsTriggerToo=false;
  
  const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
  //   std::cout << "\n == TRIGGER PATHS= " << std::endl;
  
  if( TriggerPath2_ == "" || !isTriggerToo_ ){
    flagIsTriggerToo = true;
  } else{
    flagIsTriggerToo = false;
  } 
  
  if( TriggerPath1_ == "" || !isTrigger_ ){
    flagIsTrigger = true;
  } else{
    flagIsTrigger = false;
  }
  
  
  if ( debug_ ){
    //here we can see if the trigger we ask for
    if (isTrigger_) std::cout << std::endl << std::endl << " Looking for the Trigger : " << TriggerPath1_ << std::endl;
    if (isTriggerToo_) std::cout << " And the Trigger         : " << TriggerPath2_ << std::endl; 
    if (!isTriggerToo_ && !isTrigger_)  std::cout << std::endl << "No Trigger " << std::endl; 
  }
  bool trigNameflag = true;
  for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
    if (flagIsTrigger == true && flagIsTriggerToo == true) {break;}
    
    
    /*Cut the version of the trigger path*/
    std::string nameV = names.triggerName(i);
    std::string version ="_v";
    size_t start_pos = nameV.rfind(version);
    std::string TriggerNameVersionOff;
    if (start_pos != std::string::npos){
      TriggerNameVersionOff= nameV.erase(start_pos, version.length()+4);
    }
    if ( debug_ && trigNameflag ) {
      std::cout << std::endl << "Trigger Paths: ";
      trigNameflag = false;
    }
    if ( debug_ && triggerBits->accept(i) ) std::cout << nameV << ", ";
    

    if( !flagIsTrigger && isTrigger_ ){
      //if (debug_) std::cout << "Evaluating TriggerPath 1 .................." << std::endl;
      if(TriggerNameVersionOff ==  TriggerPath1_){
	if ( debug_ && triggerBits->accept(i) ){
	  ///here we can see if the trigger was used
	  std::cout << std::endl << std::endl << "Evaluating TriggerPath 1 .................." << std::endl;
	  //std::cout << "TriggerName " << i << "          : "<< nameV << std::endl;
	  std::cout << "TriggerWithOutVersion " << i << ": "<< TriggerNameVersionOff << std::endl;
	  trigNameflag = true;
	}
	
	if(triggerBits->accept(i)){
	  flagIsTrigger = true;
	  if ( debug_ && flagIsTrigger )  std::cout << "!!!! TriggerPath 1 PASS!!!" << std::endl << std::endl;
	}
      }
    }
    if( !flagIsTriggerToo && isTriggerToo_ ){
      //if ( debug_ )  std::cout << "Evaluating TriggerPath 2.................." << std::endl;
      if(TriggerNameVersionOff ==  TriggerPath2_){
	if ( debug_ && triggerBits->accept(i) ){
	  ///here we can see if the trigger was used
	  std::cout << std::endl << std::endl << "Evaluating TriggerPath 2.................." << std::endl;
	  //std::cout << "TriggerName " << i << "          : "<< nameV << std::endl;
	  std::cout << "TriggerWithOutVersion " << i << ": "<< TriggerNameVersionOff << std::endl;
	  trigNameflag = true;
	}
	if(triggerBits->accept(i)){
	  flagIsTriggerToo = true;
	  if ( debug_ && flagIsTriggerToo )  std::cout << "!!!TriggerPath 2 PASS!!!" << std::endl << std::endl;
	}
      }
    }
  }

     
  if( flagIsTrigger && flagIsTriggerToo ){
    
    //Counting number of events which pass the triggers
    Dracarys::TriggerPathCut++;
    
    //Requiring a good primery vertex
    bool flagGoodVertex = false;
    
    Nvertices=0;
    reco::VertexCollection::const_iterator firstGoodVertex = vertices->end();
    for(reco::VertexCollection::const_iterator vtxIt = vertices->begin(); vtxIt!= vertices->end(); ++vtxIt) {
      if((vtxIt->isValid()) && !(vtxIt->isFake())) {
	if(vtxIt->ndof() < Pvtx_ndof_min_) continue; 
	if(abs(vtxIt->z()) >= Pvtx_vtx_max_) continue;
	if(sqrt((vtxIt->x()*vtxIt->x()) + (vtxIt->y()*vtxIt->y())) >= Pvtx_vtxdxy_max_) continue; 
	flagGoodVertex=true;
	if (Nvertices==0) firstGoodVertex = vtxIt;
	Nvertices++;
      }
    }
    //For the tighmuon
    
    reco::Vertex vertex = vertices->at(0);


    if (!flagGoodVertex) return;
    Dracarys::GoodVertex++;
    
    if(!is_data_) {
      Handle<std::vector< PileupSummaryInfo > > PUInfo;
      iEvent.getByToken(tok_pileup_, PUInfo); 
      
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      NObservedInTimePUVertices = 0;
      NTruePUInteractions = 0;
      for(PVI = PUInfo->begin(); PVI != PUInfo->end(); ++PVI) { //loop over pileup info
	if(PVI->getBunchCrossing() == 0) { NObservedInTimePUVertices += PVI->getPU_NumInteractions(); } // number of observed in-time pileup interactions
	if(PVI->getBunchCrossing() == 0) { NTruePUInteractions = PVI->getTrueNumInteractions(); } // number of true in-time pileup interactions
      }
    } 

    ////////////////////////////////////////
    
    if( (int)jets->size() >= MinNJets_ && jets->size() > 0 ){ //Min number apply CAS
      Dracarys::aJetatLessCut++;

      if (debug_) std::cout << "MuonSize = " << muons->size() << std::endl;


      if( (int) muons->size() >= MinNMuons_ &&  muons->size() > 0 ){ //Min number apply CAS
	
	bool flagMuonChooser=false;
	int OurMuonDefinitionCounter=0;
	
	//Temporary vector muon info container
	std::vector<XYZTLorentzVector> tempMuons;
	std::vector<double> tempMuon_charge, tempCombined_Iso;
	std::vector <bool> tempMuon_loose, tempMuon_medium, tempMuon_tight;
	int MuTigIDCount = 0;
	int MuMedIDCount = 0;
	int MuLooIDCount = 0;
	int MuISOCount = 0;
	int MuIDCount =0;

	if (debug_) std::cout << "Muon counter loop" << std::endl;
	for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
	  
	  bool flagMuIso = false;
	  bool flagMuID = false;
	  bool flagPtCut = false;
	  double MuonIso = (muon->pfIsolationR04().sumChargedHadronPt	\
			    +  max(0., muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - \
				   0.5*muon->pfIsolationR04().sumPUPt))/muon->pt();
	  ////MUON ISO////
	  if (MuonIso < MuonIso_) {
	    flagMuIso = true;
	    MuISOCount++;
	  }
	  ////MUON ID////
	  if ( muon->isLooseMuon() ){
	    MuLooIDCount++;
	    if ( MuonID_==0 ) {
	      flagMuID = true;
	      MuIDCount++;
	    }
	  }
	  if ( muon->isMediumMuon() ){
	    MuMedIDCount++;
	    if ( MuonID_==1 ) {
	      flagMuID = true;
	      MuIDCount++;
	    }
	  }
	  if ( muon->isTightMuon(vertex) ){
	    MuTigIDCount++;
	    if ( MuonID_==2 ) {
	      flagMuID = true;
	      MuIDCount++;
	    }
	  }
	  ////MUON Kinematic////

	  if( (muon->pt() > MinMuonPt_) && (muon->pt() < MaxMuonPt_) ) {
	    flagPtCut = true;
	  }else{
	    if (debug_) std::cout <<"Kinematics Fail"<< std::endl;
	  }
	  
	  
	  if(!MuonIDCut_) {//Is the Muon ID required?
	    flagMuID = true;
	   if (debug_) std::cout <<"No MuID"<< std::endl;
	  }
	  if(!MuonIsoCut_) {//Is the Muon ISO required?
	    flagMuIso = true; 
	    if (debug_) std::cout <<"No MuISO"<< std::endl;
	  }

	  if (debug_) {
	    std::cout << "Muon pt = " << muon->pt() << ", Muon Iso=" << MuonIso << ", Medium ID=" << muon->isMediumMuon();
	    std::cout << "TightID = " << muon->isTightMuon(vertex) << ", Loose ID=" << muon->isLooseMuon();
	    std::cout << ", flagPtCut =" << flagPtCut ;
	    std::cout <<", flagMuID = " << flagMuID << ", flagMuIso = "<< flagMuIso << std::endl;
	  }
	  
	  
	  
	  if ( flagMuID && flagMuIso && flagPtCut ){
	    
	    OurMuonDefinitionCounter++;
	    XYZTLorentzVector mu(muon->px(), muon->py(), muon->pz(), muon->energy());
	    tempMuons.push_back(mu);
	    tempMuon_charge.push_back(muon->charge());
	    tempCombined_Iso.push_back(MuonIso);
	    tempMuon_loose.push_back(muon->isLooseMuon()); 
	    tempMuon_medium.push_back(muon->isMediumMuon()); 
	    tempMuon_tight.push_back(muon->isTightMuon(*firstGoodVertex));
	  }
	  if (debug_) std::cout <<"OurMuonDefinitionCounter"<<OurMuonDefinitionCounter << std::endl;
	  
	}//endfor muons
	
	if ( (OurMuonDefinitionCounter>=MinNMuons_) && (OurMuonDefinitionCounter<=MaxNMuons_)) {
	  flagMuonChooser=true;
	  if (debug_) std::cout <<"MuonChooser PASS" << std::endl;
	}
	
	if ( flagMuonChooser ){
	  Dracarys::LeadingMuPtM3++;
	  //TTree Filling
	  //NMuons=OurMuonDefinitionCounter;
	  NMuons= muons->size();//The multiplicity without cuts
	  AnaMuons = tempMuons;
	  Muon_charge = tempMuon_charge;
	  Combined_Iso = tempCombined_Iso;
	  Muon_loose = 	tempMuon_loose; 
	  Muon_medium = tempMuon_medium; 
	  Muon_tight = tempMuon_tight;
	  NMuonstight = MuTigIDCount;
	  NMuonsmedium = MuMedIDCount;
	  NMuonsloose = MuLooIDCount;
	  NMuonsIso = MuISOCount;
	  NMuonsID =MuIDCount;
	  
	  //MET selection
	  if (met.pt() < MinMET_) return;
	  Dracarys::MissingETCut++;
	  MET = XYZTLorentzVector(met.px(), met.py(), met.pz(), met.energy());
	  
	  //Basic jets selection
	  bool flagAtLeastOneGoodJet=false;
	  NJets=0;
	  NbJets=0;
	  for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
	    XYZTLorentzVector je(jet->px(), jet->py(), jet->pz(), jet->energy());
	    if (debug_) std::cout <<"Jet Pt: " << je.Pt() <<std::endl;
	    if( (jet->pt() > MinJetPt_) && (abs(jet->eta()) < MaxJetEta_) )
	      {
		AnaJets.push_back(je);
		bJetDiscriminator.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
		NJets++;
	      }
	    if( (jet->pt() > MinbJetPt_) && (abs(jet->eta()) < MaxbJetEta_) && (jet->bDiscriminator("combinedSecondaryVertexBJetTags")>bJetTag_) )
	      {
		NbJets++;
	      }
	  }
	  
	  if ( (NJets >= MinNJets_) && (NJets <= MaxNJets_) ) flagAtLeastOneGoodJet=true;
	  if (!flagAtLeastOneGoodJet) return;
	  Dracarys::BasicJetsCut++;
	  
	  //bJets cut
	  bool flagbJetCut=false;
	  if ( (NbJets >= MinNbJets_) && (NbJets <= MaxNbJets_) ) flagbJetCut=true;
	  if (!flagbJetCut) return;
	  Dracarys::bJetsCut++;
	  
	  //Analysis only part
	  TLorentzVector LeadingMuon, Met;
	  LeadingMuon.SetPxPyPzE(AnaMuons[0].px(), AnaMuons[0].py(), AnaMuons[0].pz(), AnaMuons[0].energy());
	  Met.SetPxPyPzE(MET.px(), MET.py(), MET.pz(), MET.energy());
	  MT_LeadingMuon_MET = sqrt(2*LeadingMuon.Pt()*Met.Pt()*(1-cos(LeadingMuon.DeltaPhi(Met))));
	  
	  if ( (MT_LeadingMuon_MET < MinMTMuonMet_) || (MT_LeadingMuon_MET > MaxMTMuonMet_) ) return;
	  MuonMetMTCut++;
	  
	  tree_->Fill();
	  
	}else{//endif flagMuonChooser
	  if (debug_) std::cout <<"OurMuonDefinition FAIL" <<std::endl;
	}
      }else{//endif MuonSize
	if (debug_) std::cout <<"MuonSize FAIL" <<std::endl;
      }
    }else{//endif JetSize 
      if (debug_) std::cout <<"Jet Number FAIL" <<std::endl;
    }
  }//endif Triggers    

}   



void 
Dracarys::Clean()
{
  AnaMuons.clear();
  AnaJets.clear();
  MET = XYZTLorentzVector(0.0, 0.0, 0.0, 0.0);
  Nvertices=0;
  NObservedInTimePUVertices=0;
  NTruePUInteractions=0;
  Muon_charge.clear();
  Combined_Iso.clear();
  Muon_loose.clear();
  Muon_medium.clear();
  Muon_tight.clear();
  bJetDiscriminator.clear();
  NMuons=0;
  NJets=0;
  NbJets=0;
  MT_LeadingMuon_MET=0;
  NMuonstight=0;
  NMuonsmedium=0;
  NMuonsloose=0;
  NMuonsIso=0;
  NMuonsID=0;
}

// ------------ method called once each job just before starting event loop  ------------
void Dracarys::beginJob()
{

  Dracarys::NoCuts=0; 
  Dracarys::TriggerPathCut=0;
  Dracarys::GoodVertex=0;
  Dracarys::aJetatLessCut=0;
  Dracarys::LeadingMuPtM3=0;
  Dracarys::MissingETCut=0;
  Dracarys::BasicJetsCut=0;
  Dracarys::bJetsCut=0;
  Dracarys::MuonMetMTCut=0;

  //Tree Structure -> branches should be declared in decreasing size                                                                                               
  tree_->Branch("AnaMuons",&AnaMuons);
  tree_->Branch("AnaJets",&AnaJets);
  tree_->Branch("AnaMET",&MET);
  tree_->Branch("combinedSecondaryVertexbJetDiscriminator",&bJetDiscriminator);
  tree_->Branch("Combined_iso_DeltaBetaPU",&Combined_Iso);
  tree_->Branch("Muon_charge",&Muon_charge);
  tree_->Branch("MuonLooseID",&Muon_loose);
  tree_->Branch("MuonMediumID",&Muon_medium);
  tree_->Branch("MuonTightID",&Muon_tight);
  tree_->Branch("MuonMultiplicity",&NMuons);
  tree_->Branch("JetsMultiplicity",&NJets);
  tree_->Branch("bJetsMultiplicity",&NbJets);
  tree_->Branch("MT_LeadingMuon_MET",&MT_LeadingMuon_MET);
  tree_->Branch("Vertices",&Nvertices);
  tree_->Branch("InTimePU",&NObservedInTimePUVertices);
  tree_->Branch("TruePU",&NTruePUInteractions);
  tree_->Branch("MuTightIDMultiplicity",&NMuonstight);
  tree_->Branch("MuMediumIDMultiplicity",&NMuonsmedium);
  tree_->Branch("MuLooseIDMultiplicity",&NMuonsloose);
  tree_->Branch("MuIsoMultiplicity",&NMuonsIso);
  tree_->Branch("MuIDPassMultiplicity",&NMuonsID);
  
  // book histograms:
  // histContainer_["muons"  ]=fs->make<TH1F>("muons",   "muon multiplicity",     10, 0,  10);
  // histContainer_["jets"   ]=fs->make<TH1F>("jets",    "jet multiplicity",      10, 0,  10);
  // histContainer_["met"    ]=fs->make<TH1F>("met",     "missing E_{T}",         20, 0, 500);

}

// ------------ method called once each job just after ending the event loop  ------------
void Dracarys::endJob() 
{
  std::cout<< std::endl << std::endl;
  std::cout<< "++++++ CUT FLOW ++++++ " <<endl;
  std::cout<< "NoCuts         = "<< NoCuts <<endl;
  std::cout<< "TriggerPathCut = "<< TriggerPathCut <<endl;
  std::cout<< "GoodVertex     = "<< GoodVertex <<endl;
  std::cout<< "AtLeastOneJet  = "<< aJetatLessCut <<endl;
  std::cout<< "LeadingMuPtM3  = "<< LeadingMuPtM3 <<endl;
  std::cout<< "MissingET      = "<< MissingETCut <<endl;
  std::cout<< "BasicJetsCut   = "<< BasicJetsCut <<endl;
  std::cout<< "bJetsCut       = "<< bJetsCut <<endl;
  std::cout<< "MuonMetMTCut   = "<< MuonMetMTCut <<endl;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Dracarys::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we don't know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Dracarys);
