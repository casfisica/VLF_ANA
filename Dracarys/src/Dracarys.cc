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
// Original Author:  Jose Ruiz
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
  FlagTrigger_AND_ = (iConfig.getParameter<bool>("FlagTrigger_AND"));
  FlagTrigger_OR_ = (iConfig.getParameter<bool>("FlagTrigger_OR"));
  TriggerPath1_ = (iConfig.getParameter<vector<string>>("TriggerPath1"));
  TriggerPath2_ = (iConfig.getParameter<vector<string>>("TriggerPath2"));
				    //Debugging option
  debug_ = (iConfig.getParameter<bool>("debug"));
  //Cuts
  //Vertices
  FlagVertices_ = (iConfig.getParameter<bool>("FlagVertices"));
  Pvtx_ndof_min_ = (iConfig.getParameter<int>("Pvtx_ndof_min"));
  Pvtx_vtx_max_ = (iConfig.getParameter<double>("Pvtx_vtx_max"));
  Pvtx_vtxdxy_max_ = (iConfig.getParameter<double>("Pvtx_vtxdxy_max"));
  //Muons
  FlagMuonsAna_ = (iConfig.getParameter<bool>("FlagMuonsAna"));
  FlagMuonsAll_ = (iConfig.getParameter<bool>("FlagMuonsAll"));
  MinMuonPt_ = (iConfig.getParameter<double>("MinMuonPt"));
  MaxMuonPt_ = (iConfig.getParameter<double>("MaxMuonPt"));
  MuonIso_ = (iConfig.getParameter<double>("MuonIso"));
  MuonID_ = (iConfig.getParameter<int>("MuonID"));
  MinNMuons_ = (iConfig.getParameter<int>("MinNMuons"));
  MaxNMuons_ = (iConfig.getParameter<int>("MaxNMuons"));
  //MET
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
//!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
//Poner una manera de que corra en un orden espesífico (con un while y un contador para los numeros del orden con un if orden=contador ejecuta)
  
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


  //Counting events
  Dracarys::NoCuts++;
  
  //Cleaning all variables
  Clean();

  //Flags
  bool FlagSaveEvent = false;
  bool FlagPassTrigg = false;
  bool FlagPassVertex = false;
  bool FlagPassMuon = false;
  
  //////////////////////////Trigger//////////////////////////

  if ( FlagTrigger_OR_ || FlagTrigger_AND_ ){
    edm::Handle<edm::TriggerResults> triggerBits;
    edm::Handle<pat::TriggerObjectStandAlone> triggerObjects;
    edm::Handle<pat::PackedTriggerPrescales> triggerPrescales;
    
    iEvent.getByToken(triggerBits_, triggerBits);
    iEvent.getByToken(triggerObjects_, triggerObjects);
    iEvent.getByToken(triggerPrescales_, triggerPrescales);
    
    bool FlagTrigger1=false;
    bool FlagTrigger2=false;
    
    
    //To have any number of triggers
    std::vector<std::string> TriggerNamesVector;
    
    const edm::TriggerNames &names = iEvent.triggerNames(*triggerBits);
    
    int NAndTrig = TriggerPath1_.size();
    int NAndGood = 0;
    
    if ( debug_ )  std::cout<< std::endl << std::endl << "Triggers found: ";
    
    for (unsigned int i = 0, n = triggerBits->size(); i < n; ++i){
      if (!triggerBits->accept(i)) continue; //If the event does not pass trigger continue
      
      /*Cut the version of the trigger path*/
      std::string nameV = names.triggerName(i);
      std::string version ="_v";
      size_t start_pos = nameV.rfind(version);
      if (start_pos != std::string::npos){
	std::string TriggerNameVersionOff = nameV.erase(start_pos, version.length()+4);//Delete 4 chars includong _v (carefull with names like cas_varie!!!)
      }
      
      /*If the flag is true, and the And Trigger Vector is not empty compare it*/
      if( FlagTrigger_AND_ && !TriggerPath1_.empty() && !FlagTrigger1){
	for ( auto trig:TriggerPath1_ ){
	  if( nameV == trig ) {
	    NAndGood++;
	    if ( debug_ )  std::cout << nameV << " ,";
	  }
	  if ( NAndTrig == NAndGood ) {
	    FlagTrigger1=true;
	    break;
	  }
	}      
      }else {
	FlagTrigger1=true;
      }
      
      /*If the flag to use the OrTriggerVector(OTV) is true, and the OTV is not empty, and the flag is not already true, then compare it*/
      if( FlagTrigger_OR_ && !TriggerPath2_.empty() && !FlagTrigger2 ){
	for ( auto trig:TriggerPath2_ ){
	  if( nameV == trig ) {
	    FlagTrigger2=true;
	    if ( debug_ )  std::cout << nameV << " ,";
	    break;
	  }	
	}      
      }else {
	FlagTrigger2=true;
      }
      
    }//End For over Events triggers
    
    
    if( FlagTrigger1 && FlagTrigger2 ){
      FlagPassTrigg = true;
      //Counting number of events which pass the triggers
      Dracarys::TriggerPathCut++;
      if ( debug_ )  {
	if ( !FlagTrigger1 ) std::cout<< "No And Triggers by Flag" << std::endl;
	if ( !FlagTrigger2 ) std::cout<< "No Or Triggers by Flag" << std::endl;
	std::cout<< std::endl << "      !!!Event Pass Triggers!!!";
      }
    }
  }
  ////////////////////////END Trigger////////////////////////
  
  
  /////////////////////////Vertices//////////////////////////
  edm::Handle<reco::VertexCollection> vertices;
  iEvent.getByToken(tok_vertex_, vertices);

  if (FlagVertices_){
    
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
    
    if (flagGoodVertex) {
      Dracarys::GoodVertex++;
      FlagPassVertex = true;
      if ( debug_ ) std::cout<< "Number of Good Vertices: "<< Nvertices << std::endl;
    }
  }else{
    FlagPassVertex = true;
    if ( debug_ ) std::cout<< "No vertex is asked" << std::endl;
  }
    
  //For the tighmuon
  reco::Vertex vertex = vertices->at(0);

  ///////////////////////END Vertices////////////////////////
  
  ///////////////////////////Muons///////////////////////////

  /*Handling Muons*/
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(tok_muons_, muons);
  
  if (!FlagMuonsAna_) FlagPassMuon=true;//For no muon cuts   
  if( (muons->size() > 0) && FlagMuonsAna_ ){ //More than 0 and the flag
    
    int OurMuonDefinitionCounter=0;
    int MuTigIDCount = 0;
    int MuMedIDCount = 0;
    int MuLooIDCount = 0;
    int MuISOCount = 0;
    int MuIDCount =0;
    bool MuonChooser = false;
    
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
	if ( MuonID_==3 ) {
	  flagMuID = true;
	}
	
      }
      ////MUON Kinematic////
      // Could be the eta, phy, ...
      if( (muon->pt() > MinMuonPt_) && (muon->pt() < MaxMuonPt_) ) {
	flagPtCut = true;
      }else{
	if (debug_) std::cout <<"Kinematics Fail"<< std::endl;
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
	AnaMuon.push_back(mu);
	Muon_charge.push_back(muon->charge());
	Combined_Iso.push_back(MuonIso);
	Muon_loose.push_back(muon->isLooseMuon());
	Muon_medium.push_back(muon->isMediumMuon());
	Muon_tight.push_back(muon->isTightMuon(vertex));
	MuonChooser = true;
      }
      
    }//End For AnaMuons

    if ( MuonChooser ){
      NMuonstight = MuTigIDCount;
      NMuonsmedium = MuMedIDCount;
      NMuonsloose = MuLooIDCount;
      NMuonsIso = MuISOCount;
      NMuonsID = MuIDCount;
    }
    
    if (debug_) std::cout <<"OurMuonDefinitionCounter"<<OurMuonDefinitionCounter << std::endl;
    
    if ( (OurMuonDefinitionCounter>=MinNMuons_) && (OurMuonDefinitionCounter<=MaxNMuons_)) {
      FlagPassMuon = true;
      if (debug_) std::cout <<"PASS Muon Cuts" << std::endl;
    }
    	
    
  }//End If AnaMuons

   if( (muons->size() > 0) && FlagMuonsAll_ ){

     int MuLooIDCount = 0;
     int MuMedIDCount = 0;
     int MuTigIDCount = 0;
     
     for(edm::View<pat::Muon>::const_iterator muon=muons->begin(); muon!=muons->end(); ++muon){
       
       double MuonIso = (muon->pfIsolationR04().sumChargedHadronPt	\
			 +  max(0., muon->pfIsolationR04().sumNeutralHadronEt + muon->pfIsolationR04().sumPhotonEt - \
				0.5*muon->pfIsolationR04().sumPUPt))/muon->pt();
       
       if ( muon->isLooseMuon() ) MuLooIDCount++;
       if ( muon->isMediumMuon() ) MuMedIDCount++;
       if ( muon->isTightMuon(vertex) ) MuTigIDCount++;

       XYZTLorentzVector mu(muon->px(), muon->py(), muon->pz(), muon->energy());
       AllMuon.push_back(mu);
       AllMuon_charge.push_back(muon->charge());
       AllCombined_Iso.push_back(MuonIso);
       AllMuon_loose.push_back(muon->isLooseMuon());
       AllMuon_medium.push_back(muon->isMediumMuon());
       AllMuon_tight.push_back(muon->isTightMuon(vertex));
       
     }
     
     AllNMuonstight = MuTigIDCount;
     AllNMuonsmedium = MuMedIDCount;
     AllNMuonsloose = MuLooIDCount;
       
     
   }//End if AllMuons
  



   /////////////////////////END Muons/////////////////////////
   
   ////////////////////////FILLING THE TREE//////////////////

   /*************************Muons***************************/
   if ( FlagPassMuon ){
     Dracarys::LeadingMuPtM3++;
     //NMuons=OurMuonDefinitionCounter;
     
    }

   ///FlagSaveEvent Como lo conmuto!!


   if ( FlagSaveEvent ) tree_->Fill();
   
  if ( debug_ )  std::cout << std::endl << std::endl;
}


void 
Dracarys::Clean()
{
  //See what else should I put here!!!!!!!!
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
}

// ------------ method called once each job just before starting event loop  ------------
void Dracarys::beginJob()
{
  //Cuts by a map< std::string, int > name of the cut, and counter.
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
  if(FlagMuonsAna_) tree_->Branch("AnaMuons",&AnaMuons);
  if(FlagMuonsAll_) tree_->Branch("AllMuons",&AllMuons);
  tree_->Branch("AnaJets",&AnaJets);
  tree_->Branch("AnaMET",&MET);
  tree_->Branch("combinedSecondaryVertexbJetDiscriminator",&bJetDiscriminator);

  tree_->Branch("Combined_iso_DeltaBetaPU",&Combined_Iso);
  tree_->Branch("AllCombined_iso_DeltaBetaPU",&AllCombined_Iso);

  tree_->Branch("AnaMuon_charge",&Muon_charge);
  tree_->Branch("AllMuon_charge",&AllMuon_charge);
  
  tree_->Branch("AnaMuonLooseID",&Muon_loose);
  tree_->Branch("AnaMuonMediumID",&Muon_medium);
  tree_->Branch("AnaMuonTightID",&Muon_tight);

  tree_->Branch("AllMuonLooseID",&AllMuon_loose);
  tree_->Branch("AllMuonMediumID",&AllMuon_medium);
  tree_->Branch("AllMuonTightID",&AllMuon_tight);


  tree_->Branch("AnaJetsMultiplicity",&NJets);
  tree_->Branch("AnabJetsMultiplicity",&NbJets);
  tree_->Branch("MT_LeadingMuon_MET",&MT_LeadingMuon_MET);
  tree_->Branch("Vertices",&Nvertices);
  tree_->Branch("InTimePU",&NObservedInTimePUVertices);
  tree_->Branch("TruePU",&NTruePUInteractions);


  
  if ( debug_ ){
    if ( FlagTrigger_AND_ && !TriggerPath1_.empty()){
      std::cout << "============================================"<< std::endl;
      std::cout << "                   Trigger AND"<< std::endl;
      std::cout << "============================================"<< std::endl;
      
      for ( auto trig:TriggerPath1_ ){
	std::cout << trig << ", " << std::endl;
      }
    }else{
      std::cout << "============================================"<< std::endl;
      std::cout << "                NO Trigger AND"<< std::endl;
      std::cout << "============================================"<< std::endl;
    }

    if( FlagTrigger_OR_ && !TriggerPath2_.empty()){
      std::cout << "============================================"<< std::endl;
      std::cout << "                   Trigger OR"<< std::endl;
      std::cout << "============================================"<< std::endl;
      
      for ( auto trig:TriggerPath2_ ){
	std::cout << trig << ", " << std::endl;
      }
    }else{
      std::cout << "============================================"<< std::endl;
      std::cout << "                NO Trigger OR"<< std::endl;
      std::cout << "============================================"<< std::endl;
    }
  }
  


}

// ------------ method called once each job just after ending the event loop  ------------
void Dracarys::endJob() 
{
  std::cout<< "NoCuts= "<< NoCuts <<endl;
  
  std::cout<< "TriggerPathCut= "<< TriggerPathCut <<endl;
  std::cout<< "GoodVertex= "<< GoodVertex <<endl;
  std::cout<< "AtLeastOneJet= "<< aJetatLessCut <<endl;
  std::cout<< "LeadingMuPtM3= "<< LeadingMuPtM3 <<endl;
  std::cout<< "MissingET= "<< MissingETCut <<endl;
  std::cout<< "BasicJetsCut= "<< BasicJetsCut <<endl;
  std::cout<< "bJetsCut= "<< bJetsCut <<endl;
  std::cout<< "MuonMetMTCut= "<< MuonMetMTCut <<endl;
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
