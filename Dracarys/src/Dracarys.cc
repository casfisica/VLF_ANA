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
  TriggerPath1_ = (iConfig.getParameter<vector<string>>("TriggerPathAND"));
  TriggerPath2_ = (iConfig.getParameter<vector<string>>("TriggerPathOR"));
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
  MaxMET_ = (iConfig.getParameter<double>("MaxMET"));
  //Jets
  FlagJetsAna_ = (iConfig.getParameter<bool>("FlagJetsAna"));
  FlagJetsAll_ = (iConfig.getParameter<bool>("FlagJetsAll"));
  MinNJets_ = (iConfig.getParameter<int>("MinNJets"));
  MaxNJets_ = (iConfig.getParameter<int>("MaxNJets"));
  MinJetPt_ = (iConfig.getParameter<double>("MinJetPt"));
  MaxJetEta_ = (iConfig.getParameter<double>("MaxJetEta"));
  //BJets
  FlagBJets_ = (iConfig.getParameter<bool>("FlagBJets"));
  bJetTag_ = (iConfig.getParameter<double>("bJetTag"));
  MinbJetPt_ = (iConfig.getParameter<double>("MinbJetPt"));
  MaxbJetEta_ = (iConfig.getParameter<double>("MaxbJetEta"));
  MinNbJets_ = (iConfig.getParameter<int>("MinNbJets"));
  MaxNbJets_ = (iConfig.getParameter<int>("MaxNbJets"));
  //MTMuonMET
  FlagMTMuonMET_ = (iConfig.getParameter<bool>("FlagMTMuonMET"));
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
  bool FlagPassMET = false;
  bool FlagPassJets = false;
  bool FlagPassBJets = false;
  bool FlagPassMTMuonMET =false;
  
  //////////////////////////Trigger//////////////////////////

  if (  !TriggerPath1_.empty() ||  !TriggerPath2_.empty() ){
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
      
      /*If the And Trigger Vector is not empty compare it*/
      if( !TriggerPath1_.empty() && !FlagTrigger1){
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
      
      /*If OrTriggerVector(OTV) is not empty, and the flag is not already true, then compare it*/
      if( !TriggerPath2_.empty() && !FlagTrigger2 ){
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
      if ( debug_ ) std::cout<< std::endl << "Triggers cuts PASS";
    }
  
  }else{
    FlagPassTrigg = true;
    if ( debug_ ) std::cout << "No trigger was asked" << std::endl;
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
    if ( debug_ ) std::cout<< "No vertex was asked" << std::endl;
  }
    
  //For the tighmuon
  reco::Vertex vertex = vertices->at(0);

  ///////////////////////END Vertices////////////////////////
  
  ///////////////////////////Muons///////////////////////////

  /*Handling Muons*/
  edm::Handle<edm::View<pat::Muon> > muons;
  iEvent.getByToken(tok_muons_, muons);
  
  if (!FlagMuonsAna_) FlagPassMuon=true;//For no muon cuts
   
  if ( MinNMuons_ == 0 && muons->size()==0 ) FlagPassMuon=true;

  if( (FlagMuonsAna_ || FlagMuonsAll_) &&  muons->size() > 0 ){ //Min number apply CAS
    
    //bool flagMuonChooser=false;
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
    
    int MuLooIDCountAll = 0;
    int MuMedIDCountAll = 0;
    int MuTigIDCountAll = 0;
    int contadorAll = 0;



    if (debug_) std::cout << "Muon counter loop, Muon Number: " <<  muons->size() <<std::endl;
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
      if ( MuonID_==3 ) {
	flagMuID = true;
	MuIDCount++;
      }
      ////MUON Kinematic////
      
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
	tempMuons.push_back(mu);
	tempMuon_charge.push_back(muon->charge());
	tempCombined_Iso.push_back(MuonIso);
	tempMuon_loose.push_back(muon->isLooseMuon()); 
	tempMuon_medium.push_back(muon->isMediumMuon()); 
	tempMuon_tight.push_back(muon->isTightMuon(vertex));
      }
      if (debug_) std::cout <<"OurMuonDefinitionCounter"<<OurMuonDefinitionCounter << std::endl;
      



      if( FlagMuonsAll_ ){
	
	if ( muon->isLooseMuon() ) MuLooIDCountAll++;
	if ( muon->isMediumMuon() ) MuMedIDCountAll++;
	if ( muon->isTightMuon(vertex) ) MuTigIDCountAll++;
	  
	XYZTLorentzVector mu1(muon->px(), muon->py(), muon->pz(), muon->energy());
	AllMuons.push_back(mu1);
	AllMuon_charge.push_back(muon->charge());
	AllCombined_Iso.push_back(MuonIso);
	AllMuon_loose.push_back(muon->isLooseMuon());
	AllMuon_medium.push_back(muon->isMediumMuon());
	AllMuon_tight.push_back(muon->isTightMuon(vertex));
	  
	contadorAll++;
      	
      }


    }//endfor muons
    
    if( FlagMuonsAll_ ){
      AllNMuonstight = MuTigIDCount;
      AllNMuonsmedium = MuMedIDCount;
      AllNMuonsloose = MuLooIDCount;
      
      if ( debug_ )  std::cout << "Number of allmuons: "<< contadorAll << std::endl;
    }


    if ( (OurMuonDefinitionCounter>=MinNMuons_) && (OurMuonDefinitionCounter<=MaxNMuons_)) {
      FlagPassMuon =true;
      if (debug_) std::cout <<"MuonChooser PASS" << std::endl;
    }
    
    if ( FlagPassMuon ){
      Dracarys::LeadingMuPtM3++;
      //TTree Filling
      //NMuons=OurMuonDefinitionCounter;
      //NMuons= muons->size();//The multiplicity without cuts
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
      
    }  
  } else{//End If AnaMuons
    if (debug_ && FlagMuonsAna_ ) std::cout <<"No muons in the collection" << std::endl;
    if (debug_ && !FlagMuonsAna_) std::cout <<"No muon cut was asked" << std::endl;
  }
  
  
//End if AllMuons
  
  


   /////////////////////////END Muons/////////////////////////
  
  
  /////////////////////////////MET////////////////////////////
  /*Handling MET*/   
  edm::Handle<pat::METCollection> mets;
  iEvent.getByToken(tok_met_,mets);
  const pat::MET &met = mets->front();
  
  if ( debug_ )  std::cout << "MET: "<< met.pt() << std::endl;

  if (met.pt() >= MinMET_ && met.pt() <= MaxMET_){ 
    Dracarys::MissingETCut++;
    MET = XYZTLorentzVector(met.px(), met.py(), met.pz(), met.energy());
    FlagPassMET = true;
    if ( debug_ )  std::cout << "Pass MET cuts" << std::endl;
  }

  ///////////////////////////END MET//////////////////////////



  /////////////////////////////JETS///////////////////////////
  /*Handling Jets*/
  edm::Handle<edm::View<pat::Jet> > jets;
  iEvent.getByToken(tok_jets_,jets);

  if (!FlagJetsAna_) FlagPassJets = true ;
  if (!FlagBJets_)  FlagPassBJets = true ;
   
  if (MinNJets_ == 0 && jets->size() == 0 && FlagJetsAna_ ) FlagPassJets = true ;
  
  if( (jets->size() > 0) && ( FlagJetsAna_ || FlagBJets_) ){ 
    
    for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
      
      XYZTLorentzVector je(jet->px(), jet->py(), jet->pz(), jet->energy());
      
      if (debug_) std::cout <<"Jet Pt: " << je.Pt() <<std::endl;
      
      if( (jet->pt() > MinJetPt_ ) && (abs(jet->eta()) < MaxJetEta_ ) ) {
	AnaJets.push_back(je);
	bJetDiscriminator.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
      }
      if (FlagBJets_){
	if( (jet->pt() > MinbJetPt_) && (abs(jet->eta()) < MaxbJetEta_) && (jet->bDiscriminator("combinedSecondaryVertexBJetTags")>bJetTag_) ) {
	  BJets.push_back(je);
	}
      }
    }//END JETS FOR
    
    if (debug_) std::cout <<"OurJets Multiplicity: " << (int) AnaJets.size() <<std::endl;
    if (debug_) std::cout <<"OurBJets Multiplicity: " << (int) BJets.size() <<std::endl;
	 
    
    if ( ((int) AnaJets.size() >= MinNJets_) && ( (int) AnaJets.size() <= MaxNJets_) ) FlagPassJets = true;
    if ( ( (int) BJets.size() >= MinNbJets_) && ( (int) BJets.size() <= MaxNbJets_) ) FlagPassBJets=true;
    
  }//END JETS IF

    if ( debug_ && FlagPassJets ) std::cout << "Jet cuts PASS" <<std::endl;
    if ( debug_ && FlagPassBJets ) std::cout << "BJet cuts PASS" <<std::endl;
  
  if(FlagJetsAll_ && jets->size() > 0 ){
    for(edm::View<pat::Jet>::const_iterator jet=jets->begin(); jet!=jets->end(); ++jet){
      XYZTLorentzVector je(jet->px(), jet->py(), jet->pz(), jet->energy());
      //if (debug_) std::cout <<"Jet Pt: " << je.Pt() <<std::endl;
      AllJets.push_back(je);
      bAllJetDiscriminator.push_back(jet->bDiscriminator("combinedSecondaryVertexBJetTags"));
      
    }//END ALLJETS FOR
    if (debug_) std::cout <<"AllJet Multiplicity : " << AllJets.size()  <<std::endl;
  }
  
  ///////////////////////////END JETS/////////////////////////
  
  
  ///////////////////////////MTMuonMET/////////////////////////
  
  if( !FlagMTMuonMET_ ) {
    if (debug_) std::cout << "No MTMuonMET cut was asked" << std::endl;
    FlagPassMTMuonMET = true;
  }

  if( FlagMTMuonMET_ && muons->size() > 0){
    
    TLorentzVector LeadingMuon, Met;
    bool FlagNoMuon = false;
    
    //if no cuts in muon is asked, the MTMuonMET will be calculated using the leading muon in the collection
    if( !FlagMuonsAna_ ){
      const auto& mu = muons-> at(0);
      LeadingMuon.SetPxPyPzE(mu.px(), mu.py(), mu.pz(), mu.energy());
      FlagNoMuon = false;
    }else if (AnaMuons.size() > 0 ){
      LeadingMuon.SetPxPyPzE(AnaMuons[0].px(), AnaMuons[0].py(), AnaMuons[0].pz(), AnaMuons[0].energy());
    }else{
      FlagNoMuon = true;
      FlagPassMTMuonMET = false;
      if (debug_) std::cout << "No muons in our definition to calculate MTMuonMET" << std::endl;
    }
    
    if( !FlagNoMuon ){
      Met.SetPxPyPzE(MET.px(), MET.py(), MET.pz(), MET.energy());
      MT_LeadingMuon_MET = sqrt(2*LeadingMuon.Pt()*Met.Pt()*(1-cos(LeadingMuon.DeltaPhi(Met))));
      
      if ( (MT_LeadingMuon_MET < MinMTMuonMet_) || (MT_LeadingMuon_MET > MaxMTMuonMet_) ) {
	FlagPassMTMuonMET = false;
	//MuonMetMTCut++;
      }else{
	FlagPassMTMuonMET = true;
      }
    }
    if ( debug_ && FlagPassMTMuonMET ) std::cout << "MTMuonMET cut PASS" << std::endl;
  }//END if MT
   
  
  
  /////////////////////////END MTMuonMET///////////////////////



 ////////////////////////FILLING THE TREE//////////////////
    
  if ( FlagPassTrigg  && FlagPassVertex  && FlagPassMuon && FlagPassMET && FlagPassJets && FlagPassBJets && FlagPassMTMuonMET ) {
    if ( debug_ )  std::cout << "PASS ALL THE CUTS" << std::endl;
    FlagSaveEvent = true;
  }
  
  if ( FlagSaveEvent ) {
    if ( debug_ )  std::cout << "Writting the tree" << std::endl;
    tree_->Fill();
  }
  
  if ( debug_ )  std::cout << std::endl << std::endl;
}


void 
Dracarys::Clean()
{
  //See what else should I put here!!!!!!!!
  AnaMuons.clear();
  AllMuons.clear();
  AnaJets.clear();
  AllJets.clear();
  BJets.clear();
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
  MT_LeadingMuon_MET=0;
  AllNMuonstight = 0;
  AllNMuonsmedium = 0;
  AllNMuonsloose = 0;
  AllMuon_charge.clear();
  AllCombined_Iso.clear();
  AllMuon_loose.clear();
  AllMuon_medium.clear();
  AllMuon_tight.clear();

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
  if(FlagJetsAna_) tree_->Branch("AnaJets",&AnaJets);
  if(FlagJetsAll_) tree_->Branch("AllJets",&AllJets);
  if(FlagBJets_) tree_->Branch("BJets",&BJets);
  tree_->Branch("AnaMET",&MET);
  if(FlagJetsAna_) tree_->Branch("combinedSecondaryVertexbJetDiscriminator",&bJetDiscriminator);
  if(FlagJetsAll_) tree_->Branch("AllcombinedSecondaryVertexbJetDiscriminator",&bAllJetDiscriminator);

  tree_->Branch("Combined_iso_DeltaBetaPU",&Combined_Iso);
  tree_->Branch("AllCombined_iso_DeltaBetaPU",&AllCombined_Iso);

  if(FlagMuonsAna_) tree_->Branch("AnaMuon_charge",&Muon_charge);
  if(FlagJetsAll_) tree_->Branch("AllMuon_charge",&AllMuon_charge);
  
  if(FlagMuonsAna_) tree_->Branch("AnaMuonLooseID",&Muon_loose);
  if(FlagMuonsAna_) tree_->Branch("AnaMuonMediumID",&Muon_medium);
  if(FlagMuonsAna_) tree_->Branch("AnaMuonTightID",&Muon_tight);

  if(FlagMuonsAll_) tree_->Branch("AllMuonLooseID",&AllMuon_loose);
  if(FlagMuonsAll_) tree_->Branch("AllMuonMediumID",&AllMuon_medium);
  if(FlagMuonsAll_) tree_->Branch("AllMuonTightID",&AllMuon_tight);
 
  if(FlagMTMuonMET_)tree_->Branch("MT_LeadingMuon_MET",&MT_LeadingMuon_MET);

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
  
  std::cout << endl;
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
