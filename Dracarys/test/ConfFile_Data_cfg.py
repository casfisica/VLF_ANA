import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
process.options   = cms.untracked.PSet( wantSummary = cms.untracked.bool(True) )

process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

process.load('Configuration.Geometry.GeometryRecoDB_cff')
process.load("Configuration.StandardSequences.MagneticField_38T_PostLS1_cff")

process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")
# NOTE: the pick the right global tag!
#    for Spring15 50ns MC: global tag is 'auto:run2_mc_50ns'
#    for Spring15 25ns MC: global tag is 'auto:run2_mc'
#    for Run 2 data: global tag is 'auto:run2_data'
#  as a rule, find the "auto" global tag in $CMSSW_RELEASE_BASE/src/Configuration/AlCa/python/autoCond.py
#  This auto global tag will look up the "proper" global tag
#  that is typically found in the DAS under the Configs for given dataset
#  (although it can be "overridden" by requirements of a given release)
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )

process.source = cms.Source("PoolSource",
                            
                            fileNames = cms.untracked.vstring(
        'root://eoscms.cern.ch//eos/cms/store/data/Run2016C/DoubleMuon/MINIAOD/23Sep2016-v1/80000/0044DA75-708B-E611-8F8B-008CFA1974A4.root' #for data
        
        )
                            )

process.demo = cms.EDAnalyzer('Dracarys',
                              bits = cms.InputTag("TriggerResults","","HLT"),
                              prescales = cms.InputTag("patTrigger"),
                              objects = cms.InputTag("selectedPatTrigger"),
                              vertices = cms.InputTag("offlineSlimmedPrimaryVertices"),
                              pileupInfo = cms.InputTag("slimmedAddPileupInfo"),
                              obmuon=cms.InputTag("slimmedMuons"),
                              objet=cms.InputTag("slimmedJets"),
                              obmet=cms.InputTag("slimmedMETs"),
                              #Is Data boolean
                              is_data = cms.bool(True),
                              #Activate debug option
                              debug = cms.bool(True),
                              #Trigger variables
                              FlagTrigger = cms.bool(True),#If use both then the rule to be evaluated is (TriggerPathAND and TriggerPathOR)
                              TriggerPathAND = cms.vstring("HLT_PFMET110_PFMHT110_IDTight"),#leve empty to not use a trigger
                              #TriggerPathOR = cms.vstring("HLT_DoubleMu3_PFMET50","HLT_PFMET110_PFMHT110_IDTight"), #leve empty to not use a trigger
                              #TriggerPathAND = cms.vstring(),#leve empty to not use a trigger
                              TriggerPathOR = cms.vstring(),#leve empty to not use a trigger
                              #Cuts
                              #Vertices
                              FlagVertices = cms.bool(True), #What to evaluate vertices
                              Pvtx_ndof_min   = cms.int32(4), #Vertices DOF
                              Pvtx_vtx_max  = cms.double(24.),
                              Pvtx_vtxdxy_max = cms.double(24.),
                              #Muons
                              FlagMuonsAna = cms.bool(True),#Want use muons (if is False, not muon cut will be applied)
                              FlagMuonsAll = cms.bool(True),#Want to save the full colection of muon (only the events that pass the general Cut)
                              MinMuonPt = cms.double(3.0), #Min muon pt - for all muons -
                              MaxMuonPt = cms.double(60.0), #Max muon pt - for all muons -
                              MuonIso = cms.double(0.15), #(0.15)Combined isolation with delta beta PU corrections (put 100 if do not want the cut)
                              MuonID = cms.int32(1), #0: Loose, 1: Medium, 2: Tight, 3: No Cut
                              MinNMuons = cms.int32(1), #Minimal number of muons following our definition
                              MaxNMuons = cms.int32(2), #Maximum number of muons following our defintiion
                              #MET
                              FlagMET=  cms.bool(True),#Want to use MET cuts
                              MinMET = cms.double(50.0), #Min MET
                              MaxMET = cms.double(13000.0), #Max MET
                              #Jets
                              FlagJetsAna = cms.bool(True),#Want use jets (if is False, no jets cut will be applied)
                              FlagJetsAll = cms.bool(True), #Collec all jets in the collection
                              MinNJets = cms.int32(1), #Minimal number of jets following our definition
                              MaxNJets = cms.int32(5), #Maximum number of jets following our defintion
                              MinJetPt = cms.double(30.0), #Min Jet Pt (30.0)
                              MaxJetEta = cms.double(5.0), #Max Jet Eta (5.0)
                              #BJet
                              FlagBJets = cms.bool(True),#Want use bjets
                              bJetTag = cms.double(0.8484), #b-jet ID working point
                              MinbJetPt = cms.double(30.0), #Min b Jet Pt (30.0)
                              MaxbJetEta = cms.double(2.4), #Max b Jet Eta (2.4)
                              MinNbJets = cms.int32(0), #Minimal number of Bjets following our definition
                              MaxNbJets = cms.int32(0), #Maximum number of Bjets following our defintion
                              #MTMuonMET
                              FlagMTMuonMET =  cms.bool(True),#
                              MinMTMuonMet =  cms.double(50.0),
                              MaxMTMuonMet =  cms.double(1000.0),
                              )

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("Data.root"),
                                   closeFileFast = cms.untracked.bool(True)
                                   )
#
# include bad muon filter
process.load('RecoMET.METFilters.BadPFMuonFilter_cfi')
process.BadPFMuonFilter.muons = cms.InputTag("slimmedMuons")
process.BadPFMuonFilter.PFCandidates = cms.InputTag("packedPFCandidates")

# include bad charged hadron filter
process.load('RecoMET.METFilters.BadChargedCandidateFilter_cfi')
process.BadChargedCandidateFilter.muons = cms.InputTag("slimmedMuons")
process.BadChargedCandidateFilter.PFCandidates = cms.InputTag("packedPFCandidates")

# met filters
process.load("VLF_ANA.Dracarys.AdditionalFilters_cfi")

process.p = cms.Path(process.goodVerticesFilterPAT * 
                     process.EcalDeadCellTriggerPrimitiveFilterPAT *
                     process.HBHENoiseFilterPAT * 
                     process.HBHENoiseIsoFilterPAT * 
		     process.BadPFMuonFilter *
		     process.BadChargedCandidateFilter *
                     process.demo)
