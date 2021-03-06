import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.FwkReport.reportEvery = 100
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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(50000) )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(
                                #'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/DYJetsToLL_M-50_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v2/120000/02A210D6-F5C3-E611-B570-008CFA197BD4.root'
                                'root://cms-xrd-global.cern.ch//store/mc/RunIISummer16MiniAODv2/WJetsToLNu_HT-100To200_TuneCUETP8M1_13TeV-madgraphMLM-pythia8/MINIAODSIM/PUMoriond17_80X_mcRun2_asymptotic_2016_TrancheIV_v6_ext1-v1/120000/14473FEF-1ACD-E611-8C84-00266CFFBC60.root'
        
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
                              is_data = cms.bool(False),
                              #Activate debug option
                              debug = cms.bool(False),
                              #Trigger variables
                              isTrigger = cms.bool(True),
                              isTriggerToo = cms.bool(False),
                              TriggerPath1 = cms.string("HLT_PFMET110_PFMHT110_IDTight"),
#                              TriggerPath1 = cms.string("HLT_TkMu20"),
                              TriggerPath2 = cms.string("HLT_DoubleMu3_PFMET50"),
#                              TriggerPath2 = cms.string("HLT_Mu17_TrkIsoVVL"),
                              #Cuts
                              Pvtx_ndof_min   = cms.int32(4), #Vertices DOF
                              Pvtx_vtx_max  = cms.double(24.),
                              Pvtx_vtxdxy_max = cms.double(24.),
                              MinMuonPt = cms.double(3.0), #Min muon pt - for all muons - 3.0
                              MaxMuonPt = cms.double(40.0), #Max muon pt - for all muons -40.0
                              MuonIso = cms.double(0.15), #Combined isolation with delta beta PU corrections
                              MuonIsoCut = cms.bool(True), # If the cut in the Iso is taken into account
                              MuonID = cms.int32(1), #0: Loose, 1: Medium, 2: Tight
                              MuonIDCut = cms.bool(True), #If the cut in the ID is taken into account
                              MinNMuons = cms.int32(1), #Minimal number of muons following our definition
                              MaxNMuons = cms.int32(1), #Maximum number of muons following our defintiion
                              MinMET = cms.double(110.0), #Min MET 50
                              MinJetPt = cms.double(30.0), #Min Jet Pt 30
                              MaxJetEta = cms.double(5.0), #Max Jet Eta
                              bJetTag = cms.double(0.8484), #b-jet ID working point
                              MinbJetPt = cms.double(30.0), #Min b Jet Pt
                              MaxbJetEta = cms.double(2.4), #Max b Jet Eta
                              MinNJets = cms.int32(1), #Minimal number of jets following our definition
                              MaxNJets = cms.int32(10), #Maximum number of jets following our defintion
                              MinNbJets = cms.int32(0), #Minimal number of jets following our definition
                              MaxNbJets = cms.int32(0), #Maximum number of jets following our defintion
                              MinMTMuonMet =  cms.double(0.0),
                              MaxMTMuonMet =  cms.double(100.0), #100
)
process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string("DY.root"),
                                   closeFileFast = cms.untracked.bool(True)
)

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
