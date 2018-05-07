import FWCore.ParameterSet.Config as cms
######GETTING ALL ROOT FILES#########                                                                                                                                                          
import commands as cmd
##################################### 

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

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(5000) )

######GETTING ALL ROOT FILES#########
SignalContent=cmd.getoutput('ls /eos/user/j/jruizalv/VLF_Samples/MINIAODSIM/')
FullList=SignalContent.split('\n')
AllFiles=[]
for i in FullList:
    if 'root' not in i: continue
    else:
        AllFiles.append('file:/eos/user/j/jruizalv/VLF_Samples/MINIAODSIM/'+i)
print 'Files found:', AllFiles[0]
#####################################

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring(AllFiles)
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
TriggerPath2 = cms.string("HLT_DoubleMu3_PFMET50"),
                              #Cuts
Pvtx_ndof_min   = cms.int32(4), #Vertices DOF
Pvtx_vtx_max  = cms.double(24.),
                              Pvtx_vtxdxy_max = cms.double(24.),
                              MinMuonPt = cms.double(3.0), #Min muon pt - for all muons -
                              MaxMuonPt = cms.double(60.0), #Max muon pt - for all muons -
                              MuonIso = cms.double(0.15), #Combined isolation with delta beta PU corrections
                              MuonID = cms.int32(1), #0: Loose, 1: Medium, 2: Tight
                              MinNMuons = cms.int32(1), #Minimal number of muons following our definition
                              MaxNMuons = cms.int32(2), #Maximum number of muons following our defintiion
                              MinMET = cms.double(50.0), #Min MET
                              MinJetPt = cms.double(30.0), #Min Jet Pt
                              MaxJetEta = cms.double(5.0), #Max Jet Eta
                              bJetTag = cms.double(0.8484), #b-jet ID working point
                              MinbJetPt = cms.double(30.0), #Min b Jet Pt
                              MaxbJetEta = cms.double(2.4), #Max b Jet Eta
                              MinNJets = cms.int32(1), #Minimal number of jets following our definition
                              MaxNJets = cms.int32(6), #Maximum number of jets following our defintion
                              MinNbJets = cms.int32(0), #Minimal number of jets following our definition
                              MaxNbJets = cms.int32(0), #Maximum number of jets following our defintion
                              MinMTMuonMet =  cms.double(0.0),
                              MaxMTMuonMet =  cms.double(100.0),
                              )
process.TFileService = cms.Service("TFileService",
fileName = cms.string("Signal.root"),
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
