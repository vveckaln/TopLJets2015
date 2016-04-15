import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import os
import sys
options = VarParsing.VarParsing ('standard')

options.register('nevents', 1000, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Number of events")
options.register('yoda', 'run.yoda', VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.string, "YODA output file")
options.register('uselhewgt', True, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Use LHE weight")
options.register('lhewgt', 0, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.int, "Number of LHE weight used for this run")
options.register('warnings', False, VarParsing.VarParsing.multiplicity.singleton,VarParsing.VarParsing.varType.bool, "Show warnings")

# define the syntax for parsing
# you need to enter in the cfg file:
# search for arguments entered after cmsRun
if( hasattr(sys, "argv") ):
    # split arguments by comma - seperating different variables
    for args in sys.argv :
        arg = args.split(',')
        # split further by = to separate variable name and value
        for val in arg:
            val = val.split('=')
            # set variable var to value val (expected crab syntax: var=val)
            if(len(val)==2):
                setattr(options,val[0], val[1])

print options

process = cms.Process('GEN')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.MessageLogger.cerr.FwkReport.reportEvery = 1000
if not options.warnings: process.MessageLogger.cerr.threshold = 'ERROR'
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('Configuration.StandardSequences.Generator_cff')
process.load('IOMC.EventVertexGenerators.VtxSmearedNominalCollision2015_cfi')
process.load('GeneratorInterface.Core.genFilterSummary_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.nevents)
)

# Input source
process.source = cms.Source("PoolSource",
    dropDescendantsOfDroppedBranches = cms.untracked.bool(False),
    fileNames = cms.untracked.vstring(
#        '/store/mc/RunIIWinter15wmLHE/TTJets_13TeV-amcatnloFXFX-pythia8/LHE/MCRUN2_71_V1-v1/10000/00CA9EB1-C8C7-E411-B40E-02163E00F34B.root' ),
#'/store/mc/RunIISpring15DR74/TTJets_TuneCUETP8M1_13TeV-amcatnloFXFX-pythia8/AODSIM/Asympt50ns_MCRUN2_74_V9A-v1/00000/006A97CE-D301-E511-8072-0025905A60A6.root' ),
#'/store/mc/Summer12WMLHE/TTJets_TuneCUETP8M1_8TeV-amcatnloFXFX-pythia8/GEN/START53_V7C-v5/00000/00CA251D-0D78-E411-95CE-00259075D714.root'),
'/store/generator/Summer12WMLHE/TT_weights_CT10_8TeV-powheg/GEN/START53_V7C-v2/00000/00F54BA9-3464-E411-9966-0025905B85F6.root'),
    inputCommands = cms.untracked.vstring('keep *', 
        'drop LHEXMLStringProduct_*_*_*'),
    secondaryFileNames = cms.untracked.vstring()
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('Configuration/GenProduction/python/TOP-RunIIWinter15GS-00001-fragment.py nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)

# Output definition

process.RECOSIMoutput = cms.OutputModule("PoolOutputModule",
    SelectEvents = cms.untracked.PSet(
        SelectEvents = cms.vstring('generation_step')
    ),
    dataset = cms.untracked.PSet(
        dataTier = cms.untracked.string(''),
        filterName = cms.untracked.string('')
    ),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    fileName = cms.untracked.string('TOP-RunIIWinter15GS-00001-fragment_py_GEN.root'),
    outputCommands = process.RECOSIMEventContent.outputCommands,
    splitLevel = cms.untracked.int32(0)
)

# Additional output definition

# Other statements
process.genstepfilter.triggerConditions=cms.vstring("generation_step")
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:mc', '')


process.load("GeneratorInterface.RivetInterface.rivetAnalyzer_cfi")
process.rivetAnalyzer.AnalysisNames     = cms.vstring(
                                            #'MC_GENERIC',
                                            #'MC_TTBAR_HADRON',
                                            #'CMS_2015_I1370682',
                                            #'CMS_2015_I1370682_internal',
                                            #'CMS_TOP_12_042',
                                            #'CMS_LesHouches2015',
                                            #'CMS_MC_BRAD',
                                            'ATLAS_2014_I1304289',
                                            'MC_TTbar_TruthSel',
                                          )
process.rivetAnalyzer.OutputFile        = options.yoda
process.rivetAnalyzer.UseExternalWeight = True
process.rivetAnalyzer.useLHEweights     = options.uselhewgt
process.rivetAnalyzer.LHEweightNumber   = options.lhewgt
process.rivetAnalyzer.LHECollection     = cms.InputTag('externalLHEProducer')

from Configuration.Generator.HerwigppDefaults_cfi import *
from Configuration.Generator.HerwigppUE_EE_5C_cfi import *
from Configuration.Generator.HerwigppPDF_CT10_LO_cfi import *                                                                  # Import CTEQ6L PDF as shower pdf
from Configuration.Generator.HerwigppEnergy_8TeV_cfi import *
from Configuration.Generator.HerwigppLHEFile_cfi import *
from Configuration.Generator.HerwigppMECorrections_cfi import *

process.generator = cms.EDFilter("ThePEGHadronizerFilter",
        herwigDefaultsBlock,
        herwigppUESettingsBlock,
        herwigppPDFSettingsBlock,
        herwigppEnergySettingsBlock,
        herwigppLHEFileSettingsBlock,
        herwigppMECorrectionsSettingsBlock,

        configFiles = cms.vstring(),
        parameterSets = cms.vstring(
                'hwpp_cmsDefaults',
                'hwpp_ue_EE5C',
                'hwpp_cm_8TeV',
                'hwpp_pdf_CT10',                     # Shower PDF matching with the tune
                'hwpp_LHE_Powheg',                    # Showering LHE files from MadGraph5_aMC@NLO. Use the same PDF for the shower as for the hard subprocess afore
                'hwpp_MECorr_Off',                      # Switch off ME corrections while showering LHE files as recommended by Herwig++ authors
        ),

        crossSection = cms.untracked.double(-1),
        filterEfficiency = cms.untracked.double(1.0),
)

# Path and EndPath definitions
process.generation_step = cms.Path(process.generator*process.rivetAnalyzer)
process.genfiltersummary_step = cms.EndPath(process.genFilterSummary)
process.endjob_step = cms.EndPath(process.endOfProcess)

# Schedule definition
process.schedule = cms.Schedule(process.generation_step,process.genfiltersummary_step,process.endjob_step)
# filter all path with the production filter sequence
for path in process.paths:
	getattr(process,path)._seq = process.generator * getattr(process,path)._seq 


