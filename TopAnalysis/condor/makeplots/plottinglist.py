
background_xsec_systematics = {
    
    'MC13TeV_SingleTbar_tW'         : 0.0537,
    'MC13TeV_SingleT_tW'            : 0.0537,
    'MC13TeV_SingleTbar_t'          : 0.051,
    'MC13TeV_SingleT_t'             : 0.0405,
    'MC13TeV_W0Jets'                : -0.057,
    'MC13TeV_W1Jets'                : 0.101,
    'MC13TeV_W2Jets'                : 0.328,
    'MC13TeV_QCDMuEnriched30to50'   : 1.0,
    'MC13TeV_QCDMuEnriched50to80'   : 1.0,
    'MC13TeV_QCDMuEnriched80to120'  : 1.0,
    'MC13TeV_QCDMuEnriched120to170' : 1.0,
    'MC13TeV_QCDMuEnriched170to300' : 1.0,
    'MC13TeV_QCDMuEnriched300to470' : 1.0,
    'MC13TeV_QCDMuEnriched470to600' : 1.0,
    'MC13TeV_QCDMuEnriched600to800' : 1.0,
    'MC13TeV_QCDMuEnriched800to1000': 1.0,
    'MC13TeV_QCDMuEnriched1000toInf': 1.0,
    'MC13TeV_QCDEMEnriched30to50'   : 1.0,
    'MC13TeV_QCDEMEnriched50to80'   : 1.0,
    'MC13TeV_QCDEMEnriched80to120'  : 1.0,
    'MC13TeV_QCDEMEnriched120to170' : 1.0,
    'MC13TeV_QCDEMEnriched170to300' : 1.0,
    'MC13TeV_QCDEMEnriched300toInf' : 1.0,
    'MC13TeV_DY50toInf_mlm'         : 0.05,
    'MC13TeV_DY10to50'              : 0.05,
    'MC13TeV_ZZTo2L2Nu'             : 0.034,
    'MC13TeV_ZZTo2L2Q'              : 0.034,
    'MC13TeV_WWToLNuQQ'             : 0.074,
    'MC13TeV_WWTo2L2Nu'             : 0.074,
    'MC13TeV_WZTo3LNu'              : 0.047,
    'MC13TeV_TTWToQQ'               : 0.162,
    'MC13TeV_TTZToQQ'               : 0.111,
    'MC13TeV_TTZToLLNuNu'           : 0.111
}

dir2 = {
    'MC13TeV_DY50toInf_mlm'           :
    [
        ['MC13TeV_DY50toInf_mlm_xsec_up', 'MC13TeV_DY50toInf_mlm_xsec_down']
    ],
    'MC13TeV_DY10to50'                :
    [
        ['MC13TeV_DY10to50_xsec_up', 'MC13TeV_DY10to50_xsec_down']
    ],
    'MC13TeV_ZZTo2L2Nu'               :
    [
        ['MC13TeV_ZZTo2L2Nu_xsec_up', 'MC13TeV_ZZTo2L2Nu_xsec_down']
    ],
    'MC13TeV_ZZTo2L2Q'               :
    [
        ['MC13TeV_ZZTo2L2Q_xsec_up', 'MC13TeV_ZZTo2L2Q_xsec_down']
    ],
    'MC13TeV_WWToLNuQQ'               :
    [
        ['MC13TeV_WWToLNuQQ_xsec_up', 'MC13TeV_WWToLNuQQ_xsec_down']
    ],
    'MC13TeV_WWTo2L2Nu'               :
    [
        ['MC13TeV_WWTo2L2Nu_xsec_up', 'MC13TeV_WWTo2L2Nu_xsec_down']
    ],
     'MC13TeV_WZTo3LNu'               :
    [
        ['MC13TeV_WZTo3LNu_xsec_up', 'MC13TeV_WZTo3LNu_xsec_down']
    ],
     'MC13TeV_TTWToQQ'               :
    [
        ['MC13TeV_TTWToQQ_xsec_up', 'MC13TeV_TTWToQQ_xsec_down']
    ],
     'MC13TeV_TTZToQQ'               :
    [
        ['MC13TeV_TTZToQQ_xsec_up', 'MC13TeV_TTZToQQ_xsec_down']
    ],
     'MC13TeV_TTZToLLNuNu'               :
    [
        ['MC13TeV_TTZToLLNuNu_xsec_up', 'MC13TeV_TTZToLLNuNu_xsec_down']
    ],
    
    'MC13TeV_QCDEMEnriched120to170'   : 
    [
        ['MC13TeV_QCDEMEnriched120to170_xsec_up',               'MC13TeV_QCDEMEnriched120to170_xsec_down']
    ],
    'MC13TeV_QCDEMEnriched170to300'   : 
    [
        ['MC13TeV_QCDEMEnriched170to300_xsec_up',               'MC13TeV_QCDEMEnriched170to300_xsec_down']
    ],
    'MC13TeV_QCDEMEnriched300toInf'   : 
    [
        ['MC13TeV_QCDEMEnriched300toInf_xsec_up',               'MC13TeV_QCDEMEnriched300toInf_xsec_down']
    ],
    'MC13TeV_QCDEMEnriched30to50'     : 
    [
        ['MC13TeV_QCDEMEnriched30to50_xsec_up',                 'MC13TeV_QCDEMEnriched30to50_xsec_down']
    ],
    'MC13TeV_QCDEMEnriched50to80'     : 
    [
        ['MC13TeV_QCDEMEnriched50to80_xsec_up',                 'MC13TeV_QCDEMEnriched50to80_xsec_down']
    ],
    'MC13TeV_QCDEMEnriched80to120'    : 
    [
        ['MC13TeV_QCDEMEnriched80to120_xsec_up',                'MC13TeV_QCDEMEnriched80to120_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched1000toInf'  : 
    [
        ['MC13TeV_QCDMuEnriched1000toInf_xsec_up',              'MC13TeV_QCDMuEnriched1000toInf_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched120to170'   : 
    [
        ['MC13TeV_QCDMuEnriched120to170_xsec_up',               'MC13TeV_QCDMuEnriched120to170_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched170to300'   :
    [
        ['MC13TeV_QCDMuEnriched170to300_xsec_up',               'MC13TeV_QCDMuEnriched170to300_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched300to470'   :
    [
        ['MC13TeV_QCDMuEnriched300to470_xsec_up',               'MC13TeV_QCDMuEnriched300to470_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched30to50'     :
    [
        ['MC13TeV_QCDMuEnriched30to50_xsec_up',                 'MC13TeV_QCDMuEnriched30to50_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched470to600'   :
    [
        ['MC13TeV_QCDMuEnriched470to600_xsec_up',               'MC13TeV_QCDMuEnriched470to600_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched50to80'     :
    [
        ['MC13TeV_QCDMuEnriched50to80_xsec_up',                 'MC13TeV_QCDMuEnriched50to80_xsec_down']

    ],
    'MC13TeV_QCDMuEnriched600to800'   :
    [
        ['MC13TeV_QCDMuEnriched600to800_xsec_up',               'MC13TeV_QCDMuEnriched600to800_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched800to1000'  :
    [
        ['MC13TeV_QCDMuEnriched800to1000_xsec_up',              'MC13TeV_QCDMuEnriched800to1000_xsec_down']
    ],
    'MC13TeV_QCDMuEnriched80to120'    :
    [
        ['MC13TeV_QCDMuEnriched80to120_xsec_up',                'MC13TeV_QCDMuEnriched80to120_xsec_down']
    ],
    'MC13TeV_SingleT_t'               :
    [
        ['MC13TeV_SingleT_t_xsec_up',                           'MC13TeV_SingleT_t_xsec_down']
    ],
    'MC13TeV_SingleT_tW'              :
    [
        ['MC13TeV_SingleT_tW_xsec_up',                          'MC13TeV_SingleT_tW_xsec_down'],
        ['MC13TeV_SingleT_tW_m175v5',                           'MC13TeV_SingleT_tW_m169v5'],
        ['MC13TeV_SingleT_tW_isrup',                            'MC13TeV_SingleT_tW_isrdn'],
        ['MC13TeV_SingleT_tW_fsrup',                            'MC13TeV_SingleT_tW_fsrdn'],
        ['MC13TeV_SingleT_tW_meup',                             'MC13TeV_SingleT_tW_medn']
    ],
    'MC13TeV_SingleTbar_t'            :
    [
        ['MC13TeV_SingleTbar_t_xsec_up',                        'MC13TeV_SingleTbar_t_xsec_down']
    ],
    'MC13TeV_SingleTbar_tW'           :
    [
        ['MC13TeV_SingleTbar_tW_xsec_up',                       'MC13TeV_SingleTbar_tW_xsec_down'],
        ['MC13TeV_SingleTbar_tW_m175v5',                        'MC13TeV_SingleTbar_tW_m169v5'],
        ['MC13TeV_SingleTbar_tW_isrup',                         'MC13TeV_SingleTbar_tW_isrdn'],
        ['MC13TeV_SingleTbar_tW_fsrup',                         'MC13TeV_SingleTbar_tW_fsrdn'],
        ['MC13TeV_SingleTbar_tW_meup',                          'MC13TeV_SingleTbar_tW_medn']
    ],
    'MC13TeV_TTJets'                  :                  
    [
        ['MC13TeV_TTJets_pileup_up',                            'MC13TeV_TTJets_pileup_down'],
        ['MC13TeV_TTJets_trig_efficiency_correction_up',        'MC13TeV_TTJets_trig_efficiency_correction_down'],
        ['MC13TeV_TTJets_sel_efficiency_correction_up',         'MC13TeV_TTJets_sel_efficiency_correction_down'],
        ['MC13TeV_TTJets_b_fragmentation_up',                   'MC13TeV_TTJets_b_fragmentation_down'],
        ['MC13TeV_TTJets_semilep_BR_up',                        'MC13TeV_TTJets_semilep_BR_down'],
        ['MC13TeV_TTJets_jec_CorrelationGroupMPFInSitu_up',     'MC13TeV_TTJets_jec_CorrelationGroupMPFInSitu_down'],
        ['MC13TeV_TTJets_jec_RelativeFSR_up',                   'MC13TeV_TTJets_jec_RelativeFSR_down'],
        ['MC13TeV_TTJets_jec_CorrelationGroupUncorrelated_up',  'MC13TeV_TTJets_jec_CorrelationGroupUncorrelated_down'],
        ['MC13TeV_TTJets_jec_FlavorPureGluon_up',               'MC13TeV_TTJets_jec_FlavorPureGluon_down'],
        ['MC13TeV_TTJets_jec_FlavorPureQuark_up',               'MC13TeV_TTJets_jec_FlavorPureQuark_down'],
        ['MC13TeV_TTJets_jec_FlavorPureCharm_up',               'MC13TeV_TTJets_jec_FlavorPureCharm_down'],
        ['MC13TeV_TTJets_jec_FlavorPureBottom_up',              'MC13TeV_TTJets_jec_FlavorPureBottom_down'],
        ['MC13TeV_TTJets_jer_up',                               'MC13TeV_TTJets_jer_down'],
        ['MC13TeV_TTJets_btag_heavy_up',                        'MC13TeV_TTJets_btag_heavy_down'],
        ['MC13TeV_TTJets_btag_light_up',                        'MC13TeV_TTJets_btag_light_down'],
        ['MC13TeV_TTJets_csv_heavy_up',                         'MC13TeV_TTJets_csv_heavy_down'],
        ['MC13TeV_TTJets_csv_light_up',                         'MC13TeV_TTJets_csv_light_down'],
        ['MC13TeV_TTJets_tracking_up',                          'MC13TeV_TTJets_tracking_down'],
        ['MC13TeV_TTJets_m173v5',                               'MC13TeV_TTJets_m171v5'],
        ['MC13TeV_TTJets_isrup',                                'MC13TeV_TTJets_isrdn'],
        ['MC13TeV_TTJets_fsrup',                                'MC13TeV_TTJets_fsrdn'],
        ['MC13TeV_TTJets_hdampup',                              'MC13TeV_TTJets_hdampdn'],
        ['MC13TeV_TTJets_ueup',                                 'MC13TeV_TTJets_uedn']
    ],
    'MC13TeV_W0Jets'                  :
    [
        ['MC13TeV_W0Jets_xsec_up', 'MC13TeV_W0Jets_xsec_down']
    ],
    'MC13TeV_W1Jets'                  :
    [
        ['MC13TeV_W1Jets_xsec_up', 'MC13TeV_W1Jets_xsec_down']
    ],
    'MC13TeV_W2Jets'                  :
    [
        ['MC13TeV_W2Jets_xsec_up', 'MC13TeV_W2Jets_xsec_down']
    ]
} 

dir1 =  {
    'MC13TeV_SingleT_tW'              :
    ['MC13TeV_SingleT_tW_DS'],
    'MC13TeV_SingleTbar_tW'           :
    ['MC13TeV_SingleTbar_tW_DS'],
    'MC13TeV_TTJets'                  :
    ['MC13TeV_TTJets_Peterson_frag',
     'MC13TeV_TTJets_id1002muR1muF2hdampmt272.7225',
     'MC13TeV_TTJets_id1003muR1muF0.5hdampmt272.7225',
     'MC13TeV_TTJets_id1004muR2muF1hdampmt272.7225',
     'MC13TeV_TTJets_id1005muR2muF2hdampmt272.7225', 
     'MC13TeV_TTJets_id1007muR0.5muF1hdampmt272.7225',
     'MC13TeV_TTJets_id1009muR0.5muF0.5hdampmt272.7225',
     'MC13TeV_TTJets_evtgen',
     'MC13TeV_TTJets_herwig',
     'MC13TeV_TTJets_cflip',
     'MC13TeV_TTJets_erdON',
     'MC13TeV_TTJets_qcdBased',
     'MC13TeV_TTJets_gluonMove']
}
