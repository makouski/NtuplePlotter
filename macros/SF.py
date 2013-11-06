
luminosity = 19.63*1000
gSF = luminosity

# first line: Nevt from AODSIM
# second line: number of events from gg files
#Top_num        = 6923750
# ???
TTJets1l_num        = 10629902 # checked.  was off
TTJets2l_num        = 12019013 # checked.  was off

TTgamma_num      = 71598 # checked

#WJets_num      = 57709905
#WJets1_num      = 30106028 # skim
#WJets2_num      = 27603876 # skim
WJets1_num      = 27185540 # checked. new skim (???)
WJets2_num      = 30524364 # checked. new skim (???)

#ZJets_num      = 30459503
ZJets_num      = 30459504 # ?  checked. (skimmed)

#ZZ_num            = 9799908 #AODSIM inclusive?
# 2e2mu + 2e2tau + 2mu2tau + 4e + 4mu + 4tau
ZZ_2e2mu_num       = 1497445
ZZ_2e2tau_num      = 823911
ZZ_2mu2tau_num     = 823922
ZZ_4e_num          = 1499093
ZZ_4mu_num         = 1499064
ZZ_4tau_num        = 824466

#WW_num            = 10000431 #AODSIM
WW_2l2nu_num       = 1903235 # checked

#WZ_num            = 10000283 #AODSIM
#WZ_num            = 5233969 # old
WZ_3lnu_num       = 2017979 # checked
WZ_2l2q_num       = 3215990 # checked

Wgamma_num        = 4802358 # checked

WWgamma_num       = 304285 # checked

Zgamma_num        = 6588161 # checked

TTW_num           = 196046 # checked

TTZ_num           = 210160 # checked

WHIZARD_num       = 1069486 # new full sim sample 756000 # checked

#SingToptW_num     = 497658
SingToptW_num     = 497658 # checked

#SingTopbartW_num  = 493460
SingTopbartW_num  = 493460 # checked

#SingTopT_num      = 3758227
SingTopT_num      = 99876 # checked (too few!!!!!!)

#SingTopbarT_num   = 1935072
SingTopbarT_num   = 1885072 # checked (smaller!)

#SingTopS_num      = 259961
SingTopS_num      = 259961 # checked

#SingTopbarS_num   = 139974
SingTopbarS_num   = 139974 # checked

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/QGC
Zgamma_xs         = 132.6 # PREP
Wgamma_xs         = 461.6 # PREP
WWgamma_xs        = 0.528 # PREP

# https://twiki.cern.ch/twiki/bin/viewauth/CMS/StandardModelCrossSectionsat8TeV

#Top_xs            = 227 #CMS measurement
TTJets1l_xs       = 99.44 # 227*0.676*(1-0.676)*2
TTJets2l_xs       = 23.83 # 227*(1-0.676)*(1-0.676)
TTgamma_xs        = 0.9081 * 2  # https://twiki.cern.ch/twiki/bin/view/CMS/WhizardMCTeeTeeGamma
WJets_xs          = 35640.0 # CMS measurement. was: 36257.0  # 
ZJets_xs          = 3350.0 # CMS measurement. was: 3533.0 # 3503.71 #

#WZ_xs             = 33.21 # inclusive
WZ_3lnu_xs        = 0.8674 # PREP
WZ_2l2q_xs        = 1.755  # PREP

#ZZ_xs             = 8.06 # inclusive
ZZ_2e2mu_xs       = 0.018 #??? 0.1767 # PREP
ZZ_2e2tau_xs      = 0.018 #??? 0.1767 # PREP
ZZ_2mu2tau_xs     = 0.018 #??? 0.1767 # PREP
ZZ_4e_xs          = 0.009 #??? 0.07691 # PREP
ZZ_4mu_xs         = 0.009 #??? 0.07691 # PREP
ZZ_4tau_xs        = 0.009 #??? 0.07691 # PREP

# NLO? 5.995 # 4.7 PREP this is 2l2nu, inclusive is 54.8
# 69.9 CMS measurement
# 69.9*(1-0.676)*(1-0.676)
WW_2l2nu_xs       = 7.3
TTW_xs            = 0.232 #
TTZ_xs            = 0.208 #

SingToptW_xs      = 11.1 #
SingTopbartW_xs   = 11.1 #
SingTopT_xs       = 56.4 #
SingTopbarT_xs    = 30.7 #
SingTopS_xs       = 3.79 #
SingTopbarS_xs    = 1.76 #

