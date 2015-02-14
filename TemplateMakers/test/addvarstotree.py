from ROOT import TFile, TChain, TTree, TH1, gSystem, TLorentzVector, gROOT
import numpy as n
from variables import *
import sys
from time import sleep
import os
import glob

gSystem.Load('libttH-13TeVMultiLeptonsTemplateMakers.so')
gROOT.SetBatch(1)

infilestring = sys.argv[2]
outfilestring = sys.argv[3]
#infilestring += '*.root'


print infilestring
print outfilestring


#infile = TFile( 'test_200evts.root' )
filetest = False
# while True:
# 
# 	infile = TFile(infilestring)
# 	
# 	if infile.IsZombie():
# 		print "it's a zombie"
# 		infile.Close()
# 		filetest = True
# 		sleep(3)
# 	else:
# 		break

tree = TChain('OSTwoLepAna/summaryTree')
exitcount = 0

# ## This helps to negotiate an NDT3/Condor-specific issue:
# while True:
# 
# 	#fileexists = os.path.exists(infilestring)
# 	fileexists = len(glob.glob(infilestring))
# 	
# 	if not fileexists:
# 		sleep(5)
# 		print "doesnt exist"
# 		filetest = True
# 		exitcount += 1
# 		if exitcount is 20:
# 			print "been at this for 1 minute; exiting..."
# 			sys.exit(1)
# 	else:
# 		#infile = TFile(infilestring)
# 		tree.Add(infilestring)
# 		break
# 
# if filetest:
# 	print "made it out!"

tree.Add(infilestring)
#tree = infile.Get('OSTwoLepAna/summaryTree')	
	
entries = tree.GetEntries()
print " "
print entries

#copiedfile = TFile( 'newfile.root' , 'RECREATE' )
#copiedfile = TFile( 'ttH125_ss_v30_all.root' , 'RECREATE' )
copiedfile = TFile.Open( outfilestring , 'RECREATE', 'test' )
newtree = tree.CloneTree(0)

#testthing = {0.,0.,0.,0.,} ROOT.ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double> >
#testthing = TLorentzVector(0.,0.,0.,0.)

############################################################
## if the tree has these branches already, overwrite them:


SumPt_handle = n.array([-99],dtype=float)
MHT_handle = n.array([-99],dtype=float)
SumJetMass_handle = n.array([-99],dtype=float)
SumNonTaggedJetMass_handle = n.array([-99],dtype=float)
SumJetPt_handle = n.array([-99],dtype=float)
AvgBtagDiscNonBtags_handle = n.array([-99],dtype=float)
AvgBtagDiscBtags_handle = n.array([-99],dtype=float)
MinDrJets_handle = n.array([-99],dtype=float)

NumHiggsLikeDijet15_handle = n.zeros(1,dtype=float)
HiggsLikeDijetMass2_handle = n.array([-99],dtype=float)
HiggsLikeDijetMass_handle = n.array([-99],dtype=float)
GenHiggsDijetMass_handle = n.array([-99],dtype=float)
DeltaPhiMetLepLep_handle = n.array([-99],dtype=float)
DeltaRMetLepLep_handle = n.array([-99],dtype=float)
MTMetLepLep_handle = n.array([-99],dtype=float)
MassMetLepLep_handle = n.array([-99],dtype=float)
MaxDeltaPhiMetJet_fromHiggs_handle = n.array([-99],dtype=float)
MinDeltaPhiMetJet_fromHiggs_handle = n.array([-99],dtype=float)
MaxDeltaPhiMetJet_handle = n.array([-99],dtype=float)
MinDeltaPhiMetJet_handle = n.array([-99],dtype=float)
DeltaPhiMetLep2_handle = n.array([-99],dtype=float)
DeltaPhiMetLep1_handle = n.array([-99],dtype=float)
DeltaRJets_FromHiggs_handle = n.array([-99],dtype=float)
DeltaPhiJets_FromHiggs_handle = n.array([-99],dtype=float)
WLikeDijetMass81_handle = n.array([-99],dtype=float)
DeltaPhiLepLep_handle = n.array([-99],dtype=float)
DeltaRLepLep_handle = n.array([-99],dtype=float)
Zmass_handle = n.array([-99],dtype=float)
MassLepLep_handle = n.array([-99],dtype=float)
numLooseBJets_handle = n.zeros(1,dtype=float)
numMediumBJets_handle = n.zeros(1,dtype=float)

newtree.SetBranchAddress("MHT", MHT_handle)
newtree.SetBranchAddress("SumJetMass", SumJetMass_handle)
newtree.SetBranchAddress("SumNonTaggedJetMass", SumNonTaggedJetMass_handle)
newtree.SetBranchAddress("SumPt", SumPt_handle)
newtree.SetBranchAddress("SumJetPt", SumJetPt_handle)
newtree.SetBranchAddress("AvgBtagDiscNonBtags", AvgBtagDiscNonBtags_handle)
newtree.SetBranchAddress("AvgBtagDiscBtags", AvgBtagDiscBtags_handle)
newtree.SetBranchAddress("MinDrJets", MinDrJets_handle)
newtree.SetBranchAddress("NumHiggsLikeDijet15", NumHiggsLikeDijet15_handle)
newtree.SetBranchAddress("HiggsLikeDijetMass2", HiggsLikeDijetMass2_handle)
newtree.SetBranchAddress("HiggsLikeDijetMass", HiggsLikeDijetMass_handle)
newtree.SetBranchAddress("GenHiggsDijetMass", GenHiggsDijetMass_handle)
newtree.SetBranchAddress("DeltaPhiMetLepLep", DeltaPhiMetLepLep_handle)
newtree.SetBranchAddress("DeltaRMetLepLep", DeltaRMetLepLep_handle)
newtree.SetBranchAddress("MTMetLepLep", MTMetLepLep_handle)
newtree.SetBranchAddress("MassMetLepLep", MassMetLepLep_handle)
#newtree.SetBranchAddress("MaxDeltaPhiMetJet", MaxDeltaPhiMetJet_handle)
#newtree.SetBranchAddress("MinDeltaPhiMetJet", MinDeltaPhiMetJet_handle)
newtree.SetBranchAddress("MaxDeltaPhiMetJet", MaxDeltaPhiMetJet_handle)
newtree.SetBranchAddress("MinDeltaPhiMetJet", MinDeltaPhiMetJet_handle)
newtree.SetBranchAddress("DeltaPhiMetLep2", DeltaPhiMetLep2_handle)
newtree.SetBranchAddress("DeltaPhiMetLep1", DeltaPhiMetLep1_handle)
newtree.SetBranchAddress("DeltaRJets_FromHiggs", DeltaRJets_FromHiggs_handle)
newtree.SetBranchAddress("DeltaPhiJets_FromHiggs", DeltaPhiJets_FromHiggs_handle)
newtree.SetBranchAddress("WLikeDijetMass81", WLikeDijetMass81_handle)
newtree.SetBranchAddress("DeltaPhiLepLep", DeltaPhiLepLep_handle)
newtree.SetBranchAddress("DeltaRLepLep", DeltaRLepLep_handle)
newtree.SetBranchAddress("Zmass", Zmass_handle)
newtree.SetBranchAddress("MassLepLep", MassLepLep_handle)
newtree.Branch("numLooseBJets", numLooseBJets_handle, "numLooseBJets/D")
newtree.Branch("numMediumBJets", numMediumBJets_handle, "numMediumBJets/D")		  


############################################################	  
## (re)calculate vars from variables.py:
			  
for entry in tree:
	preselelectrons = entry.preselected_electrons
	looseelectrons = entry.loose_electrons
	tightelectrons = entry.tight_electrons
	
	preselmuons = entry.preselected_muons
	loosemuons = entry.loose_muons
	tightmuons = entry.tight_muons
	
	preselleptons = entry.preselected_leptons
	looseleptons = entry.loose_leptons
	tightleptons = entry.tight_leptons
		
	preseljets = entry.preselected_jets
	loosebtags = entry.loose_bJets
	met = entry.met
	
	######################
	electrons = looseelectrons
	muons = loosemuons
	leptons = looseleptons
	jets = preseljets	

	if (leptons.size()<2):
		continue
	
#	for electron in thing:
#		eleTLV = electron.obj
#		print eleTLV.Px()
	#sumjethandle[0] = 5.0
	#entry.SumJetPt = 5.0	
	
	
	# need the [0] for branch handles !!!!
	SumJetPt_handle[0] = getsumpt(jets)	# A.K.A. 'HT'
	AvgBtagDiscNonBtags_handle[0] = getAvgCSV(jets,'M',False)   
	AvgBtagDiscBtags_handle[0] = getAvgCSV(jets,'M',True)
	MinDrJets_handle[0] = getTwoObjKineExtreme(jets,jets,'min','dR')
	SumPt_handle[0] = getsumpt(jets, electrons, muons)	
	
	objs_for_mht = getsumTLV(leptons,jets)
	MHT_handle[0] = objs_for_mht.Pt()
	
	taggedjets = keepTagged(jets,'M')
	nontaggedjets = keepUnTagged(jets,'M')
	numMediumBJets_handle[0] = len(taggedjets)
	numLooseBJets_handle[0] = len(loosebtags)
	
	SumNonTaggedJetMass_handle[0] = getsumpt(taggedjets)
	HiggsLikeDijetMass2_handle[0] = pickFromSortedTwoObjKine(jets,jets,'mass', 2, 125.)
	HiggsLikeDijetMass_handle[0] = 	pickFromSortedTwoObjKine(jets,jets,'mass', 1, 125.)
	NumHiggsLikeDijet15_handle[0] = getNumTwoObjKineInRange(jets,jets,'mass',125.,15.)


#newtree.SetBranchAddress("GenHiggsDijetMass", GenHiggsDijetMass_handle)

## three obj kine:

#newtree.SetBranchAddress("DeltaPhiMetLepLep", DeltaPhiMetLepLep_handle)
#newtree.SetBranchAddress("DeltaRMetLepLep", DeltaRMetLepLep_handle)
#newtree.SetBranchAddress("MTMetLepLep", MTMetLepLep_handle)
#newtree.SetBranchAddress("MassMetLepLep", MassMetLepLep_handle)	
	
	
	MaxDeltaPhiMetJet_handle[0] = 	getTwoObjKineExtreme(met,jets,'max','dPhi')
	MinDeltaPhiMetJet_handle[0] = 	getTwoObjKineExtreme(met,jets,'min','dPhi')
	DeltaPhiMetLep2_handle[0] = 	pickFromSortedTwoObjKine(met,leptons, 'dPhi', 2, 0.0) ## probably wrong
	DeltaPhiMetLep1_handle[0] = 	pickFromSortedTwoObjKine(met,leptons, 'dPhi', 1, 0.0) ## probably wrong

	
#newtree.SetBranchAddress("DeltaRJets_FromHiggs", DeltaRJets_FromHiggs_handle)
#newtree.SetBranchAddress("DeltaPhiJets_FromHiggs", DeltaPhiJets_FromHiggs_handle)
	

	WLikeDijetMass81_handle[0] = 	pickFromSortedTwoObjKine(jets,jets,'mass',1,81.)
	DeltaPhiLepLep_handle[0] = 		getTwoObjKineExtreme(leptons,leptons,'min','dPhi')
	DeltaRLepLep_handle[0] = 		getTwoObjKineExtreme(leptons,leptons,'min','dR')	
	Zmass_handle[0] = 				pickFromSortedTwoObjKine(leptons,leptons,'mass',1,91.2) #trying loose
	MassLepLep_handle[0] = 			getTwoObjKineExtreme(leptons,leptons,'min','mass')
	
	
	
	newtree.Fill()

############################################################
## save new values to tree:

newtree.Write()
#infile.Close()
copiedfile.Close()
print '  '
print 'done'
		
############################################################	
############################################################
		

		
#		
#
#
#  TwoObjectKinematic<BNleptonCollection,BNjetCollection> myMHT("pt", "vector_sum", "mht",
#                                                               &(selectedCollections.tightLooseLeptonCollection), "leptons_by_pt", 1, 99,
#                                                               &(selectedCollections.jetCollection), "jets_by_pt", 1, 99);
#
#


#		int HiggsDecayType_intree;		
#		int checkTrig_intree;
#		
#
##		double LepTrigCorr_intree;
##		double LepIDAndIsoSF_intree;
##		double TopPt_intree;
##		double TopPtCorr_intree;
##		double PU_intree;
##		double PUCorr_intree;
#		
#		//numLooseBJets;
#		//numMediumBJets;
#		//numTightBJets;
#		
#		int numJets_fromHiggs_30_intree;
#		int numJets_fromHiggs_intree;
