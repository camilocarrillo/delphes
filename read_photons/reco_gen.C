#include "TH1.h"
#include "TSystem.h"

//------------------------------------------------------------------------------

struct MyPlots
{
    TH1 *fJetPT[2];
    TH1 *fPhotonPT[2];
    TH1 *fPhotonPTResolution[2];
    TH1 *fMissingET;
    TH1 *fAngle;
    TH1 *fElectronPT;
    TH1 *fMassResolution;
};

//------------------------------------------------------------------------------

class ExRootResult;
class ExRootTreeReader;

//------------------------------------------------------------------------------

void BookHistograms(ExRootResult *result, MyPlots *plots)
{
    THStack *stack;
    TLegend *legend;
    TPaveText *comment;
    
    // book 2 histograms for PT of 1st and 2nd leading jets
    
    plots->fJetPT[0] = result->AddHist1D("jet_pt_0", "leading jet P_{T}","jet P_{T}, GeV/c", "number of jets", 50, 0.0, 500.0);    

    plots->fJetPT[1] = result->AddHist1D("jet_pt_1", "2nd leading jet P_{T}","jet P_{T}, GeV/c", "number of jets", 50, 0.0, 500.0);

    plots->fJetPT[0]->SetLineColor(kRed);
    plots->fJetPT[1]->SetLineColor(kBlue);
    
    // book 1 stack of 2 histograms
    
    stack = result->AddHistStack("jet_pt_all", "1st and 2nd jets P_{T}");
    stack->Add(plots->fJetPT[0]);
    stack->Add(plots->fJetPT[1]);
    
  // book 2 histograms for PT of 1st and 2nd leading photons
    
    plots->fPhotonPT[0] = result->AddHist1D("photon_pt_0", "leading photon P_{T}","photon P_{T}, GeV/c", "number of photons",100, 0.0, 150.0);
    
    plots->fPhotonPT[1] = result->AddHist1D("photon_pt_1", "2nd leading photon P_{T}","photon P_{T}, GeV/c", "number of photons",100, 0.0, 100.0);
    
    plots->fPhotonPTResolution[0] = result->AddHist1D("Delta_photon_pt_0", "leading photon Delta P_{T}","#Delta photon P_{T}/P_{T}^{true}", "number of photons",100, -1, 1.0);
    
    plots->fPhotonPTResolution[1] = result->AddHist1D("Delta_photon_pt_1", "trailing photon Delta P_{T}","#Delta photon P_{T}/P_{T}^{true}", "number of photons",100, -1.0, 1.0);
    plots->fPhotonPT[0]->SetLineColor(kRed);
    plots->fPhotonPT[1]->SetLineColor(kBlue);
    
    plots->fPhotonPTResolution[0]->SetLineColor(kRed);
    plots->fPhotonPTResolution[1]->SetLineColor(kBlue);
    
    // book 1 stack of 2 histograms
    
    stack = result->AddHistStack("photon_pt_all", "1st and 2nd photons P_{T}");
    stack->Add(plots->fPhotonPT[0]);
    stack->Add(plots->fPhotonPT[1]);
    
    // book legend for stack of 2 histograms
    
    legend = result->AddLegend(0.72, 0.56, 0.98, 0.68);
    legend->AddEntry(plots->fPhotonPT[0], "leading photon", "l");
    legend->AddEntry(plots->fPhotonPT[1], "second photon", "l");
    
    // attach legend to stack (legend will be printed over stack in .eps file)
    
    result->Attach(stack, legend);
    
    // book 1 stack of 2 histograms
    
    stack = result->AddHistStack("Delta photon_pt_all", "1st and 2nd photons Delta P_{T}");
    stack->Add(plots->fPhotonPTResolution[0]);
    stack->Add(plots->fPhotonPTResolution[1]);
    
    result->Attach(stack, legend);
    
    // book more histograms
    
    plots->fElectronPT = result->AddHist1D("electron_pt", "electron P_{T}","electron P_{T}, GeV/c", "number of electrons",50, 0.0, 100.0);

    plots->fMassResolution = result->AddHist1D("Delta_m_gg", "M_{#gamma#gamma}","#Delta m_{#gamma#gamma}/m_{#gamma#gamma}^{true}", "number of events",100, -1.0, 1.0);
    
    plots->fMissingET = result->AddHist1D("missing_et", "Missing E_{T}","Missing E_{T}, GeV", "number of events",60, 0.0, 500.0);
    
    plots->fAngle = result->AddHist1D("angle", "Angle MET-{#gamma#gamma}","Angle MET-{#gamma#gamma}", "number of events",60, 0.0, 3.14159265);
    
    // book general comment

    comment = result->AddComment(0.64, 0.56, 0.98, 0.68);
    comment->AddText("H->#gamma#gamma");
    
    
    // attach comment to single histograms
    
    result->Attach(plots->fJetPT[0], comment);
    result->Attach(plots->fJetPT[1], comment);
    result->Attach(plots->fPhotonPT[0], comment);
    result->Attach(plots->fPhotonPT[1], comment);
    result->Attach(plots->fPhotonPTResolution[0], comment);
    result->Attach(plots->fPhotonPTResolution[1], comment);
    result->Attach(plots->fElectronPT, comment);
    result->Attach(plots->fMassResolution, comment);
    plots->fMissingET->SetStats();
    plots->fMassResolution->SetStats();
    plots->fPhotonPTResolution[0]->SetStats();
    plots->fPhotonPTResolution[1]->SetStats();
}

//------------------------------------------------------------------------------

void AnalyseEvents(ExRootTreeReader *treeReader, MyPlots *plots)
{
    TClonesArray *branchJet = treeReader->UseBranch("Jet");
    TClonesArray *branchPhoton = treeReader->UseBranch("Photon");
    TClonesArray *branchElectron = treeReader->UseBranch("Electron");
    TClonesArray *branchMissingET = treeReader->UseBranch("MissingET");
    TClonesArray *branchParticle = treeReader->UseBranch("Particle");
    
    
    Long64_t allEntries = treeReader->GetEntries();

    cout << "** Chain contains " << allEntries << " events" << endl;
    
    Jet *jet[2];
    Photon *photon[2];
    MissingET *met;
    Electron *electron;
    
    Long64_t entry;
    
    Int_t i;
    
    // Loop over all events
    //    for(entry = 0; entry < allEntries; ++entry)
    for(entry = 0; entry < 10000; ++entry){
	// Load selected branches with data from specified event
	treeReader->ReadEntry(entry);
	
	GenParticle *particle, *daughter[2];
	Int_t index, motherPID;
	
	float higgs_mass=0;
	float lead_pt_gen=-1;	    
	float trail_pt_gen=-1;
	
	//cout<<"Loop on the gen particles..."<<endl;
	for(index=0;index<branchParticle->GetEntriesFast();index++){
	    particle = (GenParticle *) branchParticle->At(index);
	    if(particle->PID==25){
		higgs_mass=particle->Mass;
		daughter[0] = (GenParticle *) branchParticle->At(particle->D1);
		daughter[1] = (GenParticle *) branchParticle->At(particle->D2);
		lead_pt_gen=daughter[0]->PT;
		trail_pt_gen=daughter[1]->PT;
	    }
	}
	
	float buff;
	
	if(lead_pt_gen<trail_pt_gen){
	    buff=lead_pt_gen;
	    lead_pt_gen=trail_pt_gen;
	    trail_pt_gen=buff;
	}
	
	//cout<<"Loop on gen particles finished"<<endl;
	//cout<<"higgs mass="<<higgs_mass<<endl;
	//cout<<"lead_pt_gen="<<lead_pt_gen<<endl;
	//cout<<"trail_pt_gen="<<trail_pt_gen<<endl;
	
	// Analyse two leading jets
	if(branchJet->GetEntriesFast() >= 2){
	    jet[0] = (Jet*) branchJet->At(0);
	    jet[1] = (Jet*) branchJet->At(1);
	    
	    plots->fJetPT[0]->Fill(jet[0]->PT);
	    plots->fJetPT[1]->Fill(jet[1]->PT);
	}
	
	// Analyse two leading photons
	if(branchPhoton->GetEntriesFast() >= 2){
	    //cout<<"There are "<<branchPhoton->GetEntriesFast()<<" photons in the event"<<endl;
	    for(Int_t ip=0;ip<branchPhoton->GetEntriesFast();ip++){		
		Photon *thisphoton =(Photon*) branchPhoton->At(ip);
		//cout<<"\t"<<ip<<" "<<thisphoton->PT<<" "<<endl;
	    }
	    
	    photon[0] = (Photon*) branchPhoton->At(0);
	    photon[1] = (Photon*) branchPhoton->At(1);
	    
	    plots->fPhotonPT[0]->Fill(photon[0]->PT);
	    plots->fPhotonPT[1]->Fill(photon[1]->PT);
	    
	    plots->fPhotonPTResolution[0]->Fill((lead_pt_gen-photon[0]->PT)/lead_pt_gen);
	    plots->fPhotonPTResolution[1]->Fill((trail_pt_gen-photon[1]->PT)/trail_pt_gen);
	    
	    //cout<<"diphoton mass "<<((photon[0]->P4())+(photon[1]->P4())).M()<<endl;
	    
	    plots->fMassResolution->Fill((higgs_mass- ((photon[0]->P4())+(photon[1]->P4())).M())/higgs_mass);

	    //Double_t a = v1.Angle(v2.Vect()); 

	    if(branchMissingET->GetEntriesFast() > 0){
		met = (MissingET*) branchMissingET->At(0);
		Double_t angle=((photon[0]->P4())+(photon[1]->P4())).Angle((met->P4()).Vect());
		//cout<<"angle "<<angle<<endl;
		plots->fAngle->Fill(angle);
	    }
	    
	    //cout<<"(pTLead,pTTrail)="<<photon[0]->PT<<" "<<photon[1]->PT<<endl;
	    //cout<<"(pTLead,pTTrail)="<<photon[0]->motherPID<<endl;
	}
	
	// Analyse missing ET
	if(branchMissingET->GetEntriesFast() > 0){
	    met = (MissingET*) branchMissingET->At(0);
	    plots->fMissingET->Fill(met->MET);
	}
	
	// Loop over all electrons in event
	for(i = 0; i < branchElectron->GetEntriesFast(); ++i){
	    electron = (Electron*) branchElectron->At(i);
	    plots->fElectronPT->Fill(electron->PT);
	}
    }
}

//------------------------------------------------------------------------------

void PrintHistograms(ExRootResult *result, MyPlots *plots)
{
  result->Print("png");
}

//------------------------------------------------------------------------------

void reco_gen(const char *inputFile)
{
  gSystem->Load("libDelphes");

  TChain *chain = new TChain("Delphes");
  chain->Add(inputFile);

  ExRootTreeReader *treeReader = new ExRootTreeReader(chain);
  ExRootResult *result = new ExRootResult();

  MyPlots *plots = new MyPlots;

  BookHistograms(result, plots);

  AnalyseEvents(treeReader, plots);

  PrintHistograms(result, plots);

  result->Write("results.root");

  cout << "** Exiting..." << endl;

  delete plots;
  delete result;
  delete treeReader;
  delete chain;
}

//------------------------------------------------------------------------------
