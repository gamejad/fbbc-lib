#include "TFile.h"
#include "TDirectory.h"
#include "TH1.h"
#include "TH1D.h"
#include "TH1F.h"
#include "TH2D.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include <cmath>
#include <random>

#include "../src/fbbc-lib.h"
#include "../src/fbbc-lib.cpp"

const Double_t c = 0.299792458; // speed of light, mm/ps

void example(){
    //for Z vertex gambling
    std::random_device rd{};
    std::mt19937 gen{rd()};
    std::normal_distribution<> d{0,250};
    
    vector<double> plt_coords = {-850, -630, -500, 500, 630, 850}; // MCP plates coordinates, mm
    const int NumOfPlates = plt_coords.size();
    double rin = 15; // inner radius of MCP plate rings, mm
    double rout = 25; // outer radius of MCP plate rings, mm
    int rad_sec_num = 8;
    int ang_sec_num = 8;
    double eff = 0.9; // MCP efficiency
    double t_prec = 50; //MCP readout time precision, ps

    FBBCDetector detector(plt_coords, rin, rout, rad_sec_num, ang_sec_num, eff, t_prec); // FBBC detector object

    for(const auto i : detector.GetPlatesPseudorapidity()) // write MCP plate's covered pseudorapidity regions
        cout << i[0] << ' ' << i[1] << endl;

    vector <TH1D *> fHistTimeDistr(NumOfPlates);        //create time distribution histograms for particles,
    for (int i = 0; i < NumOfPlates; i++)                //registered by different MCPs
        fHistTimeDistr[i] = new TH1D (Form("Time distr. %d, %3.1f mm", i, plt_coords.at(i)),
                                      Form("Time distr. %d, %3.1f mm; T, ps; counts", i, plt_coords.at(i)), 2000, 0, 40000);
    TH1D *fHistZ = new TH1D("Z dist.", "Z dist.", 20, -50, 50);

    TFile file("Particles.root");     // output data of MC simulation in .root format
    TTree* tree = (TTree*)file.Get("particles");

    const int NMaxTrack = 15000;

    Int_t ntracks;
    Double_t impact;
    Double_t p0[NMaxTrack];
    Double_t px[NMaxTrack];
    Double_t py[NMaxTrack];
    Double_t pz[NMaxTrack]; //
    Int_t pdgcode[NMaxTrack];
    Int_t charge[NMaxTrack];

    tree->SetBranchAddress("npart",    &ntracks );
    tree->SetBranchAddress("impact_b", &impact  );
    tree->SetBranchAddress("p0",       p0       );
    tree->SetBranchAddress("px",       px       );
    tree->SetBranchAddress("py",       py       );
    tree->SetBranchAddress("pz",       pz       );
    tree->SetBranchAddress("pdgcode",  pdgcode  );
    tree->SetBranchAddress("charge",   charge   );

    int nEvents = tree->GetEntries();
    cout<< "\nNEvents = " << nEvents << endl;


    // EVENT LOOP STARTS HERE :
    cout << "Start of Event loop:\n";

    for ( int iEvent = 0; iEvent < nEvents; iEvent++)
    {
        tree->GetEntry(iEvent);
        if ( iEvent % 1 == 0 ) cout << "processing  " << (int)iEvent << "\r"; cout.flush();

        double z_bias = d(gen); // Z coord. of event, mm
        fHistZ->Fill(z_bias);
        vector<ParticleFBBC> particles_fbbc; // vector of particles in event, converted into ParticleFBBC format

        for ( int iTrack = 0; iTrack < ntracks; iTrack++ )
        {
            if(charge[iTrack] == 0) continue;
            double p = sqrt(px[iTrack]*px[iTrack]+py[iTrack]*py[iTrack]+pz[iTrack]*pz[iTrack]);
            double phi = atan2(py[iTrack], px[iTrack]);
            ParticleFBBC part = {iTrack, p0[iTrack], p/c, pz[iTrack]/c, phi, z_bias }; // create particle in ParticleFBBC format
            particles_fbbc.push_back(part);                                        // we devide momenta on speed of light
        }                                                                           // to convert it from GeV to GeV/c

        detector.SetParticlesFBBC(particles_fbbc); // Set the event to our FBBC detector
        auto plates_out = detector.GetOutputVector(); //output of FBBC detector processing in format
                                                        // vector<vector<pair<double, int>>> (plates<particles<[time, id]>>)
        for(int i = 0; i < NumOfPlates; i++)
        {   for(const auto [t, id] : plates_out[i]) // fill time distr. histos for all MCP in detector
                fHistTimeDistr[i]->Fill(t);
        }

    }
    
    cout << endl;

    TFile *fileOutput = new TFile( "example_output.root", "RECREATE" ); //write our histos in output file
    for(const auto hist : fHistTimeDistr)
        hist->Write();
    fHistZ->Write();
    fileOutput->Close();
    
    TCanvas *c1 = new TCanvas("c1","Canvas Example",200,10,1600,1200);
    c1->Divide(2,1);
    c1->cd(1);
    fHistTimeDistr.at(0)->Draw();
    c1->cd(2);
    fHistZ->Draw();



}
