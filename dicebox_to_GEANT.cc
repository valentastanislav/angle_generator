#include <iostream> // for cin & cout - requires "namespace std"
#include <fstream> // for ifstream - requires "namespace std"
#include <sstream>
#include <string>
#include <vector>
#include <dirent.h> 
#include <stdio.h> 
#include <TFile.h>
#include <TNtuple.h>
#include <TTree.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TMath.h>
#include <Math/IFunction.h>
#include <TProfile.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "Fk.h"

using namespace std;

int main(int argc, char **argv) {
  if (argc!=3){
    cout << endl;
    cout << "need directory as 1st argument" << endl;
    cout << "need cascade file name as 2nd argument" << endl;
    exit(-1);
  }
  DIR *d;
  struct dirent *dir;
  char* path = argv[1];
  d = opendir(path);
  if (!d) {
    cout << "directory " << dir->d_name << " can not be opened" << endl;
    exit(-2);
  }

  std::ifstream file(Form("%s/%s",argv[1],argv[2]),ios::in);
  if(file.fail()){
    cout << "file " << argv[1] << "/" << argv[2] << " not found" << endl;
    exit(-3);
  }
  std::ofstream outfile(Form("%s/%s.out",argv[1],argv[2]),ios::out);
  if(!outfile.is_open()){
    cout << "unable to open output file " << argv[1] << "/" << argv[2] << ".out" << endl;
    exit(-4);
  }  

  gStyle->SetPadTickX(1);  gStyle->SetPadTickY(1);
  const double RTG_K=(0.0028224-0.0003); //Zn(0.0097-0.0012) //Nb (0.006539-0.00065);//Mn (0.010367-0.0011432);//Ga (0.007113-0.00072); //Fe (0.008979-0.00095); //Cu 
  int iline=0;
  vector<double> temp, Eg, Ee, Epair, Spin, Delta, IC, Angle;
  vector<int> L;
  TH1F *h_g_sums, *h_gammas;
  string line;
  float Ak,Bk;
  TF1 *W = new TF1("W(theta)","1+[0]*ROOT::Math::legendre(2,TMath::Cos(x))+[1]*ROOT::Math::legendre(4,TMath::Cos(x))",-TMath::Pi(),+TMath::Pi());
  TFile *root_file = new TFile(Form("%s/%s.root",argv[1],argv[2]),"update");

  // Tree and its branch declaration
  TTree *cascades = new TTree("cascades","cascades");
  cascades->Branch("gammas",&Eg);
  cascades->Branch("electrons",&Ee);
  cascades->Branch("pairs",&Epair);
  cascades->Branch("spins",&Spin);
  cascades->Branch("delta",&Delta);
  cascades->Branch("angle",&Angle);

  cout << "reading file " << argv[1] << "/"<< argv[2] << endl;
  while (getline(file,line)){
    ++iline;
      // if (iline%10000==0) cout << iline << endl;
      // cout << "in " << " EVENTS.S" << ireal+1 << ".R001 on line " << iline << endl;
    stringstream ss(line);
    double n;
    ss >> n;
    while (ss.good()){
      temp.push_back(n);
      ss >> n;
    }
    temp.push_back(n);      //IC_type !0 = gamma, 1 = K-shell, 2 = higher-shell, 3 = pair
    if(iline%5==0){ //3,6,9,...
      int multiplicity=temp.at(2);
      Spin.push_back(temp.at(0));
      for(int itemp=0; itemp<multiplicity; ++itemp){
        IC.push_back(temp.at(4+2*multiplicity+itemp));
        if(IC.at(itemp)==0){ //transition is gamma
	        // cout << (temp.at(4+itemp-1)-temp.at(4+itemp)) << " ";
	        Eg.push_back(temp.at(4+itemp-1)-temp.at(4+itemp));
          Spin.push_back(temp.at(4+multiplicity+itemp));
          Delta.push_back(temp.at(4+3*multiplicity+itemp));
        }
        else if(IC.at(itemp)==3) { //transition is pair
	        Epair.push_back(temp.at(4+itemp-1)-temp.at(4+itemp)); //situation complicated, probably can be approximated by a pair of 511 keV gammas (from the positron anihilation)
          Eg.push_back(.511); Eg.push_back(.511);
          Spin.push_back(temp.at(4+multiplicity+itemp));
          Delta.push_back(temp.at(4+3*multiplicity+itemp));
        }
        else if(IC.at(itemp)==2){ //transition is higher-shell electron
	        // cout << (temp.at(4+itemp-1)-temp.at(4+itemp)) << " ";
	        Ee.push_back(temp.at(4+itemp-1)-temp.at(4+itemp));
          Spin.push_back(temp.at(4+multiplicity+itemp));
          Delta.push_back(temp.at(4+3*multiplicity+itemp));
        }
        else if(IC.at(itemp)==1){ //transition is K-shell electron
	        // cout << (temp.at(4+itemp-1)-temp.at(4+itemp)) << " ";
	        Ee.push_back(temp.at(4+itemp-1)-temp.at(4+itemp)-RTG_K); //write Ee-RTG to Ee
          Eg.push_back(RTG_K); //write RTG to Eg
          Spin.push_back(temp.at(4+multiplicity+itemp));
          Delta.push_back(temp.at(4+3*multiplicity+itemp));
        }
        else{
          cout << "WTF? what is this internal conversion label?" << endl;
        }
      }
      if((multiplicity==1)&&(IC.at(0)==3)){
        Angle.push_back(TMath::Pi());
      }
      else if((multiplicity==2)&&(IC.at(0)==3)&&(IC.at(1)==3)){
        Angle.push_back(TMath::Pi());
        Angle.push_back(-1e3);
        Angle.push_back(TMath::Pi());
      }
      else if((multiplicity==2)&&(IC.at(0)==3)){
        Angle.push_back(TMath::Pi());
        Angle.push_back(-1e3);
      }
      else if((multiplicity==2)&&(IC.at(1)==3)){
        Angle.push_back(-1e3);
        Angle.push_back(TMath::Pi());
      }
      else if(multiplicity>1){
        for(int itemp=0; itemp<(multiplicity-1); ++itemp){
          if((IC.at(itemp)==0)&&(IC.at(itemp+1)==0)){
            L.push_back(abs(Spin.at(itemp)-Spin.at(itemp+1))); //first multipolarity
            L.push_back(abs(Spin.at(itemp+1)-Spin.at(itemp+2))); //second multipolarity
            // cout << "got M=2 cascade with L1=" << L.at(0) << " L2=" << L.at(1) << " made by spins " << Spin.at(itemp) << " -> " << Spin.at(itemp+1) << " -> " << Spin.at(itemp+2) << endl;
            // cout << Fk(2,L.at(0),L.at(0),Spin.at(itemp),Spin.at(itemp+1)) << endl;

            //set parameters for TF1 W(theta)
            Bk = (Fk(2,L.at(0),L.at(0),Spin.at(itemp),Spin.at(itemp+1))-2*Delta.at(itemp)*Fk(2,L.at(0),L.at(0)+1,Spin.at(itemp),Spin.at(itemp+1))+Delta.at(itemp)*Delta.at(itemp)*Fk(2,L.at(0)+1,L.at(0)+1,Spin.at(itemp),Spin.at(itemp+1)))/(1+Delta.at(itemp)*Delta.at(itemp));
            Ak = (Fk(2,L.at(1),L.at(1),Spin.at(itemp+2),Spin.at(itemp+1))-2*Delta.at(itemp+1)*Fk(2,L.at(1),L.at(1)+1,Spin.at(itemp+2),Spin.at(itemp+1))+Delta.at(itemp+1)*Delta.at(itemp+1)*Fk(2,L.at(1)+1,L.at(1)+1,Spin.at(itemp+2),Spin.at(itemp+1)))/(1+Delta.at(itemp+1)*Delta.at(itemp+1));
            W -> SetParameter(0,Bk*Ak);
            W -> SetParameter(1,Fk(4,L.at(0),L.at(0),Spin.at(itemp),Spin.at(itemp+1))*Fk(4,L.at(1),L.at(1),Spin.at(itemp+2),Spin.at(itemp+1)));
            Angle.push_back(W -> GetRandom());
          // cout << "angle = " << Angle.at(0) << endl;
          // cout << "coefficients are " << W->GetParameter(0) << "  " << W->GetParameter(1) << endl;
            L.clear();
          }
          else if ((IC.at(itemp)==3))
          {
            Angle.push_back(TMath::Pi());
            Angle.push_back(-1e3);
          }
          else if ((IC.at(itemp+1)==3))
          {
            Angle.push_back(-1e3);
            // Angle.push_back(TMath::Pi());
          }
          else{
            Angle.push_back(-1e3);
          }
        }
        // if(multiplicity>2){cout << "angles = "; for(int i=0; i<Angle.size(); ++i){cout << Angle.at(i) << "  ";} cout << endl;}
      }
      cascades->Fill();
      temp.clear();
      if (Eg.size()>0){
        outfile << Eg.size() << " ";
        for(int itemp=0; itemp<Eg.size(); ++itemp){
          outfile << Eg.at(itemp) << " ";}
        for(int itemp=0; itemp<Angle.size(); ++itemp){
          outfile << Angle.at(itemp) << " ";}
        outfile << endl;
      }
      // clean up the present cascade to get ready for the next one
      Eg.clear(); Ee.clear(); Epair.clear(); Spin.clear(); Delta.clear(); IC.clear(); Angle.clear();
    }
  }
  cout << "done reading " << argv[1] << "/"<< argv[2] << endl;
  file.close();
  outfile.close();

  cascades->Write("cascades",TObject::kOverwrite);
 //___________
 //here some basic plots
  TCanvas *canv = new TCanvas();
  canv->Divide(3,2);

  canv->cd(3);
  // gPad->SetLogy();
  cascades->Draw("angle>>angle_distribution(360,-TMath::Pi(),+TMath::Pi())");

  canv->cd(4);
  gPad->SetLogy();
  cascades->Draw("gammas>>h_gammas_temp(1000,0,10)");
  h_gammas = (TH1F*)gDirectory->Get("h_gammas_temp");
  h_gammas->SetDirectory(0);
  h_gammas->Write("h_gammas",TObject::kOverwrite);

  canv->cd(6);
  gPad->SetLogy();
  cascades->Draw("pairs>>h_pairs(1000,0,10)");

  canv->cd(1);
  gPad->SetLogy();
  cascades->Draw("Sum$(gammas)>>h_g_sums_temp(1000,0,10)");
  h_g_sums =(TH1F*)gDirectory->Get("h_g_sums_temp");
  h_g_sums->SetDirectory(0);
  h_g_sums->Write("h_g_sums",TObject::kOverwrite);

  canv->cd(2);
  cascades->Draw("Length$(gammas)>>h_g_M_temp(10,-0.5,9.5)");
  TH1F* h_g_M = (TH1F*)gDirectory->Get("h_g_M_temp");
  h_g_M->Write("h_g_M",TObject::kOverwrite);

  canv->cd(5);
  gPad->SetLogy();
  cascades->Draw("electrons>>h_elns(1000,0,10)");

  canv->SaveAs(Form("%s/%s.pdf",argv[1],argv[2]));

  root_file->Close();
  cout << "file " << root_file->GetName() << " was written" << endl;
  
  closedir(d);
  return 0;
}
