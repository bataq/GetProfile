// C++ file
// Auther KAWABATA Tomoki
// Date 2017-03-01

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TRandom.h>
#include <TMath.h>

#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <boost/program_options.hpp>
using namespace boost::program_options;

//#define __CINT__

void WriteLogFile(int enem,int eneM,int w,int type);

#ifdef __CINT__
int GetProfile(const TString rootFile){
#else
  int GetProfile(const TString rootFile,Int_t energy_m,Int_t energy_M,Int_t width,Int_t ty);
  int GetProfilePave(const TString rootFile,const TString onepixProfile,Int_t energy_m,Int_t energy_M,Int_t width,Int_t ty);

int main(int argc, char* argv[]){

  // Option Analysis //
  TString dataFile,onepixFile;
  Int_t enem=0,eneM=0,w=0,type=0;
  options_description opt("Options");
  opt.add_options()
    ("file,f",value<std::string>(),"Event List ROOT File")
    ("onepix,o",value<std::string>(),"One Pixel Period File")
    //("subpix,s","Use Subpix Event List")
    //("onepix,o","Use Onepix Event List")
    ("type,t",value<int>(),"Event Type (1:Single & Side Double, 2:Side Double)")
    ("enem,m",value<int>()->default_value(50),"Couning Energy Minimum (PH)")
    ("eneM,M",value<int>()->default_value(100),"Couning Energy Maximum (PH)")
    ("width,w",value<int>()->default_value(40),"Width of Profile Range (pix)")
    ("help,h","HELP");
  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);
  if( vm.count("help") || !vm.count("file")){
    std::cout << "Usage  :./GetProfile -f <Event List Root File> -s (OR -o) -m <Energy Min(PH)> -M <Energy Max(PH)> -w <Width of Profile Range>" << std::endl;
    std::cout << opt << std::endl;
    return 0;
  }else{
    dataFile = vm["file"].as<std::string>();
    onepixFile  = vm["onepix"].as<std::string>();
    type  = vm["type"].as<int>();
    enem  = vm["enem"].as<int>();
    eneM  = vm["eneM"].as<int>();
    w  = vm["width"].as<int>();
    std::cout << "dataFile is "<<dataFile<<std::endl;
    //    if(vm.count("subpix")){flag=1}
  }
  
  if( !vm.count("onepix") ){
    GetProfile(dataFile,enem,eneM,w,type);
  }else{
    GetProfilePave(dataFile,onepixFile,enem,eneM,w,type);
  }    
  return 0;
}

 int GetProfile(const TString rootFile,Int_t energy_m,Int_t energy_M,Int_t width,Int_t ty)
 {
#endif

  if(!rootFile){
    std::cerr<<"Usage: GetProfile(root file name)"<<std::endl;
    return 0;
  }

  std::cout<<"rootFile name is "<<rootFile<<std::endl;
  TFile *f = new TFile(rootFile,"read");
  if(!f){
    std::cerr<<"Error: no such file." << std::endl;
    return 0;
  }
  
  TTree *str=(TTree*)f->Get("tree_subpix");
  Double_t        sra;
  Double_t        sca;
  Int_t        stype;
  Float_t         ph_merge;
  TBranch        *b_sra;
  TBranch        *b_sca;
  TBranch        *b_stype;
  TBranch        *b_ph_merge;
  str->SetBranchAddress("sca",&sca,&b_sca);
  str->SetBranchAddress("sra",&sra,&b_sra);
  str->SetBranchAddress("stype",&stype,&b_stype);
  str->SetBranchAddress("ph_merge",&ph_merge,&b_ph_merge);

  Int_t allEntries=str->GetEntries();
    
  const Int_t ANGLE = 0;  
  Double_t rad = ANGLE * TMath::Pi()/ 180;
  Int_t wh;
  Int_t xx_m = 71*cos(rad)-71*sin(rad)+143*sin(rad)-60;
  Int_t xx_M = 71*cos(rad)-71*sin(rad)+143*sin(rad)+60;
  Int_t yy_m;
  Int_t yy_M;
Int_t count[2048] = {};
TString filename;
if(ty==1){filename="profileSingleSideDouble.txt";}
if(ty==2){filename="profileSideDouble.txt";}

  if(width%2==0){
    wh=width/2;
    yy_m = 71*sin(rad)+71*cos(rad)-wh;
    yy_M = 71*sin(rad)+71*cos(rad)+wh;
  }else{
    wh=(width+1)/2;
    yy_m = 71*sin(rad)+71*cos(rad)-wh;
    yy_M = 71*sin(rad)+71*cos(rad)+wh-1;
  }


  for(int i=0; i<allEntries;i++){
    str->GetEntry(i);
    int m=0;
    if(ty==1){
      if(stype==10||stype==11||stype==22){
	if (sca >= xx_m-0.5 && sca < xx_M+0.5){
	  if(sra >= yy_m && sra <= yy_M){
	    if(ph_merge >=energy_m && ph_merge<=energy_M){
	      m=int(10*sca);// xxを切り捨て
	      m=m-10*(xx_m-0.5);// count[]のindex作成
	      count[m]++;
	    }
	  }
	}
      }
    }if(ty==2){
      if(stype==22){
	if (sca >= xx_m-0.5 && sca < xx_M+0.5){
	  if(sra >= yy_m && sra <= yy_M){
	    if(ph_merge >=energy_m && ph_merge<=energy_M){
	      m=int(10*sca);// xxを切り捨て
	      m=m-10*(xx_m-0.5);// count[]のindex作成
	      count[m]++;
	    }
	  }
	}
      }
    }
  }
  std::ofstream fout(filename);
  std::cout<<"make "<<filename<<std::endl;
  fout << "sca"<< std::setw(10) <<"intensity"<< std::endl;

  for (int m=0; m<10*(xx_M-xx_m+1);m++) {	
  Float_t x=Float_t(m+10*(xx_m-0.5))/10;
   fout <<std::fixed<< std::setprecision(1) << x << std::setw(10) << count[m] << std::endl;
  }
  fout.close();
  WriteLogFile(energy_m,energy_M,width,ty);
  delete str;
  return 0;
}

 int GetProfilePave(const TString rootFile,const TString onepixProfile,Int_t energy_m,Int_t energy_M,Int_t width,Int_t ty)
 {

  if(!rootFile){
    std::cerr<<"Usage: GetProfile(root file name)"<<std::endl;
    return 0;
  }

  std::cout<<"rootFile name is "<<rootFile<<std::endl;
  TFile *f = new TFile(rootFile,"read");
  if(!f){
    std::cerr<<"Error: no such file." << std::endl;
    return 0;
  }
  
  TTree *str=(TTree*)f->Get("tree_subpix");
  Double_t        sra;
  Double_t        sca;
  Int_t        stype;
  Float_t         ph_merge;
  TBranch        *b_sra;
  TBranch        *b_sca;
  TBranch        *b_stype;
  TBranch        *b_ph_merge;
  str->SetBranchAddress("sca",&sca,&b_sca);
  str->SetBranchAddress("sra",&sra,&b_sra);
  str->SetBranchAddress("stype",&stype,&b_stype);
  str->SetBranchAddress("ph_merge",&ph_merge,&b_ph_merge);

  Int_t allEntries=str->GetEntries();
    
  const Int_t ANGLE = 0;  
  Double_t rad = ANGLE * TMath::Pi()/ 180;
  Int_t wh;
  Int_t xx_m = 71*cos(rad)-71*sin(rad)+143*sin(rad)-60;
  Int_t xx_M = 71*cos(rad)-71*sin(rad)+143*sin(rad)+60;
  Int_t yy_m;
  Int_t yy_M;
  double count[2048] = {};
  TString filename;
  if(ty==1){filename="profileSingleSideDoublePave.txt";}
  if(ty==2){filename="profileSideDoublePave.txt";}

  std::ifstream onepixFile(onepixProfile);
  double a;
  double onepix[10];
  int i=0;
  while(onepixFile>>a){
    onepix[i]=a;
    i++;
  }

  if(width%2==0){
    wh=width/2;
    yy_m = 71*sin(rad)+71*cos(rad)-wh;
    yy_M = 71*sin(rad)+71*cos(rad)+wh;
  }else{
    wh=(width+1)/2;
    yy_m = 71*sin(rad)+71*cos(rad)-wh;
    yy_M = 71*sin(rad)+71*cos(rad)+wh-1;
  }


  for(int i=0; i<allEntries;i++){
    str->GetEntry(i);
    int m=0;
    if(ty==1){
      if(stype==10||stype==11||stype==22){
	if (sca >= xx_m-0.5 && sca < xx_M+0.5){
	  if(sra >= yy_m && sra <= yy_M){
	    if(ph_merge >=energy_m && ph_merge<=energy_M){
	      m=int(10*sca);// xxを切り捨て
	      m=m-10*(xx_m-0.5);// count[]のindex作成
	      count[m]++;
	    }
	  }
	}
      }
    }else if(ty==2){
      if(stype==22){
	if (sca >= xx_m-0.5 && sca < xx_M+0.5){
	  if(sra >= yy_m && sra <= yy_M){
	    if(ph_merge >=energy_m && ph_merge<=energy_M){
	      m=int(10*sca);// xxを切り捨て
	      m=m-10*(xx_m-0.5);// count[]のindex作成
	      count[m]++;
	    }
	  }
	}
      }
    }
  }
  std::ofstream fout(filename);
  std::cout<<"make "<<filename<<std::endl;
  fout << "sca"<< std::setw(12) <<"intensity/pane"<< std::setw(10) <<"err" << std::endl;

  for (int m=0; m<10*(xx_M-xx_m+1);m++) {	
   Float_t x=Float_t(m+10*(xx_m-0.5))/10; 
   i=m%10;
   fout <<std::fixed<< std::setprecision(1) << x << std::setw(12);
   fout<<std::fixed<< std::setprecision(7)  << count[m]/onepix[i] << std::setw(12)<<  sqrt(count[m])/onepix[i] << std::endl;
  }

  fout.close();
  WriteLogFile(energy_m,energy_M,width,ty);
  delete str;
  return 0;
}

 void WriteLogFile(int enem,int eneM,int w,int type){
   std::ofstream f("log_GetProfile.txt");
   time_t t =time(NULL);
   f<<"Date :"<<ctime(&t)<<std::endl;
   f<<"Energy Minimum : "<<enem<<std::endl;
   f<<"Energy Maxmum : "<<eneM<<std::endl;
   f<<"width : "<<w<<"(RA)"<<std::endl;
   f<<"type : "<<type<<std::endl;
   f.close();
   return;
 }

 
