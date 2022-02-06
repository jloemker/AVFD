#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <map>
#include "TFile.h"
#include "TTree.h"
#include "Event.h"
#include "Particle.h"
const int nPid = 29;//changed from 29
//https://github.com/chunshen1987/iEBE/blob/54435e159f8cc0d05a0b32b78902e2569fcf4016/EBE-Node/EbeCollector/EbeCollector.py
const int UrQMDpid[nPid] = {2101, -1899, 101 , 1106, -894, -1106, 894, 1001, -999, -1001, 999, 2040, -1960, 40, -2040, 1960, -40, 1049, -951, -1049, 951, 27, -27, 55, -55, 109, 102, 107, 100};
                  
const int PDGpid[nPid] = {211, -211, 111, 321, 311, -321, -311, 2212, 2112, -2212, -2112, 3222, 3112, 3212, -3222, -3112, -3212, 3322, 3312, -3322, -3312, 3122, -3122, 3334, -3334, 333, 221, 331, 22};                    
/*
class Particle { public:
 public:
  Particle(int pid, int charge, float energy, float rapidity, float pt,
           float phi)
      : pid_(pid),
        charge_(charge),
        energy_(energy),
        rapidity_(rapidity),
        pt_(pt),
        phi_(phi){};
  Particle() = default;

  int getPid() { return pid_; }
  int getCharge() { return charge_; }
  double getEnergy() { return energy_; }
  double getRapidity() { return rapidity_; }
  double getPt() {return pt_; }
  double getPhi() {return phi_; }

 private:
  int pid_{0};
  int charge_{0};
  float energy_{-999.};
  float rapidity_{-999.};
  float pt_{-999.};
  float phi_{-999};
};

class Event {
 private:
  std::vector<Particle> particles_;
  int nTrack_;

 public:
  Event(int nTrack, std::vector<Particle> particles) : nTrack_(nTrack), particles_(particles){};
  Event() = default;
  Particle getParticle(int ith) { return particles_.at(ith); };
  int getnTrack() { return nTrack_; };
};

*/
void convert_tree_splitFiles(Int_t centID, Int_t dirID) 
{
//for(Int_t centID = 0; centID < 7; centID++){
string dir;
if(dirID == 0){
dir = "Baseline";
} 
if(dirID ==1){
dir = "a-0.1";
}
 std::cout <<"Centrality ID: " << centID << " and directory: " << dir << std::endl;
 //Int_t centID = {0,1,2,3,4,5,6,7}; //0: 0-5%, 1: 5-10%, 2: 10-20%, 3: 20-30%, 4: 30-40%, 5: 40-50%, 6: 50-60%, 7: 60-70%

  int nJob = 999;//was 2000 but panos has only 999 jobs in baseline tir
  int nFile = 50;//was 201
  int nFinalFiles =20;//was 10
  // map to convert UrQMD pid to PDG pid
  std::map<int, int> pid_conversion_map;
  for (int ithPid = 0; ithPid < nPid; ithPid++) {  
	pid_conversion_map[UrQMDpid[ithPid]] = PDGpid[ithPid];
  }
  
  int jobNum = 0;
  
  int nPosChargeTotal = 0;
  int nNegChargeTotal = 0;
  int nNeutralChargeTotal = 0;
  int nPosPion{0};
  int nNegPion{0};
  int nPosKaon{0};
  int nNegKaon{0};
  int nPosProton{0};
  int nNegProton{0};
  string directory;
  
  int val1 = (centID-1)*10;
  int val2 = (centID)*10;
  if (centID == 0 || centID == 1) {
  val1 = (centID)*5;
  val2 = (centID+1)*5;
  }

  for (int nthOutput = 0; nthOutput < nFinalFiles; nthOutput++) {
	  TFile file(Form("tree_5.02TeV_Cent%d_%d_%d.root", val1, val2, nthOutput),"RECREATE");//submit puts me into dirID/...
	  TTree tree("events", "event");
	  Event ev;
	  tree.Branch("event", &ev);
	  for (int ithJob = jobNum; ithJob <= nJob; ithJob++){ //jobNum+(nJob/nFinalFiles)+(nJob%nFinalFiles); ithJob++) {
		for (int ithFile = 1; ithFile <= nFile; ithFile++) {
		  std::ifstream file_dat;
		if(dirID == 0){//only 30_40 exists!
		  directory = Form("/dcache/alice/jlomker/sim/NoBField/5.02TeV/Centrality%d_%d/job-%d/Result/event-1/particle_distribution_final/%d.dat",val1,val2,ithJob,ithFile);
		  //directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.02TeV/Centrality%d-%d/Baseline/job-%d/particle_distribution_final/%d.dat",val1, val2, ithJob,ithFile);
                  }else{
		  directory = Form("/dcache/alice/panosch/alice/sim/2020/AVFD/5.02TeV/Centrality%d-%d/a-0.1/job-%d/particle_distribution_final/%d.dat",val1, val2, ithJob,ithFile);
  		  }	
		  file_dat.open(directory.c_str());
		  
		  if (!file_dat.is_open()) {
			  std::cout<<"Cannot find "<<directory<<std::endl;
			  continue;
		  } else {
			  if (ithJob%10==0 && ithFile == 1) std::cout<<"Processing job #"<<ithJob<<std::endl;
		  }

		  std::string line;
		  std::vector<Particle> p;

		  int pid{0};
		  int charge{0};
		  float energy{-999.};
		  float rapidity{-999.};
		  float pt{-999.};
		  float phi{-999};
		  int nTrack{0};
		  

		  
		  //while (std::getline(file_dat, line)) {
		  while (file_dat >> pid >> charge >> energy >> rapidity >> pt >> phi) {
		    if (charge == -1) {
		      nNegChargeTotal++;
            } else if (charge == 0) {
			  nNeutralChargeTotal++;
		    } else if (charge == 1) {
			  nPosChargeTotal++;
            }
			nTrack++;
			
			pid = pid_conversion_map[pid];
			
			if (pid == 211) {//PosPion
              nPosPion++;
            } else if (pid == -211) {
              nNegPion++;
            } else if (pid == 321) {//PosKaon
              nPosKaon++;
            } else if (pid == -321) {
              nNegKaon++;
            } else if (pid == 2212) {
              nPosProton++;
            } else if (pid == -2212) {
              nNegProton++;
            }	
			// std::cout << line << std::endl;
			//std::cout << "pid ==== "<< pid <<", charge = "<<charge<<", energy = "<<energy<< std::endl;
			auto part = Particle(pid, charge, energy, rapidity, pt, phi);
			p.push_back(part);
		  }

		  ev = Event(nTrack, p);

		  tree.Fill();
		  file_dat.close();
		}
		
	  }
	  jobNum = jobNum+(nJob/nFinalFiles)+(nJob%nFinalFiles)+1;
	  tree.Write();
  }
  cout<<"===> nNegChargeTotal === "<<nNegChargeTotal<<endl;
  cout<<"===> nNeutralChargeTotal === "<<nNeutralChargeTotal<<endl;
  cout<<"===> nPosChargeTotal === "<<nPosChargeTotal<<endl;
  cout<<"===> nPosPion ======= "<<nPosPion<<endl;
  cout<<"===> nNegPion ======= "<<nNegPion<<endl;
  cout<<"===> nPosKaon ======= "<<nPosKaon<<endl;
  cout<<"===> nNegKaon ======= "<<nNegKaon<<endl;
  cout<<"===> nPosProton ======= "<<nPosProton<<endl;
  cout<<"===> nNegProton ======= "<<nNegProton<<endl;
//}
}
