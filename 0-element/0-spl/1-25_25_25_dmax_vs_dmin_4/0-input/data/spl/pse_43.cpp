//
//  main.cpp
//  pse
//
//  Created by Mahdi Taiebat on 2016-06-21.
//  Copyright © 2016 UBC. All rights reserved.
//

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <cstdlib>
#include <set>
#include <string>
#include <stdio.h>
#include <random>
#include <math.h>
#include <algorithm>
#include <iterator>
#include <stdlib.h>
#include <stack>
#include <ctime>

#define MAXIT 100
#define EPS 3.0e-7
#define FPMIN 1.0e-30

std::stack<clock_t> tictoc_stack;
using namespace std;

// Function declarations
void header();
void readFlow();
void readParameters();
void setParticleSize();
void setParticlePosition();
void writeSample();
void writeConfig();
double betacf(double a, double b, double x);
double gammln(double xx);
double betai(double a, double b, double x);
void initCompaction1();
void initCompaction2();
void initCompaction3();
void VerletList();
void CompactionLoop1();
void CompactionLoop2();
int CompactionLoop3();
void writevtk(const long int ivtk);
void tic();
void toc();
void CalcPackingFraction();
void CalcCoordinationNumber();
int sit();
void WriteScreen();
void WriteFile();
vector<double> linspace();
//void writeMCprogress(const long int& iStop, const long int& iMove, const long int& iRejectMove, const long int& iAcceptMove, const double& RejectionRate, const double& phi, const double& z, int flag);

// Type declaration for global parameters
// Parameters for DEM based compaction
size_t iFlow,iWriteScreen,iWriteFile,nWriteScreen,nWriteFile,iGrading,iWall,iConfig;
size_t Npar;
double Dmin,SizeRatio,Dtopwall;
double a,b,amin,amax,bmin,bmax;
size_t Nx, Ny;
double Gap;
double Gap_imperfection;
vector<double> D, rx, ry, rz;
string IutputFile_GradingCurve, OutputFile_ClassGradingCurve, OutputFile_ClassSizeDistribution, OutputFile_SampleBoxSize,OutputFile_sample,OutputFile_sample_radius,OutputFile_config,IutputFile_Config1,IutputFile_Config2;
size_t NClass, Nmin_pc;
// Additional parameters for Monte Carlo Compaction
size_t iCompaction,iTopWallControl;
string IutputFile_Sample;
double Kd, KMao;
double DverFactor,ROverlap;
long int Nver,NMove,NStop,Nvtk;
vector<size_t> vlist,ipoint,vlistWall6;
vector<int> vlistVx,vlistVy;
double BoxDimLow_x,BoxDimLow_y,BoxDimLow_z,BoxDimHigh_x,BoxDimHigh_y,BoxDimHigh_z;
double CVDimHigh_z;
double Ld,LMao,Dver,BetaRef,rand_Ld,Ldmin,Ldmax;
double RepulsionFactor,GravityFactor,ForceTopWall;
double kcv;
double phi,phi_last,delta_phi,sum_delta_phi;
double z;
double DilationFactor = 0.0;
string OutputFile_MCprogress;
double Lx,Ly,Lz;
double DperFactor,Dper;
double dmin,dmax;
double MeanOverlapRatio;
long int NOverlap;
double TopWallOverlapRatio;
long int NTopWallOverlap;
double UFactorOld,UFactorNew;
ofstream MCprogress;
long int iRejectMove;
long int iAccept = 0; // cumulative number of success
double RejectionRateOld = 0;
double RejectionRate = 0;
long int iAcceptMove;
long int iStop,iMove;
int Nangle,Nrelax,dphiCheckN;
double dphiCheckTol;
size_t iBetaFn;
string OutputFile_BetaFnSensitivity;
double B_beta;

int myrandom (int i) { return std::rand()%i;}

// Write a header line
void header(){
  cout << "Generate file sample.spl" << endl;
}

// Read flow parameters
void readFlow() {
  string line;
  ifstream myFile;  //buffer
  myFile.open("Flow.dat");
  
  // Error check
  if(myFile.fail()){
    cerr << "* File can not be found or opened.";
    exit(1);
  }
  
  istringstream iss;
  while(myFile.good()) {
    getline(myFile, line);  // header of the input file
    
    getline(myFile, line);
    iss.str(line);
    iss >> iFlow;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iWriteScreen >> nWriteScreen;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iWriteFile >> nWriteFile;
    
    if (iFlow == 0){
      sit();
    }
  }
  myFile.close();
  return ;
}

// Read input parameters
void readParameters() {
  
  string line;
  ifstream myFile;  //buffer
  myFile.open("InputParameters.dat");
  
  // Error check
  if(myFile.fail()){
    cerr << "* File can not be found or opened.";
    exit(1);
  }
  
  istringstream iss;
  
  while(myFile.good()) {
    getline(myFile, line);  // header of the input file
    getline(myFile, line);  // header of the input file
    getline(myFile, line);  // header of the input file
    getline(myFile, line);  // header of the input file
    
    getline(myFile, line);
    iss.str(line);
    iss >> iGrading;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Npar;
    
    getline(myFile, line);
    iss.str(line);
    iss >> NClass >> Nmin_pc;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Dmin >> SizeRatio;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Dtopwall;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iBetaFn;
    
    getline(myFile, line);
    iss.str(line);
    iss >> amin >> amax >> bmin >> bmax;

    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_BetaFnSensitivity;
    
    getline(myFile, line);
    iss.str(line);
    iss >> a >> b;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Nx >> Ny;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Gap;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Gap_imperfection;
    
    getline(myFile, line);
    iss.str(line);
    iss >> IutputFile_GradingCurve;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_ClassGradingCurve;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_ClassSizeDistribution;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_SampleBoxSize;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_sample;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_sample_radius;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_config;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iWall;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iConfig;
    
    getline(myFile, line);
    iss.str(line);
    iss >> IutputFile_Config1;
    
    getline(myFile, line);
    iss.str(line);
    iss >> IutputFile_Config2;

    getline(myFile, line);  // Parameters for Monte Carlo Compaction
    
    getline(myFile, line);
    iss.str(line);
    iss >> iCompaction;
    
    getline(myFile, line);
    iss.str(line);
    iss >> iTopWallControl;
    
    getline(myFile, line);
    iss.str(line);
    iss >> IutputFile_Sample;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Kd >> KMao;
    
    getline(myFile, line);
    iss.str(line);
    iss >> DverFactor >> Nver;
    
    getline(myFile, line);
    iss.str(line);
    iss >> ROverlap;
    
    getline(myFile, line);
    iss.str(line);
    iss >> NMove >> NStop;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Nvtk;
    
    getline(myFile, line);
    iss.str(line);
    iss >> kcv;
    
    getline(myFile, line);
    iss.str(line);
    iss >> OutputFile_MCprogress;
    
    getline(myFile, line);
    iss.str(line);
    iss >> DperFactor;
    
    getline(myFile, line);
    iss.str(line);
    iss >> Nangle >> Nrelax;
    
    getline(myFile, line);
    iss.str(line);
    iss >> dphiCheckN >> dphiCheckTol;
  }
  
  myFile.close();
  return ;
}

// Define the array of particle radii
void setParticleSize() {
  
  double Dmax=SizeRatio*Dmin;
  double x_value, h_value, x1=0., x2=0., h1=0., h2=0.;
  vector <double> x, h, xClass, hClass, nTemp1Class, nTemp2Class, alphaClass, nClass;
  //istringstream iss;
  
  if (iConfig == 1) {
    string line;
    istringstream iss;
    ifstream myFile;
    myFile.open (IutputFile_Config1);
    getline(myFile, line);  // reading Dmax in the first line
    iss.str(line);
    string temp1;
    double temp2;
    iss >> temp1 >> temp2;
    Dmax = temp2;
    Dmin=Dmax/SizeRatio;
    myFile.close();
  }

  if (iConfig == 2) {
    string particle;
    int grp;
    double Rval, rxval, ryval, rzval;
    string line;
    istringstream iss;
    ifstream myFile;
    myFile.open (IutputFile_Config2);
    getline(myFile, line);  // reading "sample{" in the first line
    Npar = 0;
    while(myFile.good()) {
      getline(myFile, line);
      iss.str(line);
      iss >> particle;
      if (particle == "sphere"){
        iss >> grp >> Rval >> rxval >> ryval >> rzval;
        if (grp==0) {
          Npar++;
        }
      }
    }
    myFile.close();
  }
  
  if (iGrading == 1) {
    for (size_t ipar = 0; ipar < Npar; ipar++){
      double s= (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
      double Dval = Dmin+s*(Dmax-Dmin);
      D.push_back(Dval);
    }
  }
  else {    // for iGrading of 0 and 2
    //======================================================
    if (iGrading == 0) {   // read x and h from GradingCurveFile.dat
      string line;
      ifstream myFile;  //buffer
      myFile.open(IutputFile_GradingCurve);
      
      // Error check
      if(myFile.fail()){
        cerr << "* GradingCurveFile can not be found or opened.";
        exit(1);
      }
      
      getline(myFile, line); // header line of GradingCurveFile
      cout << line << '\n';
      
      while(getline(myFile, line)) {
        stringstream iss(line);
        //iss.str(line);
        iss >> x_value >> h_value;
        x.push_back(x_value);
        h.push_back(h_value);
        cout << x_value << "  " << h_value << '\n';
      }
      cout << "--------Done with creating the x and h vectores!" << '\n';
      myFile.close();
    }
    else {
      cout << "The h(x) prescribed by a function!" << '\n';
    }
    
    // Create xClass
//    cout << "xClass values:" << '\n';
    double xClassHalfStep = 0.5 * (Dmax - Dmin)/(NClass-1);
    double DminCalc = Dmin - xClassHalfStep;
    double DmaxCalc = Dmax + xClassHalfStep;
    if (DminCalc == 0.){
      DminCalc = 1e-10; // This is to avoid division by zero when Dmax/Dmin = 2*NClass - 1 which results in DminCalc=0
    }
    double SizeRatioCalc = DmaxCalc/DminCalc;
    for (int i = 0; i < NClass; i++) {
      double xClass_val = Dmin + (Dmax-Dmin)/(NClass-1)*i;
      xClass.push_back(xClass_val);
//      cout << xClass_val << '\n';
    }
//    cout << "--------Done with creating the xClass!" << '\n';
    
    // calculate hClass and NClass
    for (unsigned i=0; i < xClass.size(); i++){     // march over 0 to size(xClass)-1
      
      if (iGrading == 0) {   // Find x1 & x2 and h1 & h2 and then hClass for a given xClass
        cout << "--------Generating the x1 and x2 values for xClass_val = " << xClass.at(i) << '\n';
        for (unsigned j=0; j < x.size(); j++){ // march over 0 to size(x)-1
          if (x.at(j) > xClass.at(i)) {
            x1  = x.at(j-1);
            x2  = x.at(j);
            h1  = h.at(j-1);
            h2  = h.at(j);
//            cout << "x1= " << x1 << ",   x2= " <<  x2 << '\n';
//            cout << "h1= " << h1 << ",   h2= " <<  h2 << '\n';
            break;
          }
        }
        double hClass_val = h1 + (xClass.at(i)-x1)/(x2-x1)*(h2-h1);
        hClass.push_back(hClass_val);
      } else if  (iGrading == 2) {  // Find hClass for a given xClass
        double hClass_val = betai(a,b,((xClass.at(i)+ xClassHalfStep) -DminCalc)/(DminCalc*(SizeRatioCalc-1.))) * 100.;
        hClass.push_back(hClass_val);
      }
      
//      cout << "xClass= " << xClass.at(i) << ",   hClass= " <<  hClass.at(i) << '\n';
      
      if (i==0){
        double nTempClass_val = (hClass.at(i) - 0) / (M_PI/6.*pow(xClass.at(i), 3.));
        nTemp1Class.push_back(nTempClass_val);
      }
      else{
        double nClass_val = (hClass.at(i) - hClass.at(i-1)) / (M_PI/6.*pow(xClass.at(i), 3.));
        nTemp1Class.push_back(nClass_val);
      }
    }
    
//    cout << "--------Here are the temporary number of particles in various xClass categories:" << '\n';
    double Nmin_pc_Coeff = Nmin_pc / (*min_element(nTemp1Class.begin(), nTemp1Class.end()));
    for (int i=0; i < nTemp1Class.size(); i++){
      double nClass_val = floor(Nmin_pc_Coeff * nTemp1Class.at(i));
      nTemp2Class.push_back(nClass_val);
//      cout << "xClass= " << xClass.at(i) << ",  nTemp2Class= " << nTemp2Class.at(i) << '\n';
    }
    
    double sum_nTemp2Class=0;
    for (int i=0; i < nTemp2Class.size(); i++){
      sum_nTemp2Class += nTemp2Class.at(i);
    }
//    cout << "--------Here are the final number of particles in various xClass categories:" << '\n';
    
    for (int i=0; i < nTemp2Class.size(); i++){
      if (Npar > sum_nTemp2Class ) {
        double nClass_val = floor(Npar/sum_nTemp2Class * nTemp2Class.at(i));
        nClass.push_back(nClass_val);
      } else {
        double nClass_val = nTemp2Class.at(i);
        nClass.push_back(nClass_val);
      }
//      cout << "xClass= " << xClass.at(i) << ",  nClass= " << nClass.at(i) << '\n';
    }
    
    
    double Npar_New=0;
    for (int i=0; i < nClass.size(); i++){
      Npar_New += nClass.at(i);
    }

    if (Npar_New < Npar) {
      nClass.at(nClass.size()-1) = nClass.at(nClass.size()-1) + (Npar - Npar_New);
    }
    
    Npar_New=0;
    for (int i=0; i < nClass.size(); i++){
      Npar_New += nClass.at(i);
    }

    Npar = Npar_New;
    
    cout << "Npar_New= " << Npar_New << '\n';
//    cout << "xClass.size= " << xClass.size() << '\n';
//    cout << "nClass.size= " << nClass.size() << '\n';
    
    if (iBetaFn == 0){
      for (int i = 0; i < nClass.size(); i++){
        for (int j = 0; j < nClass.at(i); j++){
          double Dval = xClass.at(i);
          D.push_back(Dval);
        }
      }
    }
    std::random_shuffle ( D.begin(), D.end() , myrandom);
    
    ofstream newFile;
    //cout << OutputFile_ClassGradingCurve;
    newFile.open (OutputFile_ClassGradingCurve);
    newFile << "# OutputFile_ClassGradingCurve: sieve openning (m), percentage finer (%)" << endl;
    for (size_t ipar = 0; ipar < nClass.size(); ipar++){
      newFile << xClass.at(ipar) << "  " << hClass.at(ipar) << endl;
    }
    newFile.close();
    
    //cout << OutputFile_ClassSizeDistribution;
    newFile.open (OutputFile_ClassSizeDistribution);
    newFile << "# OutputFile_ClassSizeDistribution: particle diameter (m), percentage (%)" << endl;
    for (size_t ipar = 0; ipar < nClass.size(); ipar++){
      newFile << xClass.at(ipar) << "  " << nClass.at(ipar) << endl;
    }
    newFile.close();
    //======================================================
  }
  
}

// Define the arrays of particle coordinates
void setParticlePosition() {
  
  size_t Nz = ceil(Npar/(Nx*Ny*1.));
  double LatticeStep=SizeRatio*Dmin*(1.+Gap);
  if (iConfig == 0){
    for (size_t iz = 0; iz < Nz; iz++){
      for (size_t iy = 0; iy < Ny; iy++){
        for (size_t ix = 0; ix < Nx; ix++){
          size_t ipar = ix + Nx * iy + Nx * Ny * iz;
          if (ipar < Npar){    // this can be improved!
            double randx= (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            double randy= (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            double randz= (static_cast <float> (rand()) / static_cast <float> (RAND_MAX));
            
            double imperfection_x = (2.*randx - 1.)*Gap*Gap_imperfection*SizeRatio*Dmin;
            double imperfection_y = (2.*randy - 1.)*Gap*Gap_imperfection*SizeRatio*Dmin;
            double imperfection_z = (2.*randz - 1.)*Gap*Gap_imperfection*SizeRatio*Dmin;
            
            rx.push_back(ix*LatticeStep+imperfection_x);
            ry.push_back(iy*LatticeStep+imperfection_y);
            rz.push_back(iz*LatticeStep+imperfection_z);
          }
        }
      }
    }
  }
  if (iConfig == 1) {
    double rxval, ryval, rzval;
    string line;
    istringstream iss;
    ifstream myFile;
    myFile.open (IutputFile_Config1);
    getline(myFile, line);  // reading dmax in the first line
    getline(myFile, line);  // header of the input file
    for(size_t i = 0 ; i < Npar ; i++){
      getline(myFile, line);
      //find a command to report error message if there was nothing in "line"
      iss.str(line);
      iss >> rxval >> ryval >> rzval;
      rx.push_back(rxval);
      ry.push_back(ryval);
      rz.push_back(rzval);
    }
    myFile.close();
  }
  if (iConfig == 2) {
    string particle;
    int grp;
    double Rval, rxval, ryval, rzval;
    string line;
    istringstream iss;
    ifstream myFile;
    myFile.open (IutputFile_Config2);
    getline(myFile, line);  // reading "sample{" in the first line
    Npar = 0;
    while(myFile.good()) {
      getline(myFile, line);
      iss.str(line);
      iss >> particle;
      if (particle == "sphere"){
        iss >> grp >> Rval >> rxval >> ryval >> rzval;
        if (grp==0) {
          Npar++;
          rx.push_back(rxval);
          ry.push_back(ryval);
          rz.push_back(rzval);
        }
      }
    }
    myFile.close();
  }
  
  //  BoxDimLow_x     = -1/2.*LatticeStep;
  //  BoxDimLow_y     = -1/2.*LatticeStep;
  //  BoxDimLow_z  = -1/2.*LatticeStep;
  //  BoxDimHigh_x = (Nx-1/2.)*LatticeStep;
  //  BoxDimHigh_y = (Ny-1/2.)*LatticeStep;
  //  BoxDimHigh_z = (Nz+1/2.)*LatticeStep;
  
  double xmin = *min_element(rx.begin(), rx.end());
  double xmax = *max_element(rx.begin(), rx.end());
  double ymin = *min_element(ry.begin(), ry.end());
  double ymax = *max_element(ry.begin(), ry.end());
  double zmin = *min_element(rz.begin(), rz.end());
  dmin = *min_element(D.begin(), D.end());
  dmax = *max_element(D.begin(), D.end());
  
  BoxDimLow_x = xmin - dmax/2.;
  BoxDimHigh_x = xmax + dmax/2.;
  BoxDimLow_y = ymin - dmax/2.;
  BoxDimHigh_y = ymax + dmax/2.;
  BoxDimLow_z = zmin - dmax/2.;
  
  if (iWall == 1){
      if (iConfig == 2) {
        double zmax = *max_element(rz.begin(), rz.end());
        double LatticeSteptopWall = Dtopwall;
        double NxWall = int((BoxDimHigh_x - BoxDimLow_x + 1.0*Dtopwall)/LatticeSteptopWall);
        double NyWall = int((BoxDimHigh_y - BoxDimLow_y + 1.0*Dtopwall)/LatticeSteptopWall);
        double GapWall_x = ((BoxDimHigh_x - BoxDimLow_x) - NxWall*LatticeSteptopWall)/2.;
        double GapWall_y = ((BoxDimHigh_y - BoxDimLow_y) - NyWall*LatticeSteptopWall)/2.;
        for (size_t iy = 0; iy < NyWall; iy++){
          for (size_t ix = 0; ix < NxWall; ix++){
            rx.push_back(BoxDimLow_x+GapWall_x+LatticeSteptopWall/2.+ix*LatticeSteptopWall);
            ry.push_back(BoxDimLow_y+GapWall_y+LatticeSteptopWall/2.+iy*LatticeSteptopWall);
            rz.push_back(zmax + dmax/2. + LatticeSteptopWall/2.); // corrected on 2016.07.15
          }
        }
      } else {
//    if (iConfig == 2) {
//      string particle;
//      int grp;
//      double Rval, rxval, ryval, rzval;
//      string line;
//      istringstream iss;
//      ifstream myFile;
//      myFile.open (IutputFile_Config2);
//      getline(myFile, line);  // reading "sample{" in the first line
//      double NparTopWall = 0;
//      while(myFile.good()) {
//        getline(myFile, line);
//        iss.str(line);
//        iss >> particle;
//        if (particle == "sphere"){
//          iss >> grp >> Rval >> rxval >> ryval >> rzval;
//          if (grp==1) {
//            NparTopWall++;
//            rx.push_back(rxval);
//            ry.push_back(ryval);
//            rz.push_back(rzval);
//          }
//        }
//      }
//      myFile.close();
//    } else {
      double LatticeSteptopWall = Dtopwall;
      double NxWall = int((BoxDimHigh_x - BoxDimLow_x + 1.0*Dtopwall)/LatticeSteptopWall);
      double NyWall = int((BoxDimHigh_y - BoxDimLow_y + 1.0*Dtopwall)/LatticeSteptopWall);
      double GapWall_x = ((BoxDimHigh_x - BoxDimLow_x) - NxWall*LatticeSteptopWall)/2.;
      double GapWall_y = ((BoxDimHigh_y - BoxDimLow_y) - NyWall*LatticeSteptopWall)/2.;
      for (size_t iy = 0; iy < NyWall; iy++){
        for (size_t ix = 0; ix < NxWall; ix++){
          rx.push_back(BoxDimLow_x+GapWall_x+LatticeSteptopWall/2.+ix*LatticeSteptopWall);
          ry.push_back(BoxDimLow_y+GapWall_y+LatticeSteptopWall/2.+iy*LatticeSteptopWall);
          rz.push_back(Nz*LatticeStep);
        }
      }
//      double LatticeSteptopWall = Dtopwall;
//      double NxWall = int((BoxDimHigh_x - BoxDimLow_x + 0.0*dmin)/LatticeSteptopWall); //leaving at least 0.05*dmin distance from each end wall!
//      double NyWall = int((BoxDimHigh_y - BoxDimLow_y + 0.0*dmin)/LatticeSteptopWall);
//      double GapWall_x = ((BoxDimHigh_x - BoxDimLow_x) - NxWall*LatticeSteptopWall)/2.;
//      double GapWall_y = ((BoxDimHigh_y - BoxDimLow_y) - NyWall*LatticeSteptopWall)/2.;
//      for (size_t iy = 0; iy < NyWall; iy++){
//        for (size_t ix = 0; ix < NxWall; ix++){
//          rx.push_back(BoxDimLow_x+GapWall_x+LatticeSteptopWall/2.+ix*LatticeSteptopWall);
//          ry.push_back(BoxDimLow_y+GapWall_y+LatticeSteptopWall/2.+iy*LatticeSteptopWall);
//          rz.push_back(Nz*LatticeStep);
//        }
//      }
    }//iConfig
  }//iWall
  
  double zmax = *max_element(rz.begin(), rz.end());
  BoxDimHigh_z = zmax + dmax/2.;//*20.;
  
  ofstream newFile;
  //cout << OutputFile_SampleBoxSize;
  newFile.open (OutputFile_SampleBoxSize);
  newFile << "SampleBoxSize" << endl;
  newFile << BoxDimLow_x << "  " << BoxDimLow_y << "  " << BoxDimLow_z << endl;
  newFile << BoxDimHigh_x << "  " << BoxDimHigh_y << "  " << BoxDimHigh_z << endl;
  newFile.close();
}

// Read samples from an input file
void initCompaction1(){
  string dummy;
  size_t idummy;
  double Rval, rxval, ryval, rzval;
  string line;
  istringstream iss;
  
  ifstream myFile;
  myFile.open (IutputFile_Sample);
  getline(myFile, line);  // header of the input file
  Npar = 0;
  while(myFile.good()) {
    getline(myFile, line);
    iss.str(line);
    if (line == "}"){
      break;
    }
    iss >> dummy >> idummy >> Rval >> rxval >> ryval >> rzval;
    D.push_back(Rval*2.); // read radius of sample but store the diameter
    rx.push_back(rxval);
    ry.push_back(ryval);
    rz.push_back(rzval);
    Npar++;
    //cout << "rz = " << rzval << endl;
  }
  cout << "Npar = " << Npar << endl;
  myFile.close();
  
  // Find dmin (used for scaling in compaction routine)
  // Wall position
  
  double xmin = *min_element(rx.begin(), rx.end());
  double ymin = *min_element(ry.begin(), ry.end());
  double zmin = *min_element(rz.begin(), rz.end());
  double xmax = *max_element(rx.begin(), rx.end());
  double ymax = *max_element(ry.begin(), ry.end());
  double zmax = *max_element(rz.begin(), rz.end());
  dmin = *min_element(D.begin(), D.end());
  dmax = *max_element(D.begin(), D.end());
  
  BoxDimLow_x = xmin - dmax/2.;
  BoxDimLow_y = ymin - dmax/2.;
  BoxDimLow_z = zmin - dmax/2.;
  BoxDimHigh_x = xmax + dmax/2.;
  BoxDimHigh_y = ymax + dmax/2.;
  BoxDimHigh_z = zmax + dmax/2.;//*20.;
  CVDimHigh_z = kcv*dmin;
  
  Lx = BoxDimHigh_x - BoxDimLow_x;
  Ly = BoxDimHigh_y - BoxDimLow_y;
  Lz = BoxDimHigh_z - BoxDimLow_z;
  
  cout << "BoxDimLow_x = " << BoxDimLow_x << "   " << "BoxDimLow_y = " << BoxDimLow_y << "   " << "BoxDimLow_z = " << BoxDimLow_z << endl;
  cout << "BoxDimHigh_x = " << BoxDimHigh_x << "   " << "BoxDimHigh_y = " << BoxDimHigh_y << "   " << "BoxDimHigh_z = " << BoxDimHigh_z << endl;
  
  double E = 1e3;// kPa
  double nu = 0.3;
  double Etild = E/(1.-pow(nu,2));
  double epsilon = 0.01;// the smaller, the stiffer
  
  Ld = Kd * dmin;
  Ldmin = .2*Ld;
  Ldmax = Ld;
  LMao = KMao * dmin;
  Dver = DverFactor * dmin;
  
  double rmin = dmin/2;
  RepulsionFactor =  sqrt(rmin)/5.*Etild;
  GravityFactor = 0.;
  double p = 5.*pow(epsilon,1.5) / 2 / sqrt(2) * Etild;
  double S = Lx * Ly;
  ForceTopWall = p * S;
  
  double BetaFactor = 1.;
  double alpha = 1./2./sqrt(2) * pow(epsilon,1.5) * pow(dmin,2.) * Etild;
  BetaRef = BetaFactor * 0.7 / alpha / Ld;
  
  
  Dper = DperFactor * dmax;
  
  cout << "dmin = " << dmin << "  Dver = " << Dver << "   Ld = " << Ld << endl;
  
  
  //system("rmdir visualisation"); // only if the directory is empty :(
  system("mkdir visualisation");
  writevtk(0); // Initial position. Necessary if the first i loop is successful!
  VerletList();
  CalcPackingFraction();
  //CalcCoordinationNumber();
  phi_last = phi;
  ofstream newFile;
  newFile.open (OutputFile_MCprogress);
  newFile << "iStop" << "   " << "iMove" << "   " << "iRejectMove" << "   " << "iAcceptMove" << "   " << "RejectionRate" << "   " << "phi" << "   " << "z" << "   " << "MeanOverlapRatio" << "   " << "Ld" << endl;
  newFile.close();
}


// Read samples from an input file
void initCompaction2(){
  string dummy;
  size_t idummy;
  double Rval, rxval, ryval, rzval;
  string line;
  istringstream iss;
  
  ifstream myFile;
  myFile.open (IutputFile_Sample);
  getline(myFile, line);  // header of the input file
  Npar = 0;
  while(myFile.good()) {
    getline(myFile, line);
    iss.str(line);
    if (line == "}"){
      break;
    }
    iss >> dummy >> idummy >> Rval >> rxval >> ryval >> rzval;
    D.push_back(Rval*2.); // read radius of sample but store the diameter
    rx.push_back(rxval);
    ry.push_back(ryval);
    rz.push_back(rzval);
    Npar++;
    //cout << "rz = " << rzval << endl;
  }
  cout << "Npar = " << Npar << endl;
  myFile.close();
  
  // Find dmin (used for scaling in compaction routine)
  // Wall position
  
  double xmin = *min_element(rx.begin(), rx.end());
  double ymin = *min_element(ry.begin(), ry.end());
  double zmin = *min_element(rz.begin(), rz.end());
  double xmax = *max_element(rx.begin(), rx.end());
  double ymax = *max_element(ry.begin(), ry.end());
  double zmax = *max_element(rz.begin(), rz.end());
  dmin = *min_element(D.begin(), D.end());
  dmax = *max_element(D.begin(), D.end());
  
  BoxDimLow_x = xmin - dmax/2.;
  BoxDimLow_y = ymin - dmax/2.;
  BoxDimLow_z = zmin - dmax/2.;
  BoxDimHigh_x = xmax + dmax/2.;
  BoxDimHigh_y = ymax + dmax/2.;
  BoxDimHigh_z = zmax + dmax/2.;//*20.;
  CVDimHigh_z = kcv*dmin;
  
  Lx = BoxDimHigh_x - BoxDimLow_x;
  Ly = BoxDimHigh_y - BoxDimLow_y;
  Lz = BoxDimHigh_z - BoxDimLow_z;
  
  cout << "BoxDimLow_x = " << BoxDimLow_x << "   " << "BoxDimLow_y = " << BoxDimLow_y << "   " << "BoxDimLow_z = " << BoxDimLow_z << endl;
  cout << "BoxDimHigh_x = " << BoxDimHigh_x << "   " << "BoxDimHigh_y = " << BoxDimHigh_y << "   " << "BoxDimHigh_z = " << BoxDimHigh_z << endl;
  
  //-----------------
  //The difference between init/compaction1 and initCompaction2
  //needs gravity for both topwall and particles
  Ld = Kd * dmin;
  Ldmin = .2*Ld;
  Ldmax = Ld;
  LMao = KMao * dmin;
  Dver = DverFactor * dmin;
  
  GravityFactor = 1.;
  ForceTopWall = 1.;//p * S;
  
  double BetaFactor = 1.e10;
  BetaRef = BetaFactor * 0.7 / ForceTopWall / Ld;
  //-----------------
  
  Dper = DperFactor * dmax;
  
  cout << "dmin = " << dmin << "  Dver = " << Dver << "   Ld = " << Ld << endl;
  
  
  //system("rmdir visualisation"); // only if the directory is empty :(
  system("mkdir visualisation");
  writevtk(0); // Initial position. Necessary if the first i loop is successful!
  VerletList();
  CalcPackingFraction();
  //CalcCoordinationNumber();
  phi_last = phi;
  ofstream newFile;
  newFile.open (OutputFile_MCprogress);
  newFile << "iStop" << "   " << "iMove" << "   " << "iRejectMove" << "   " << "iAcceptMove" << "   " << "RejectionRate" << "   " << "phi" << "   " << "z" << "   " << "MeanOverlapRatio" << "   " << "Ld" << endl;
  newFile.close();
}


// Read samples from an input file
void initCompaction3(){
  string dummy;
  size_t idummy;
  double Rval, rxval, ryval, rzval;
  string line;
  istringstream iss;
  double ri,rzi,distance;
  size_t ipar;
  
  ifstream myFile;
  myFile.open (IutputFile_Sample);
  getline(myFile, line);  // header of the input file
  Npar = 0;
  while(myFile.good()) {
    getline(myFile, line);
    iss.str(line);
    if (line == "}"){
      break;
    }
    iss >> dummy >> idummy >> Rval >> rxval >> ryval >> rzval;
    D.push_back(Rval*2.); // read radius of sample but store the diameter
    rx.push_back(rxval);
    ry.push_back(ryval);
    rz.push_back(rzval);
    Npar++;
    //cout << "rz = " << rzval << endl;
  }
  cout << "Npar = " << Npar << endl;
  myFile.close();
  
  // Find dmin (used for scaling in compaction routine)
  // Wall position
  
  double xmin = *min_element(rx.begin(), rx.end());
  double ymin = *min_element(ry.begin(), ry.end());
  double zmin = *min_element(rz.begin(), rz.end());
  double xmax = *max_element(rx.begin(), rx.end());
  double ymax = *max_element(ry.begin(), ry.end());
  double zmax = *max_element(rz.begin(), rz.end());
  dmin = *min_element(D.begin(), D.end());
  dmax = *max_element(D.begin(), D.end());
  
  BoxDimLow_x = xmin - dmax/2.;
  BoxDimLow_y = ymin - dmax/2.;
  BoxDimLow_z = zmin - dmax/2.;
  BoxDimHigh_x = xmax + dmax/2.;
  BoxDimHigh_y = ymax + dmax/2.;
  BoxDimHigh_z = zmax + dmax/2.;//*20.;
  CVDimHigh_z = kcv*dmin;
  
  Lx = BoxDimHigh_x - BoxDimLow_x;
  Ly = BoxDimHigh_y - BoxDimLow_y;
  Lz = BoxDimHigh_z - BoxDimLow_z;
  
  cout << "BoxDimLow_x = " << BoxDimLow_x << "   " << "BoxDimLow_y = " << BoxDimLow_y << "   " << "BoxDimLow_z = " << BoxDimLow_z << endl;
  cout << "BoxDimHigh_x = " << BoxDimHigh_x << "   " << "BoxDimHigh_y = " << BoxDimHigh_y << "   " << "BoxDimHigh_z = " << BoxDimHigh_z << endl;
  
  double E = 1e3;// kPa
  double nu = 0.3;
  double Etild = E/(1.-pow(nu,2));
  double epsilon = 0.01;// the smaller, the stiffer
  
  Ld = Kd * dmin;
  Ldmax = Ld;
  Ldmin = .2*Ld;
  LMao = KMao * dmin;
  Dver = DverFactor * dmin;
  
  //  double rmin = dmin/2;
  //  RepulsionFactor =  sqrt(rmin)/5.*Etild;
  RepulsionFactor =  Etild/5.;
  GravityFactor = 0.;
  double p = 30.*pow(epsilon,1.5) / 2 / sqrt(2) * Etild;
  double S = Lx * Ly;
  ForceTopWall = p * S;
  
  double BetaFactor = 1.;
  double alpha = 1./2./sqrt(2) * pow(epsilon,1.5) * pow(dmin,2.) * Etild;
  BetaRef = BetaFactor * 0.7 / alpha / Ld;
  
  Dper = DperFactor * dmax;
  
  cout << "dmin = " << dmin << "  Dver = " << Dver << "   Ld = " << Ld << endl;
  
  
  //system("rmdir visualisation"); // only if the directory is empty :(
  system("mkdir visualisation");
  writevtk(0); // Initial position. Necessary if the first i loop is successful!
  VerletList();
  CalcPackingFraction();
  //CalcCoordinationNumber();
  phi_last = phi;
  //initializing UFactorOld
  double vlistWall6_size = vlistWall6.size();
  UFactorOld = 0.;
  for (size_t k = 0; k < vlistWall6_size; k++){
    ipar = vlistWall6.at(k);
    ri = D.at(ipar)/2.;
    rzi = rz.at(ipar);
    distance = abs(rzi - BoxDimHigh_z);
    if (distance < ri){
      UFactorOld += pow(ri-distance,2.5);
    }
  }
  MCprogress.open (OutputFile_MCprogress);
  MCprogress << "iStop" << "   " << "iMove" << "   " << "iRejectMove" << "   " << "iAcceptMove" << "   " << "RejectionRate" << "   " << "phi" << "   " << "z" << "   " << "MeanOverlapRatio" << "   " << "Ld" << "   " << "TopWallOverlapRatio" << endl;
}


// Create Verlet list
void VerletList(){
  double distance,distance1;
  double ri,rj,rxi,ryi,rzi,rxj0,ryj0,rxj,ryj,rzj;
  
  ipoint.clear();
  vlist.clear();
  vlistVx.clear();
  vlistVy.clear();
  vlistWall6.clear();
  
  // bring back into the main cell all real particles outside the cell
  for(size_t i = 0 ; i < Npar ; i++){
    ri = D.at(i)/2.;
    rxi = rx.at(i);
    ryi = ry.at(i);
    if (rxi < BoxDimLow_x){
      rx.at(i) = rxi+Lx;
    }
    if (ryi < BoxDimLow_y){
      ry.at(i) = ryi+Ly;
    }
    if (rxi > BoxDimHigh_x){
      rx.at(i) = rxi-Lx;
    }
    if (ryi > BoxDimHigh_y){
      ry.at(i) = ryi-Ly;
    }
  }
  
  for (size_t i = 0 ; i < Npar ; i++) {
    ri = D.at(i)/2.;
    size_t k = vlist.size();
    ipoint.push_back(k);
    rxi = rx.at(i);
    ryi = ry.at(i);
    rzi = rz.at(i);
    for (size_t j = 0 ; j < Npar ; j++) {
      rj = D.at(j)/2.;
      rxj0 = rx.at(j);
      ryj0 = ry.at(j);
      rzj  = rz.at(j);
      /*
       plan view:
       -------------
       | 7 | 4 | 8 |
       -------------
       | 1 | 0 | 2 |
       -------------
       | 5 | 3 | 6 |
       -------------
       */
      if (j != i) {
        //==================
        rxj = rxj0;
        ryj = ryj0;
        distance = sqrt(pow(rxi - rxj , 2.) +
                        pow(ryi - ryj , 2.) +
                        pow(rzi - rzj , 2.) ) -ri -rj;
        if (distance < Dver){
          vlist.push_back(j);
          vlistVx.push_back(0);
          vlistVy.push_back(0);
        } else {
          if((rxi>BoxDimLow_x+Dper && rxi<BoxDimHigh_x-Dper) && (ryi>BoxDimLow_y+Dper && ryi<BoxDimHigh_y-Dper)){//i in zone 0
            continue;
          } else {//i in one of zones 1-8
            if (rxi<BoxDimLow_x+Dper && (ryi>BoxDimLow_y+Dper && ryi<BoxDimHigh_y-Dper)){//i in zone 1; kx = -1, ky = 0
              rxj = rxj0 - Lx;
              ryj = ryj0;
              distance1 = sqrt(pow(rxi - rxj , 2.) +
                               pow(ryi - ryj , 2.) +
                               pow(rzi - rzj , 2.) ) -ri -rj;
              if (distance1 < Dver){
                vlist.push_back(j);
                vlistVx.push_back(-1);
                vlistVy.push_back(0);
              }
            } else if (rxi>BoxDimHigh_x-Dper && (ryi>BoxDimLow_y+Dper && ryi<BoxDimHigh_y-Dper)){//i in zone 2; kx = 1, ky = 0
              rxj = rxj0 + Lx;
              ryj = ryj0;
              distance1 = sqrt(pow(rxi - rxj , 2.) +
                               pow(ryi - ryj , 2.) +
                               pow(rzi - rzj , 2.) ) -ri -rj;
              if (distance1 < Dver){
                vlist.push_back(j);
                vlistVx.push_back(1);
                vlistVy.push_back(0);
              }
            } else if ((rxi>BoxDimLow_x+Dper && rxi<BoxDimHigh_x-Dper) && ryi<BoxDimLow_y+Dper){//i in zone 3; kx = 0, ky = -1
              rxj = rxj0;
              ryj = ryj0 - Ly;
              distance1 = sqrt(pow(rxi - rxj , 2.) +
                               pow(ryi - ryj , 2.) +
                               pow(rzi - rzj , 2.) ) -ri -rj;
              if (distance1 < Dver){
                vlist.push_back(j);
                vlistVx.push_back(0);
                vlistVy.push_back(-1);
              }
            } else if ((rxi>BoxDimLow_x+Dper && rxi<BoxDimHigh_x-Dper) && ryi>BoxDimHigh_y-Dper){//i in zone 4; kx = 0, ky = 1
              rxj = rxj0;
              ryj = ryj0 + Ly;
              distance1 = sqrt(pow(rxi - rxj , 2.) +
                               pow(ryi - ryj , 2.) +
                               pow(rzi - rzj , 2.) ) -ri -rj;
              if (distance1 < Dver){
                vlist.push_back(j);
                vlistVx.push_back(0);
                vlistVy.push_back(1);
              }
            } else {//i in zone 5-8
              for (int kx=-1;kx<2;kx++){
                for (int ky=-1;ky<2;ky++){
                  if (kx==0 && ky==0){
                    continue;
                  } else {
                    rxj = rxj0 + kx*Lx;
                    ryj = ryj0 + ky*Ly;
                    distance1 = sqrt(pow(rxi - rxj , 2.) +
                                     pow(ryi - ryj , 2.) +
                                     pow(rzi - rzj , 2.) ) -ri -rj;
                    if (distance1 < Dver){
                      vlist.push_back(j);
                      vlistVx.push_back(kx);
                      vlistVy.push_back(ky);
                    }
                  }//if (kx !=0 && ky !=0){
                }//for (ky=-1;ky<2;ky++){
              }//for (kx=-1;kx<2;kx++){
            }//else (i in zones 5-8)
          }//else (i in zones 1-8)
        }//else for Dver
        //==================
      }//if (j != i)
    }//j
    distance = abs(rzi - BoxDimLow_z) - ri;
    if (distance < Dver){
      vlist.push_back(Npar);
      vlistVx.push_back(0);
      vlistVy.push_back(0);
    }
    distance = abs(rzi - BoxDimHigh_z) - ri;
    if (distance < Dver){
      vlist.push_back(Npar+1);
      vlistVx.push_back(0);
      vlistVy.push_back(0);
    }
  }
  size_t k = vlist.size();
  ipoint.push_back(k);
  // creat vlistWall6 for the top wall neighbors
  for (size_t i = 0 ; i < Npar ; i++) {
    ri = D.at(i)/2.;
    rzi = rz.at(i);
    distance = abs(rzi - BoxDimHigh_z) - ri;
    if (distance < Dver){
      vlistWall6.push_back(i);
    }
  }
}


// MC main compaction loop
void  CompactionLoop1(){
  long int ipar,jpar;
  double theta_angle,phi_angle;//
  double nx,ny,nz;
  double rxnew,rynew,rznew,BoxDimHigh_znew;
  double distance,distancePrime;//
  double DeltaU,DeltaUFactor_Hertz;
  double ri,rj,rxi,ryi,rzi,rxj,ryj,rzj,riPrj;
  size_t bad_move = 0;
  double r_,pb; // Boltzman weight
  ldiv_t divresult;
  long int ivtk;
  long int iRejectMove;
  long int iAccept = 0; // cumulative number of success
  double RejectionRate;
  long int iAcceptMove;
  size_t Niter_phi;
  size_t vlistWall6_size;
  int kx,ky;
  double UFactor,UFactorPrime;
  int NiparLoop;
  
  ofstream newFile;
  newFile.open (OutputFile_MCprogress,ios::app);
  for (long int iStop = 0 ; iStop < NStop ; iStop++){
    Niter_phi = 0;
    sum_delta_phi = 0;
    for (long int iMove = 0 ; iMove < NMove ; iMove++){
      //==== Moving Wall6
      vlistWall6_size = vlistWall6.size();
      DeltaUFactor_Hertz = 0.;
      BoxDimHigh_znew = BoxDimHigh_z;
      //rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (-Ld); // in the range -Ld to 0
      rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2.*Ld - Ld; // in the range -Ld to Ld
      BoxDimHigh_znew = BoxDimHigh_z + rand_Ld;
      for (size_t k = 0; k < vlistWall6_size; k++){
        ipar = vlistWall6.at(k);
        ri = D.at(ipar)/2.;
        rzi = rz.at(ipar);
        distance = abs(rzi - BoxDimHigh_z);
        distancePrime = abs(rzi - BoxDimHigh_znew);
        if (distance > ri){
          UFactor = 0.;
        } else {
          UFactor = pow(ri-distance,2.5);
        }
        if (distancePrime > ri){
          UFactorPrime = 0.;
        } else {
          UFactorPrime = pow(ri-distancePrime,2.5);
        }
        DeltaUFactor_Hertz += UFactorPrime - UFactor;
      }
      r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
      //DeltaU = ForceTopWall * rand_Ld * pow((BoxDimHigh_z - BoxDimLow_z)/Lz,4) + RepulsionFactor * DeltaUFactor_Hertz;
      DeltaU = ForceTopWall * rand_Ld + RepulsionFactor * DeltaUFactor_Hertz;
      if (DeltaU < 0. ){
        BoxDimHigh_z = BoxDimHigh_znew;
      } else {
        pb = exp(-DeltaU*BetaRef);
        if (r_ < pb) {
          BoxDimHigh_z = BoxDimHigh_znew;
        }
      }
      //==== End of moving Wall6
      
      NiparLoop = 0;
      do {
        iRejectMove = 0;
        iAcceptMove = 0; // number of success for this iStop over all iMoves and all particles and all angles
        for (ipar = 0; ipar < Npar; ipar++){
          // find a particle i
          rxi = rx.at(ipar);
          ryi = ry.at(ipar);
          rzi = rz.at(ipar);
          ri = D.at(ipar)/2.;
          rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (Ld-Ldmin) + Ldmin; // in the range Ldmin to Ld
          for (long int iangle = 0; iangle < Nangle; iangle++){
            phi_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *M_PI; // in the range 0 to pi
            theta_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *2.*M_PI; // in the range 0 to 2pi
            nx = cos(theta_angle) * sin(phi_angle);
            ny = sin(theta_angle) * sin(phi_angle);
            nz = cos(phi_angle);
            // displace the particle
            rxnew = rxi + nx * rand_Ld;
            rynew = ryi + ny * rand_Ld;
            rznew = rzi + nz * rand_Ld;
            bad_move = 0;
            DeltaUFactor_Hertz = 0.;
            for (size_t k = ipoint.at(ipar) ; k < (ipoint.at(ipar+1)) ; k++){
              jpar = vlist.at(k);
              kx = vlistVx.at(k);
              ky = vlistVy.at(k);
              // calculate overlaps
              if (jpar < Npar){
                rxj = rx.at(jpar) + kx*Lx;
                ryj = ry.at(jpar) + ky*Ly;
                rzj = rz.at(jpar);
                rj = D.at(jpar)/2.;
                riPrj = ri + rj;
                distance = sqrt ( pow(rxi - rxj , 2.) + pow(ryi - ryj , 2.) + pow(rzi - rzj , 2.) );
                distancePrime = sqrt ( pow(rxnew - rxj , 2.) + pow(rynew - ryj , 2.) + pow(rznew - rzj , 2.) );
                if (distance > riPrj){
                  UFactor = 0.;
                } else {
                  UFactor = pow(riPrj-distance,2.5);
                }
                if (distancePrime > riPrj - LMao){
                  UFactorPrime = 0.;
                } else {
                  UFactorPrime = pow(riPrj-distancePrime,2.5);
                }
                DeltaUFactor_Hertz += UFactorPrime - UFactor;
              } else if (jpar == Npar) {
                distance = abs(rzi - BoxDimLow_z);
                distancePrime = abs(rznew - BoxDimLow_z);
                if (distance > ri){
                  UFactor = 0.;
                } else {
                  UFactor = pow(ri-distance,2.5);
                }
                if (distancePrime > ri){
                  UFactorPrime = 0.;
                } else {
                  UFactorPrime = pow(ri-distancePrime,2.5);
                }
                DeltaUFactor_Hertz += UFactorPrime - UFactor;
              } else if (jpar == Npar+1) {
                distance = abs(rzi - BoxDimHigh_z);
                distancePrime = abs(rznew - BoxDimHigh_z);
                if (distance > ri){
                  UFactor = 0.;
                } else {
                  UFactor = pow(ri-distance,2.5);
                }
                if (distancePrime > ri){
                  UFactorPrime = 0.;
                } else {
                  UFactorPrime = pow(ri-distancePrime,2.5);
                }
                DeltaUFactor_Hertz += UFactorPrime - UFactor;
              } else {
                cout << "=================== jpar is out of range!" << endl;
              }
            }
            DeltaU = GravityFactor * rand_Ld*nz + RepulsionFactor * DeltaUFactor_Hertz;
            if (DeltaU < 0. ){
              rx.at(ipar) = rxnew;
              ry.at(ipar) = rynew;
              rz.at(ipar) = rznew;
              break; // out of iangle loop
            } else {
              r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
              pb = exp(-DeltaU*BetaRef);
              if (r_ < pb) {
                rx.at(ipar) = rxnew;
                ry.at(ipar) = rynew;
                rz.at(ipar) = rznew;
                break; // out of iangle loop
              } else {
                bad_move = 1;
              }
            }
          } //iangle
          if (bad_move == 1){
            iRejectMove++;
          } else {
            iAccept++;
            iAcceptMove++;
          }
          // if necessary update VerletList
          divresult = ldiv(iAccept,Nver);
          if (divresult.rem == 0){
            VerletList();
          }
          divresult = ldiv(iAccept,Nvtk);
          if (divresult.rem == 0){
            ivtk = divresult.quot;
            writevtk(ivtk);
          }
        } //ipar
        NiparLoop++;
      } while (NiparLoop < Nrelax);// do while loop
      
      RejectionRate = iRejectMove*1./(iAcceptMove+iRejectMove);
      CalcPackingFraction();
      CalcCoordinationNumber();
      cout << "iStop = " << iStop << ",  "  <<
      "iMove = " << iMove << ",  " <<
      "iAccept = " << iAccept << ",  " <<
      "iRejMove = " << iRejectMove << ",  " <<
      //"iAcceptMove = " << iAcceptMove << ",  " <<
      "RejRate = " << RejectionRate  << ",  "  <<
      "phi = " << phi << ",  "  <<
      "z = " << z << ",  "  <<
      "MOR = " << MeanOverlapRatio << ",  "  <<
      //"TWOR = " << TopWallOverlapRatio << ",  "  <<
      //"OmegaTot = " << OmegaTot << ",  "  <<
      //"Ldmin/Ld = " << Ldmin/Ld <<
      "Ld = " << Ld <<
      endl;
      newFile << iStop << " ; " << iMove << " ; " << iRejectMove << " ; " << iAcceptMove << " ; " << RejectionRate << " ; " << phi << " ; " << z << " ; " << MeanOverlapRatio  << " ; " << Ld << endl;
      //      // Terminate if dphi is too samll ============
      //      delta_phi = phi - phi_last;
      //      if (delta_phi < 0.001){
      //        Niter_phi++;
      //        sum_delta_phi = sum_delta_phi + delta_phi;
      //      } else {
      //        Niter_phi = 0;
      //        sum_delta_phi = 0;
      //      }
      //      if (Niter_phi > 10 && (sum_delta_phi / Niter_phi) < 0.0001){
      //        cout << "delta_phi is too small -->" << endl;
      //        break;
      //      }
      //      phi_last = phi;
      //============================================
      // Terminate if dphi is too samll ============
      divresult = ldiv(iMove,dphiCheckN);//50); // check increase of phi every 10 times
      if (divresult.rem == 0){
        delta_phi = phi - phi_last;
        //        if (delta_phi < 0.000001){//< 0.00001){
        //          cout << "delta_phi is too small -->" << endl;
        //          break;
        //        }
        phi_last = phi;
      }
      //============================================
      
      if (RejectionRate > (ROverlap)){
        cout << "Overlap rate is too big -->" << endl;
        break;
      }
      // if necessary update Ldmin
      if (iRejectMove > 0.99*Npar){
        break;
      }
    } // iMove
    cout << "--> reduce the Ldmin" << endl;
    Ldmin = Ldmin/1.1;
    cout << "Ldmin is now reduced" << endl;
    //cout << "--> reduce the Ld" << endl;
    //Ld = Ld/1.1;//Ld/2.;
    //cout << "Ld is now reduced" << endl;
    Nver = Nver*1.1;//Nver*2;
    //BetaRef = .9*Ld;
    //BetaRef = BetaRef/1.1;
    // breathing!=====
    //for (ipar = 0; ipar<Npar; ipar++){
    //  rz.at(ipar) = DilationFactor * (rz.at(ipar) - BoxDimLow_z) + rz.at(ipar);
    //}
    //VerletList();
    //CalcPackingFraction();
    //CalcCoordinationNumber();
    //phi_last = phi;
    //================
    if (Ld < LMao){
      cout << "Ultimate achievable precision is reached by MC algorithm!" << endl;
      sit();
    }
  } // iStop
  newFile.close();
}



// MC main compaction loop
void  CompactionLoop2(){
  long int ipar,jpar;
  double theta_angle,phi_angle;//
  double nx,ny,nz;
  double rxnew,rynew,rznew,BoxDimHigh_znew;
  double distancePrime;//
  double DeltaU;
  double ri,rj,rxi,ryi,rzi,rxj,ryj,rzj,riPrj;
  size_t bad_move = 0;
  double r_,pb; // Boltzman weight
  ldiv_t divresult;
  long int ivtk;
  long int iRejectMove;
  long int iAccept = 0; // cumulative number of success
  double RejectionRate;
  long int iAcceptMove;
  size_t Niter_phi;
  size_t vlistWall6_size;
  int kx,ky;
  int NiparLoop;
  
  ofstream newFile;
  newFile.open (OutputFile_MCprogress,ios::app);
  for (long int iStop = 0 ; iStop < NStop ; iStop++){
    Niter_phi = 0;
    sum_delta_phi = 0;
    for (long int iMove = 0 ; iMove < NMove ; iMove++){
      //==== Moving Wall6
      vlistWall6_size = vlistWall6.size();
      BoxDimHigh_znew = BoxDimHigh_z;
      //rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (-Ld); // in the range -Ld to 0
      rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2.*Ld - Ld; // in the range -Ld to Ld
      BoxDimHigh_znew = BoxDimHigh_z + rand_Ld;
      bad_move = 0;
      for (size_t k = 0; k < vlistWall6_size; k++){
        ipar = vlistWall6.at(k);
        ri = D.at(ipar)/2.;
        rzi = rz.at(ipar);
        distancePrime = abs(rzi - BoxDimHigh_znew);
        if (distancePrime - ri < - LMao){
          bad_move = 1;
          break;
        }
      }
      if (bad_move == 0) {
        DeltaU = ForceTopWall * rand_Ld;
        if (DeltaU < 0. ){
          BoxDimHigh_z = BoxDimHigh_znew;
        } else {
          r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
          pb = exp(-DeltaU*BetaRef);
          if (r_ < pb) {
            BoxDimHigh_z = BoxDimHigh_znew;
          }
        }
      }
      //==== End of moving Wall6
      
      NiparLoop = 0;
      do {
        iRejectMove = 0;
        iAcceptMove = 0; // number of success for this iStop over all iMoves and all particles and all angles
        MeanOverlapRatio = 0.;
        NOverlap = 0;
        for (ipar = 0; ipar < Npar; ipar++){
          // find a particle i
          rxi = rx.at(ipar);
          ryi = ry.at(ipar);
          rzi = rz.at(ipar);
          ri = D.at(ipar)/2.;
          rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (Ld-Ldmin) + Ldmin; // in the range Ldmin to Ld
          for (long int iangle = 0; iangle < Nangle; iangle++){
            phi_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *M_PI; // in the range 0 to pi
            theta_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *2.*M_PI; // in the range 0 to 2pi
            nx = cos(theta_angle) * sin(phi_angle);
            ny = sin(theta_angle) * sin(phi_angle);
            nz = cos(phi_angle);
            // displace the particle
            rxnew = rxi + nx * rand_Ld;
            rynew = ryi + ny * rand_Ld;
            rznew = rzi + nz * rand_Ld;
            bad_move = 0;
            for (size_t k = ipoint.at(ipar) ; k < (ipoint.at(ipar+1)) ; k++){
              jpar = vlist.at(k);
              kx = vlistVx.at(k);
              ky = vlistVy.at(k);
              // calculate overlaps
              if (jpar < Npar){
                rxj = rx.at(jpar) + kx*Lx;
                ryj = ry.at(jpar) + ky*Ly;
                rzj = rz.at(jpar);
                rj = D.at(jpar)/2.;
                riPrj = ri + rj;
                distancePrime = sqrt ( pow(rxnew - rxj , 2.) + pow(rynew - ryj , 2.) + pow(rznew - rzj , 2.) );
                if (distancePrime - riPrj < -LMao){
                  bad_move=1;
                }
              } else if (jpar == Npar) {
                distancePrime = abs(rznew - BoxDimLow_z);
                if (distancePrime - ri < -LMao){
                  bad_move=1;
                }
              } else if (jpar == Npar+1) {
                distancePrime = abs(rznew - BoxDimHigh_z);
                if (distancePrime - ri < -LMao){
                  bad_move=1;
                }
              } else {
                cout << "=================== jpar is out of range!" << endl;
              }
            }
            if (bad_move==0){
              DeltaU = GravityFactor * rand_Ld*nz;
              if (DeltaU < 0. ){
                rx.at(ipar) = rxnew;
                ry.at(ipar) = rynew;
                rz.at(ipar) = rznew;
                break; // out of iangle loop
              } else {
                r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
                pb = exp(-DeltaU*BetaRef);
                if (r_ < pb) {
                  rx.at(ipar) = rxnew;
                  ry.at(ipar) = rynew;
                  rz.at(ipar) = rznew;
                  break; // out of iangle loop
                } else {
                  bad_move = 1;
                }
              }
            }
          } //iangle
          if (bad_move == 1){
            iRejectMove++;
          } else {
            iAccept++;
            iAcceptMove++;
          }
          // if necessary update VerletList
          divresult = ldiv(iAccept,Nver);
          if (divresult.rem == 0){
            VerletList();
          }
          divresult = ldiv(iAccept,Nvtk);
          if (divresult.rem == 0){
            ivtk = divresult.quot;
            writevtk(ivtk);
          }
        } //ipar
        NiparLoop++;
      } while (NiparLoop < Nrelax);// do while loop
      
      RejectionRate = iRejectMove*1./(iAcceptMove+iRejectMove);
      MeanOverlapRatio = MeanOverlapRatio/(NOverlap*dmin);
      CalcPackingFraction();
      CalcCoordinationNumber();
      cout << "iStop = " << iStop << ",  "  <<
      "iMove = " << iMove << ",  " <<
      "iAccept = " << iAccept << ",  " <<
      "iRejMove = " << iRejectMove << ",  " <<
      //"iAcceptMove = " << iAcceptMove << ",  " <<
      "RejRate = " << RejectionRate  << ",  "  <<
      "phi = " << phi << ",  "  <<
      "z = " << z << ",  "  <<
      "MOR = " << MeanOverlapRatio << ",  "  <<
      //"TWOR = " << TopWallOverlapRatio << ",  "  <<
      //"OmegaTot = " << OmegaTot << ",  "  <<
      //"Ldmin/Ld = " << Ldmin/Ld <<
      "Ld = " << Ld <<
      endl;
      newFile << iStop << " ; " << iMove << " ; " << iRejectMove << " ; " << iAcceptMove << " ; " << RejectionRate << " ; " << phi << " ; " << z << " ; " << MeanOverlapRatio  << " ; " << Ld << endl;
      //      // Terminate if dphi is too samll ============
      //      delta_phi = phi - phi_last;
      //      if (delta_phi < 0.001){
      //        Niter_phi++;
      //        sum_delta_phi = sum_delta_phi + delta_phi;
      //      } else {
      //        Niter_phi = 0;
      //        sum_delta_phi = 0;
      //      }
      //      if (Niter_phi > 10 && (sum_delta_phi / Niter_phi) < 0.0001){
      //        cout << "delta_phi is too small -->" << endl;
      //        break;
      //      }
      //      phi_last = phi;
      //============================================
      // Terminate if dphi is too samll ============
      divresult = ldiv(iMove,500);//50); // check increase of phi every 10 times
      if (divresult.rem == 0){
        delta_phi = phi - phi_last;
        //        if (delta_phi < 0.000001){//< 0.00001){
        //          cout << "delta_phi is too small -->" << endl;
        //          break;
        //        }
        phi_last = phi;
      }
      //============================================
      
      if (RejectionRate > (ROverlap)){
        cout << "Overlap rate is too big -->" << endl;
        break;
      }
      // if necessary update Ldmin
      if (iRejectMove > 0.99*Npar){
        break;
      }
    } // iMove
    cout << "--> reduce the Ldmin" << endl;
    Ldmin = Ldmin/1.1;
    cout << "Ldmin is now reduced" << endl;
    //cout << "--> reduce the Ld" << endl;
    //Ld = Ld/1.1;//Ld/2.;
    //cout << "Ld is now reduced" << endl;
    Nver = Nver*1.1;//Nver*2;
    //BetaRef = .9*Ld;
    //BetaRef = BetaRef/1.1;
    // breathing!=====
    //for (ipar = 0; ipar<Npar; ipar++){
    //  rz.at(ipar) = DilationFactor * (rz.at(ipar) - BoxDimLow_z) + rz.at(ipar);
    //}
    //VerletList();
    //CalcPackingFraction();
    //CalcCoordinationNumber();
    //phi_last = phi;
    //================
    if (Ld < LMao){
      cout << "Ultimate achievable precision is reached by MC algorithm!" << endl;
      sit();
    }
  } // iStop
  newFile.close();
}


// MC main compaction loop
int  CompactionLoop3(){
  long int ipar,jpar;
  double theta_angle,phi_angle;//
  double nx,ny,nz;
  double rxnew,rynew,rznew,BoxDimHigh_znew;
  double distance,distancePrime;//
  double DeltaU,DeltaUFactor_Hertz;
  double ri,rj,rxi,ryi,rzi,rxj,ryj,rzj,riPrj;
  size_t bad_move = 0;
  double r_,pb; // Boltzman weight
  ldiv_t divresult;
  long int ivtk;
  size_t Niter_phi;
  size_t vlistWall6_size;
  int kx,ky;
  double UFactor,UFactorPrime;
  int NiparLoop;
  double DeltaLd,DeltaR;
  
  for (iStop = 0 ; iStop < NStop ; iStop++){
    Niter_phi = 0;
    sum_delta_phi = 0;
    for (iMove = 0 ; iMove < NMove ; iMove++){
      //==== Moving Wall6
      if (iTopWallControl==0){ //Force control
        //Ld = Ldmax * exp(-TopWallOverlapRatio/(Kd));
        //Ld = Ldmax * exp(-RejectionRate/(4.*Kd));
        //Ld = Ldmax * (1.-RejectionRate);
        //=======================
        DeltaR = (RejectionRate - RejectionRateOld);
        if (DeltaR > 0){
          DeltaLd = - Ldmax * DeltaR * 2.;//Ld reaches 0 at RejRate=0.5
          RejectionRateOld = RejectionRate;
          Ld = Ld + DeltaLd;
          Ldmin = 0.2*Ld;
        } else {
          DeltaLd = 0;
        }
        //=======================
        Nver = Nver * (1.-DeltaLd/Ld);
        //=======================
        if (Ld < LMao){
          cout << "Ultimate achievable precision is reached by MC algorithm!" << endl;
          sit();
          //return 0;//end program
        }
        vlistWall6_size = vlistWall6.size();
        DeltaUFactor_Hertz = 0.;
        BoxDimHigh_znew = BoxDimHigh_z;
        //rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (-Ld); // in the range -Ld to 0
        rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * 2.*Ld - Ld; // in the range -Ld to Ld
        //rand_Ld  = Ld;
        BoxDimHigh_znew = BoxDimHigh_z + rand_Ld;
        UFactorNew = 0.;
        for (size_t k = 0; k < vlistWall6_size; k++){
          ipar = vlistWall6.at(k);
          ri = D.at(ipar)/2.;
          rzi = rz.at(ipar);
          distancePrime = abs(rzi - BoxDimHigh_znew);
          if (distancePrime < ri){
            UFactorNew += pow(ri-distancePrime,2.5)*sqrt(ri);
          }
        }
        r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
        //DeltaU = ForceTopWall * rand_Ld * pow((BoxDimHigh_z - BoxDimLow_z)/Lz,4) + RepulsionFactor * DeltaUFactor_Hertz;
        DeltaU = ForceTopWall * rand_Ld + RepulsionFactor * (UFactorNew - UFactorOld);
        if (DeltaU <= 0. ){
          BoxDimHigh_z = BoxDimHigh_znew;
          UFactorOld = UFactorNew;
        } else {
          pb = exp(-DeltaU*BetaRef);
          if (r_ < pb) {
            BoxDimHigh_z = BoxDimHigh_znew;
            UFactorOld = UFactorNew;
          }
        }
      } else {  //iTopWallControl==1 : displacement control
        //BoxDimHigh_z = BoxDimHigh_z - Ld;
        BoxDimHigh_z = BoxDimHigh_z - Ld * exp(-TopWallOverlapRatio/Kd + 1);
      }
      //==== End of moving Wall6
      
      NiparLoop = 0;
      do {
        iRejectMove = 0;
        iAcceptMove = 0; // number of success for this iStop over all iMoves and all particles and all angles
        for (ipar = 0; ipar < Npar; ipar++){
          // find a particle i
          rxi = rx.at(ipar);
          ryi = ry.at(ipar);
          rzi = rz.at(ipar);
          ri = D.at(ipar)/2.;
          //rand_Ld  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) * (Ld-Ldmin) + Ldmin; // in the range Ldmin to Ld
          //rand_Ld  = Ld;
          rand_Ld  = Ld*dmin/(2.*ri);
          for (long int iangle = 0; iangle < Nangle; iangle++){
            phi_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *M_PI; // in the range 0 to pi
            theta_angle = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)) *2.*M_PI; // in the range 0 to 2pi
            nx = cos(theta_angle) * sin(phi_angle);
            ny = sin(theta_angle) * sin(phi_angle);
            nz = cos(phi_angle);
            // displace the particle
            rxnew = rxi + nx * rand_Ld;
            rynew = ryi + ny * rand_Ld;
            rznew = rzi + nz * rand_Ld;
            bad_move = 0;
            DeltaUFactor_Hertz = 0.;
            for (size_t k = ipoint.at(ipar) ; k < (ipoint.at(ipar+1)) ; k++){
              jpar = vlist.at(k);
              kx = vlistVx.at(k);
              ky = vlistVy.at(k);
              // calculate overlaps
              if (jpar < Npar){
                rxj = rx.at(jpar) + kx*Lx;
                ryj = ry.at(jpar) + ky*Ly;
                rzj = rz.at(jpar);
                rj = D.at(jpar)/2.;
                riPrj = ri + rj;
                distancePrime = sqrt ( pow(rxnew - rxj , 2.) + pow(rynew - ryj , 2.) + pow(rznew - rzj , 2.) );
                if (distancePrime < riPrj - LMao){
                  bad_move=1;
                  break;
                }
              } else if (jpar == Npar) {
                distancePrime = abs(rznew - BoxDimLow_z);
                if (distancePrime < ri - LMao){
                  bad_move=1;
                  break;
                }
              } else if (jpar == Npar+1) {
                distance = abs(rzi - BoxDimHigh_z);
                distancePrime = abs(rznew - BoxDimHigh_z);
                if (distance > ri){
                  UFactor = 0.;
                } else {
                  UFactor = pow(ri-distance,2.5);
                }
                if (distancePrime > ri){
                  UFactorPrime = 0.;
                } else {
                  UFactorPrime = pow(ri-distancePrime,2.5);
                }
                DeltaUFactor_Hertz += (UFactorPrime - UFactor)*sqrt(ri);
              } else {
                cout << "=================== jpar is out of range!" << endl;
              }
            }
            if (bad_move==0){
              DeltaU = RepulsionFactor * DeltaUFactor_Hertz;
              if (DeltaU <= 0. ){
                rx.at(ipar) = rxnew;
                ry.at(ipar) = rynew;
                rz.at(ipar) = rznew;
                UFactorOld += DeltaUFactor_Hertz;
                break; // out of iangle loop
              } else {
                r_  = (static_cast <float> (rand()) / static_cast <float> (RAND_MAX)); // in the range 0 to 1
                pb = exp(-DeltaU*BetaRef);
                if (r_ < pb) {
                  rx.at(ipar) = rxnew;
                  ry.at(ipar) = rynew;
                  rz.at(ipar) = rznew;
                  UFactorOld += DeltaUFactor_Hertz;
                  break; // out of iangle loop
                } else {
                  bad_move = 1;
                }
              }
            }
          } //iangle
          if (bad_move == 1){
            iRejectMove++;
          } else {
            iAccept++;
            iAcceptMove++;
          }
          // if necessary update VerletList
          divresult = ldiv(iAccept,Nver);
          if (divresult.rem == 0){
            VerletList();
          }
          divresult = ldiv(iAccept,Nvtk);
          if (divresult.rem == 0){
            ivtk = divresult.quot;
            writevtk(ivtk);
          }
        } //ipar
        NiparLoop++;
      } while (NiparLoop < Nrelax);// do while loop
      
      RejectionRate = iRejectMove*1./(iAcceptMove+iRejectMove);
      MeanOverlapRatio = MeanOverlapRatio/(NOverlap*dmin);
      CalcPackingFraction();
      CalcCoordinationNumber();
      WriteScreen();
      WriteFile();
      readFlow();
      // Terminate if dphi is too samll ============
      //delta_phi = phi - phi_last;
      //if (delta_phi < 0.001){
      //  Niter_phi++;
      //  sum_delta_phi = sum_delta_phi + delta_phi;
      //} else {
      //  Niter_phi = 0;
      //  sum_delta_phi = 0;
      //}
      //if (Niter_phi > 10 && (sum_delta_phi / Niter_phi) < 0.0001){
      //  cout << "delta_phi is too small -->" << endl;
      //  break;
      //}
      //phi_last = phi;
      //      //============================================
      //      divresult = ldiv(iMove,dphiCheckN); // check increase of phi every dphiCheckN iMoves
      //      if (divresult.rem == 0){
      //        delta_phi = phi - phi_last;
      //        if (delta_phi < dphiCheckTol){
      //          cout << "* delta_phi is too small!" << endl;
      //          cout << "Ultimate achievable precision is reached by MC algorithm!" << endl;
      //          sit();
      //        }
      //        phi_last = phi;
      //      }
      //      //============================================
      
      //      if (RejectionRate > (ROverlap)){
      //        cout << "Overlap rate is too big -->" << endl;
      //        break;
      //      }
      //      // if necessary update Ldmin
      //      if (iRejectMove > 0.99*Npar){
      //        break;
      //      }
    } // iMove
    //    cout << "--> reduce the Ldmin" << endl;
    //    Ldmin = Ldmin/1.1;
    //    cout << "Ldmin is now reduced" << endl;
    //cout << "--> reduce the Ld" << endl;
    //Ld = Ld/1.1;//Ld/2.;
    //cout << "Ld is now reduced" << endl;
    //    Nver = Nver*1.1;//Nver*2;
    //BetaRef = .9*Ld;
    //BetaRef = BetaRef/1.1;
    // breathing!=====
    //for (ipar = 0; ipar<Npar; ipar++){
    //  rz.at(ipar) = DilationFactor * (rz.at(ipar) - BoxDimLow_z) + rz.at(ipar);
    //}
    //VerletList();
    //CalcPackingFraction();
    //CalcCoordinationNumber();
    //phi_last = phi;
    //================
    //    if (Ld < LMao){
    //      cout << "Ultimate achievable precision is reached by MC algorithm!" << endl;
    //      break;
    //    }
  } // iStop
  return 0;
}


// Write samples in an output file
void writeSample(){
  
  ofstream newFile;
  //cout << OutputFile_sample;
  newFile.open (OutputFile_sample);
  newFile << "sample{" << endl;
  for (size_t ipar = 0; ipar < Npar; ipar++){
      newFile << "sphere" << " " << 0 << " " << D.at(ipar)/2. << " "
              << rx.at(ipar) << " " << ry.at(ipar) << " " << rz.at(ipar) << " "
              << 0 << " " << 0 << " " << 0 << " "
              << 0 << " " << 0 << " " << 0 << " "
              << 0 << " " << 0 << " " << 0 << " "
              << 0 << " " << 0 << " " << 0 << " "
              << 0 << " " << 0 << " " << 0 << endl;
  }
  if (iWall == 1){
      for (size_t ipar = Npar; ipar < rx.size(); ipar++){
          newFile << "sphere" << " " << 1 << " " << Dtopwall/2. << " "
                  << rx.at(ipar) << " " << ry.at(ipar) << " " << rz.at(ipar) << " "
                  << 0 << " " << 0 << " " << 0 << " "
                  << 0 << " " << 0 << " " << 0 << " "
                  << 0 << " " << 0 << " " << 0 << " "
                  << 0 << " " << 0 << " " << 0 << " "
                  << 0 << " " << 0 << " " << 0 << endl;
      //        cout << "rx=" << rx.at(ipar) << "    " << endl;
    }
  }
  newFile << "}" << endl;
  newFile.close();
//
  //cout << OutputFile_sample_radius;
  newFile.open (OutputFile_sample_radius);
  for (size_t ipar = 0; ipar < Npar; ipar++){
    newFile << D.at(ipar)/2. << endl;
  }
  newFile.close();
}

// Write config in an output file
void writeConfig(){
  ofstream newFile;
  //cout << OutputFile_sample;
  newFile.open (OutputFile_config);
  newFile << "Max_Dist= " << Dmin*SizeRatio << endl;
  newFile << "rx   ry    rz" << endl;
  for (size_t ipar = 0; ipar < Npar; ipar++){
    newFile << rx.at(ipar) << " " << ry.at(ipar) << " " << rz.at(ipar) << " " << endl;
  }
  newFile.close();
}

void writevtk(
              const long int ivtk
              )
{
  char filename[100];
  sprintf(filename,"visualisation/sphere_%.6li.vtk",ivtk);
  ofstream vtk_spl(filename,ios::out);
  
  vtk_spl.precision(15);
  vtk_spl.setf(ios::fixed,ios::floatfield); //fixes the default format for
  
  vtk_spl << "# vtk DataFile Version 4.0\nSortie Particles\nASCII" << endl;
  vtk_spl << "DATASET UNSTRUCTURED_GRID" << endl;
  vtk_spl << "POINTS" << " " << D.size() << " " << "float" << endl;
  
  for (unsigned long int i = 0; i < D.size(); i++)
    vtk_spl << rx.at(i) << " " << ry.at(i) << " " << rz.at(i) << endl;
  
  vtk_spl << "POINT_DATA" << " " << D.size() << endl;
  vtk_spl << "VECTORS Radius float" << endl;
  
  for (unsigned long int i = 0 ; i < D.size(); i++)
    vtk_spl << 0. << " " << 0. << " " << D.at(i)/2.0 << endl; // writing radius
  
  vtk_spl.close();
  
  
  sprintf(filename,"visualisation/wall_%.6li.vtk",ivtk);
  ofstream vtk_wall(filename,ios::out);
  
  vtk_wall << "# vtk DataFile Version 4.0" << endl << "Wall" << endl << "ASCII" << endl;
  vtk_wall << "DATASET POLYDATA" << endl;
  vtk_wall << "POINTS 8 float" << endl;
  
  vtk_wall << BoxDimLow_x << " " << BoxDimLow_y << " " << BoxDimLow_z << endl;
  vtk_wall << BoxDimHigh_x << " " << BoxDimLow_y << " " << BoxDimLow_z << endl;
  vtk_wall << BoxDimHigh_x << " " << BoxDimHigh_y << " " << BoxDimLow_z << endl;
  vtk_wall << BoxDimLow_x << " " << BoxDimHigh_y << " " << BoxDimLow_z << endl;
  
  vtk_wall << BoxDimLow_x << " " << BoxDimLow_y << " " << BoxDimHigh_z << endl;
  vtk_wall << BoxDimHigh_x << " " << BoxDimLow_y << " " << BoxDimHigh_z << endl;
  vtk_wall << BoxDimHigh_x << " " << BoxDimHigh_y << " " << BoxDimHigh_z << endl;
  vtk_wall << BoxDimLow_x << " " << BoxDimHigh_y << " " << BoxDimHigh_z << endl;
  
  vtk_wall << "POLYGONS 5 25" << endl;
  
  vtk_wall << "4 0 1 2 3" << endl;
  vtk_wall << "4 0 4 5 1" << endl;
  vtk_wall << "4 1 5 6 2" << endl;
  vtk_wall << "4 2 6 7 3" << endl;
  vtk_wall << "4 0 3 7 4" << endl;
  //vtk_wall << "4 4 5 6 7" << endl;
  
  vtk_wall.close();
}


//// Write MCprogress in an output file
//void writeMCprogress(const long int& iStop, const long int& iMove, const long int& iRejectMove, const long int& iAcceptMove, const double& RejectionRate, const double& phi, const double& z,int flag){
//
//  ofstream newFile;
//  if (flag == 0){
//    newFile.open (OutputFile_MCprogress);
//    newFile << "iStop" << "   " << "iMove" << "   " << "iRejectMove" << "   " << "iAcceptMove" << "   " << "RejectionRate" << "   " << "phi" << "   " << "z" << endl;
//  } else {
//    newFile << iStop << "   " << iMove << "   " << iRejectMove << "   " << iAcceptMove << "   " << RejectionRate << "   " << phi << "   " << z << endl;
//    newFile.close();
//  }
//}


// Calculation of beta function
double betacf(double a, double b, double x){
  //void nrerror(char error_text[]);
  int m,m2;
  double aa,c,d,del,h,qab,qam,qap;
  qab=a+b;
  qap=a+1.0;
  qam=a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN; d=1.0/d;
  h=d;
  for (m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if (fabs(c) < FPMIN) c=FPMIN; d=1.0/d;
    del=d*c;
    h *= del;
    if (fabs(del-1.0) < EPS) break;
  }
  if (m > MAXIT) perror("a or b too big, or MAXIT too small in betacf");
  return h;
}

double gammln(double xx){
  double x,y,tmp,ser;
  static double cof[6]={76.18009172947146,-86.50532032941677,
    24.01409824083091,-1.231739572450155,
    0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  y=x=xx;
  tmp=x+5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for (j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

double betai(double a, double b, double x){
  double betacf(double a, double b, double x);
  double gammln(double xx);
  //void nrerror(char error_text[]);
  double bt;
  if (x < 0.0 || x > 1.0) perror("Bad x in routine betai");
  if (x == 0.0 || x == 1.0) bt=0.0;
  else
    bt=exp(gammln(a+b)-gammln(a)-gammln(b)+a*log(x)+b*log(1.0-x));
  if (x < (a+1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

// clock tic toc functions
void tic() {
  tictoc_stack.push(clock());
}
void toc() {
  std::cout << "Time elapsed: "
  << ((double)(clock() - tictoc_stack.top())) / CLOCKS_PER_SEC
  << " Sec." << std::endl;
  tictoc_stack.pop();
}

//// Postprocessing function for calculation of packing fraction
//void CalcPackingFraction(){
//  double Vcv;
//  double Vs = 0;
//  double radius;
//  double h;
//
//  // calculate the Vcv = total volume of the CV
//  //Vcv = (BoxDimHigh_x - BoxDimLow_x) * (BoxDimHigh_y - BoxDimLow_y) * (CVDimHigh_z - BoxDimLow_z);
//  Vcv = (BoxDimHigh_x - BoxDimLow_x) * (BoxDimHigh_y - BoxDimLow_y) * (BoxDimHigh_z - BoxDimLow_z);
//
//  for (size_t ipar = 0; ipar < Npar; ipar++){
//    radius = D.at(ipar)/2.;
//    Vs = Vs + pow(radius,3.);
//  } // ipar loop
//  phi = 4.*M_PI/3.*Vs /Vcv;
//}


// Postprocessing function for calculation of packing fraction
void CalcPackingFraction(){
  double Vcv;
  double Vs = 0;
  double Vs_overlap = 0;
  size_t jpar;
  double rxi,ryi,rzi,rxj,ryj,rzj;
  double ri,rj,riPrj;
  double h,d;
  int kx,ky;
  
  // calculate the Vcv = total volume of the CV
  Vcv = (BoxDimHigh_x - BoxDimLow_x) * (BoxDimHigh_y - BoxDimLow_y) * (BoxDimHigh_z - BoxDimLow_z);
  
  MeanOverlapRatio = 0.;
  NOverlap = 0;
  TopWallOverlapRatio = 0.;
  NTopWallOverlap = 0.;
  for (size_t ipar = 0; ipar < Npar; ipar++){
    ri = D.at(ipar)/2.;
    rxi = rx.at(ipar);
    ryi = ry.at(ipar);
    rzi = rz.at(ipar);
    Vs += pow(ri,3.);
    for (size_t k = ipoint.at(ipar) ; k < (ipoint.at(ipar+1)) ; k++){
      jpar = vlist.at(k);
      kx = vlistVx.at(k);
      ky = vlistVy.at(k);
      // calculate overlaps
      if (jpar < Npar){
        rj = D.at(jpar)/2.;
        rxj = rx.at(jpar) + kx*Lx;
        ryj = ry.at(jpar) + ky*Ly;
        rzj = rz.at(jpar);
        riPrj = ri+rj;
        d = sqrt(pow(rxi - rxj , 2.) +
                 pow(ryi - ryj , 2.) +
                 pow(rzi - rzj , 2.) );
        if (d < riPrj){
          Vs_overlap += pow(ri+rj-d,2.) * ( pow(d,2.) + 2.*d*ri + 2.*d*rj - 3.*pow(ri,2.) - 3.*pow(rj,2.) + 6*ri*rj ) / 24. / d; // already devided by 2
          MeanOverlapRatio += riPrj-d;
          NOverlap++;
        }
      } else if (jpar == Npar) {
        d = abs(rzi - BoxDimLow_z);
        if (d < ri) {
          h = ri-d;
          Vs_overlap += 1/3*pow(h,2)*(3*ri-h);
          MeanOverlapRatio += ri-d;
          NOverlap++;
        }
        if (rzi < BoxDimLow_z){
          cerr << "Particle passed the bottom wall!" << endl;
          //          sit();
        }
      } else if (jpar == Npar+1) {
        d = abs(rzi - BoxDimHigh_z);
        if (d < ri) {
          h = ri-d;
          Vs_overlap += 1/3*pow(h,2)*(3*ri-h);
          //MeanOverlapRatio += ri-d;
          //NOverlap++;
          TopWallOverlapRatio += h;
          NTopWallOverlap++;
        }
        if (rzi > BoxDimHigh_z){
          cerr << "Particle passed the top wall!" << endl;
          //          sit();
        }
      } else {
        cout << "=================== jpar is out of range!" << endl;
      }
    }
  } // ipar loop
  phi = ((4.*M_PI/3.*Vs) - (M_PI*Vs_overlap)) / Vcv;
  if(NOverlap==0){
    MeanOverlapRatio = 0.;
  } else {
    MeanOverlapRatio = MeanOverlapRatio/(NOverlap*dmin);
  }
  if(NTopWallOverlap==0){
    TopWallOverlapRatio = 0.;
  } else {
    TopWallOverlapRatio = TopWallOverlapRatio/(NTopWallOverlap*dmin);
  }
  
}

//// Postprocessing function for calculation of packing fraction
//void CalcPackingFraction(){
//  double Vcv;
//  double Vs = 0;
//  double dz,z1;
//  double radius;
//  double h;
//
//  // calculate the Vcv = total volume of the CV
//  //Vcv = (BoxDimHigh_x - BoxDimLow_x) * (BoxDimHigh_y - BoxDimLow_y) * (CVDimHigh_z - BoxDimLow_z);
//  Vcv = (BoxDimHigh_x - BoxDimLow_x) * (BoxDimHigh_y - BoxDimLow_y) * (BoxDimHigh_z - BoxDimLow_z);
//
//  for (size_t ipar = 0; ipar < Npar; ipar++){
//    radius = D.at(ipar)/2.;
//    z1 = rz.at(ipar);
//    dz = CVDimHigh_z - z1;
//    if (dz >= 0 ){
//      if (dz - radius >=0){
//        h = 2*radius;
//      } else {
//        h = dz + radius;
//      }
//    } else {
//      if (dz + radius <=0){
//        h = 0.;
//      } else {
//        h = dz + radius;
//      }
//    }
//    Vs = Vs + M_PI*pow(h,2.)/3*(3*radius-h);
//  } // ipar loop
//  phi = Vs /Vcv;
//}

// Postprocessing function for calculation of coordination number
void CalcCoordinationNumber(){
  size_t Ncontact_cv = 0;
  size_t jpar;
  double rxi,ryi,rzi,rxj,ryj,rzj;
  double ri,rj;
  double distance;
  int kx,ky;
  
  for (size_t ipar = 0; ipar < Npar; ipar++){
    ri = D.at(ipar)/2.;
    rxi = rx.at(ipar);
    ryi = ry.at(ipar);
    rzi = rz.at(ipar);
    for (size_t k = ipoint.at(ipar) ; k < (ipoint.at(ipar+1)) ; k++){
      jpar = vlist.at(k);
      kx = vlistVx.at(k);
      ky = vlistVy.at(k);
      // calculate overlaps
      if (jpar < Npar){
        rj = D.at(jpar)/2.;
        rxj = rx.at(jpar) + kx*Lx;
        ryj = ry.at(jpar) + ky*Ly;
        rzj = rz.at(jpar);
        distance = sqrt(pow(rxi - rxj , 2.) +
                        pow(ryi - ryj , 2.) +
                        pow(rzi - rzj , 2.) ) - ri - rj;
        if (distance < Ld/2.) {Ncontact_cv++;}
        //if (distance < 0.) {Ncontact_cv++;}
      } else if (jpar == Npar) {
        distance = abs(rzi - BoxDimLow_z) - ri;
        if (distance < Ld/2.) {Ncontact_cv++;}
        //if (distance < 0.) {Ncontact_cv++;}
      } else if (jpar == Npar+1) {
        distance = abs(rzi - BoxDimHigh_z) - ri;
        if (distance < Ld/2.) {Ncontact_cv++;}
        //if (distance < 0.) {Ncontact_cv++;}
      } else {
        cout << "=================== jpar is out of range!" << endl;
      }
    }//k loop
  } //ipar loop
  z = double(Ncontact_cv) / double(Npar);
}

vector<double> linspace(double a, double b, int n) {
  vector<double> array;
  double step = (b-a) / (n-1);
  
  while(a <= b) {
    array.push_back(a);
    a += step;           // could recode to better handle rounding errors
  }
  return array;
}

void WriteScreen(){
  if (iWriteScreen == 1){
    ldiv_t divresult = ldiv(iMove,nWriteScreen);
    if (divresult.rem == 0){
      cout << "iStop = " << iStop << ",  "  <<
      "iMove = " << iMove << ",  " <<
      "iAccept = " << iAccept << ",  " <<
      "iRejMove = " << iRejectMove << ",  " <<
      //"iAcceptMove = " << iAcceptMove << ",  " <<
      "RejRate = " << RejectionRate  << ",  "  <<
      "phi = " << phi << ",  "  <<
      "z = " << z << ",  "  <<
      "MOR = " << MeanOverlapRatio << ",  "  <<
      "TWOR = " << TopWallOverlapRatio << ",  "  <<
      //"OmegaTot = " << OmegaTot << ",  "  <<
      //"Ldmin/Ld = " << Ldmin/Ld <<
      "Ld = " << Ld <<
      endl;
    }
  }
}

void WriteFile(){
  if (iWriteFile == 1){
    ldiv_t divresult = ldiv(iMove,nWriteFile);
    if (divresult.rem == 0){
      MCprogress << iStop << " ; " << iMove << " ; " << iRejectMove << " ; " << iAcceptMove << " ; " << RejectionRate << " ; " << phi << " ; " << z << " ; " << MeanOverlapRatio  << " ; " << Ld << " ; " << TopWallOverlapRatio << endl;
    }
  }
}


int sit(){
  writeSample();
  writeConfig();
  MCprogress.close();
  toc();
  cerr << "Program ends in sit." << endl;
  exit(1);
  return 0;
}

// Main program
int main (){
  tic();
  header();
  cout << "--------- Reading input parameters ----------" << endl;
  srand((unsigned)time(0));
  readParameters();
  readFlow();
  if (iCompaction == 0){
    if (iBetaFn == 1){
    size_t Npar_init = Npar;
    ofstream newFile;
    newFile.open (OutputFile_BetaFnSensitivity);
    newFile << "# OutputFile_BetaFnSensitivity: a, b, N_p^min" << endl;
    vector<double> aVector = linspace(amin,amax,200);
    vector<double> bVector = linspace(bmin,bmax,200);
      for (unsigned i=0; i < aVector.size(); i++){
        for (unsigned j=0; j < bVector.size(); j++){
          a = aVector.at(i);
          b = bVector.at(j);
          setParticleSize();
          newFile << a << ";" << b << ";" << Npar << endl;
          Npar = Npar_init;
        }
      }
      newFile.close();
      exit(1);
    }
    setParticleSize();
    setParticlePosition();
  }
  else if (iCompaction == 1){
    initCompaction1();
    CompactionLoop1();
  }
  else if (iCompaction == 2){
    initCompaction2();
    CompactionLoop2();
  }
  else if (iCompaction == 3){
    initCompaction3();
    CompactionLoop3();
  }
  sit();
}

