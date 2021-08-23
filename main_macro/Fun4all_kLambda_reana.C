#include <G4_UserkLambda.C>
#include <TROOT.h>
#include <fun4all/Fun4AllDstInputManager.h>
#include <fun4all/Fun4AllDstOutputManager.h>
#include <fun4all/Fun4AllOutputManager.h>
#include <fun4all/Fun4AllServer.h>
#include <phool/recoConsts.h>

R__LOAD_LIBRARY(libfun4all.so)

int Fun4all_kLambda_reana(){

  gSystem->Load("libfun4all.so");
  gSystem->Load("libg4dst.so");

  Enable::USER = true;
    
  Fun4AllServer *se = Fun4AllServer::instance();

  se->Verbosity(0); 
    
  Fun4AllInputManager *hitsin= new Fun4AllDstInputManager("DSTin");

  hitsin->AddListFile("mykLambdaFileList.txt");

  // #Add your analysis modules here
    
  se->registerInputManager(hitsin);

  if (Enable::USER) UserAnalysisInit();

  Int_t nEvents = 25000;

  se->run(nEvents);

  se->End();
    
  return 0;

}
