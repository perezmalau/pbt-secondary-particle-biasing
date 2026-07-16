//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file GB01/exampleGB01.cc
/// \brief Main program of the GB01 example
//
//

#include "FTFP_BERT.hh"
#include "QGSP_BIC.hh"
#include "G4EmStandardPhysics_option4.hh"
#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGeneratorAction.hh"

#include "G4StepLimiterPhysics.hh"

#include "G4ScoringManager.hh"

#include "G4GenericBiasingPhysics.hh"
#include "G4RunManagerFactory.hh"
#include "G4Types.hh"
#include "G4UIExecutive.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#ifdef G4MULTITHREADED
#include "G4MTRunManager.hh"
#else
#include "G4RunManager.hh"
#endif

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

namespace
{
void PrintUsage()
{
  G4cerr << " Usage: " << G4endl;
  G4cerr << " ./biased_pbt [-m macro ] "
         << " [-b biasing {'on','off'}]"
         << "\n or\n ./biased_pbt [macro.mac]" << G4endl;
}
}  // namespace

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc, char** argv)
{
  // Evaluate arguments
  //
  if (argc > 7) {
    PrintUsage();
    return 1;
  }


  G4String macro("");
  G4String onOffBiasing("");
  G4int nThreads = 4;  // default

  if (argc == 2)
    macro = argv[1];
  else {
    for (G4int i = 1; i < argc; i = i + 2) {
      if (G4String(argv[i]) == "-m")
        macro = argv[i + 1];
      else if (G4String(argv[i]) == "-b")
        onOffBiasing = argv[i + 1];
      else if (G4String(argv[i]) == "-t")
        nThreads = G4UIcommand::ConvertToInt(argv[i + 1]);
      else {
        PrintUsage();
        return 1;
      }
    }
  }

  if (onOffBiasing == "") onOffBiasing = "off";

  // Instantiate G4UIExecutive if interactive mode
  G4UIExecutive* ui = nullptr;
  if (macro == "") {
    ui = new G4UIExecutive(argc, argv);
  }

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine());
  auto seed = (G4int) time(nullptr);
  // Possibility to set an offset to random number via environment variable "JOB_ID"
  char *jobID = getenv("JOB_ID");
  if (jobID != nullptr)
    seed += atoi(jobID) * 4852234 * G4UniformRand();
  G4cout << "Seed:" << seed << G4endl;
  G4Random::setTheSeed(seed);


  // -- Construct the run manager : MT or sequential one
#ifdef G4MULTITHREADED
  auto *runManager = new G4MTRunManager;
  runManager->SetNumberOfThreads(nThreads);
#else
  G4RunManager *runManager = new G4RunManager;
#endif

  // Activate UI-command base scorer
  G4ScoringManager * scManager = G4ScoringManager::GetScoringManager();
  scManager->SetVerboseLevel(1);

  // -- Set mandatory initialization classes
  DetectorConstruction* detector = new DetectorConstruction();
  runManager->SetUserInitialization(detector);
  // -- Select a physics list:
  //FTFP_BERT* physicsList = new FTFP_BERT;
  QGSP_BIC* physicsList = new QGSP_BIC;
  physicsList->ReplacePhysics(new G4EmStandardPhysics_option4());

  physicsList->SetVerboseLevel(0);
  physicsList->RegisterPhysics(new G4StepLimiterPhysics());
  runManager->SetUserInitialization(physicsList);

  // -- and augment it with biasing facilities:
  G4GenericBiasingPhysics* biasingPhysics = new G4GenericBiasingPhysics();
  if (onOffBiasing == "on") {
    biasingPhysics->PhysicsBias("proton");
    physicsList->RegisterPhysics(biasingPhysics);
    G4cout << "      ********************************************************* " << G4endl;
    G4cout << "      ***** proton inelastic processes wrapped for biasing **** " << G4endl;
    G4cout << "      ********************************************************* " << G4endl;
  }
  else {
    G4cout << "      ************************************************* " << G4endl;
    G4cout << "      ********** processes are not wrapped ************ " << G4endl;
    G4cout << "      ************************************************* " << G4endl;
  }
  runManager->SetUserInitialization(physicsList);
  // -- Action initialization:
  runManager->SetUserInitialization(new ActionInitialization);

  physicsList->SetVerboseLevel(0);
  // Initialize G4 kernel
  runManager->Initialize();

  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  visManager->Initialize();

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (!ui)  // batch mode
  {
    G4String command = "/control/execute ";
    UImanager->ApplyCommand(command + macro);
  }
  else {  // interactive mode
    UImanager->ApplyCommand("/control/execute vis.mac");
    ui->SessionStart();
    delete ui;
  }

  delete visManager;
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
