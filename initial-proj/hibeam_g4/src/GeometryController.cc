#include "GeometryController.hh"

#include "G4RunManager.hh"
#include "ConvertGeo.hh"
#include "G4VUserParallelWorld.hh"
#include "WasaDetectorConstruction.hh"

GeometryController::GeometryController(Config& _par) : Par(_par) {}

GeometryController::~GeometryController() {}

int GeometryController::SetGeometry(G4String nameGeo)
{
	G4cout << "Activating geometry " << nameGeo << G4endl;
	if(nameGeo == "Wasa")
	{
		detectorBuilder = new WasaDetectorConstruction(Par);
		registerGeometry(detectorBuilder);
	}
	else
	{
		std::cout << "E> No Geometry selected !" << nameGeo << "\n";
		return -2;
	}

	nameGeometry = nameGeo;
	return 0;
}

void GeometryController::registerGeometry(G4VUserDetectorConstruction* detector)
{
  G4RunManager* runManager = G4RunManager::GetRunManager();
  runManager->SetUserInitialization(detector);
  runManager->GeometryHasBeenModified();
}

const std::vector<G4String>& GeometryController::GetNameDetectors() { return detectorBuilder->GetNameDetectors(); }

void GeometryController::ConvertG4toRoot(const std::string& nameConvertRoot)
{

  auto* logi = detectorBuilder->experimentalHall_physOutRoot->GetLogicalVolume();
  if(logi != nullptr)
    {
      std::cout << logi->GetName() << std::endl;
      ConvertGeo convertor(detectorBuilder->experimentalHall_physOutRoot);
      convertor.Convert(nameConvertRoot, nameGeometry, detectorBuilder->GetNameDetectors(), Par);
    }
  else
    {
      std::cout << "E> Convertion G4toRoot failed!\n";
      return;
    }
}
