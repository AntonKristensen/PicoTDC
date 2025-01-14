// Scorer for LindhardDamage
//
// ********************************************************************
// *                                                                  *
// * Copyright 2024 The TOPAS Collaboration                           *
// * Copyright 2022 The TOPAS Collaboration                           *
// *                                                                  *
// * Permission is hereby granted, free of charge, to any person      *
// * obtaining a copy of this software and associated documentation   *
// * files (the "Software"), to deal in the Software without          *
// * restriction, including without limitation the rights to use,     *
// * copy, modify, merge, publish, distribute, sublicense, and/or     *
// * sell copies of the Software, and to permit persons to whom the   *
// * Software is furnished to do so, subject to the following         *
// * conditions:                                                      *
// *                                                                  *
// * The above copyright notice and this permission notice shall be   *
// * included in all copies or substantial portions of the Software.  *
// *                                                                  *
// * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,  *
// * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES  *
// * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND         *
// * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT      *
// * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,     *
// * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING     *
// * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR    *
// * OTHER DEALINGS IN THE SOFTWARE.                                  *
// *                                                                  *
// ********************************************************************
//

#include "LindhardDamage.hh"

#include "TsVGeometryComponent.hh"
#include "TsVScorer.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"

#include <cmath> // To be able to take the power of numbers

LindhardDamage::LindhardDamage(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM, G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer) : TsVBinnedScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer)
{
	SetUnit("MeV");
}


LindhardDamage::~LindhardDamage() {;}


/////////////////////////////////////////////////////////////////////////////////////////////////////
// Function for calculating the Lindhard partition. A1/Z1 are for the projectile, A2/Z2 for the medium
G4double LindhardQ(G4double E, G4double A1, G4double Z1, G4double A2, G4double Z2) {
	if ((A1==0) || (Z1==0) || (A2==0) || (Z2==0)) {return 0;} // Returns a value of zero if the projectile or medium is not hadronic. This makes electrons and gammas do no damage which is not 100% correct.
	
	// The coefficients are from M. T. Robinson, “The dependence of radiation effects on the primary recoil energy,” in Proc. Int. Conf. Radiation-Induced Voids in Metal, Albany, NY, 1972, pp. 397–429. Values are retrieved from I. Jun, W. Kim and R. Evans, "Electron Nonionizing Energy Loss for Device Applications," in IEEE Transactions on Nuclear Science, vol. 56, no. 6, pp. 3229-3235, Dec. 2009, doi: 10.1109/TNS.2009.2033692.
	// Just as a side note, it really irritates me that these coefficients are so precisely stated, yet come with no uncertainty.
	G4double e =  E / (30.73547 * Z1 * Z2 * pow((pow(Z1,2./3.) + pow(Z2, 2./3.)),1./2.) * (A1 + A2) / A1);
	
	//G4double g = e + 0.40244 * pow(e, 3./4.) + 3.4008 * pow(e, 1./6.);
	G4double g = 0.74422 * e + 1.6812 * pow(e, 3./4.) + 0.90565 * pow(e, 1./6.); // A different parametrisation from A. Akkerman and J. Barak, "New Partition Factor Calculations for Evaluating the Damage of Low Energy Ions in Silicon," in IEEE Transactions on Nuclear Science, vol. 53, no. 6, pp. 3667-3674, Dec. 2006, doi: 10.1109/TNS.2006.884382.
	
	G4double k = 0.079524 * pow(Z1, 2./3.) * pow(Z2, 1./2.) * pow(A1 + A2, 3./2.) / (pow(pow(Z1, 2./3.) + pow(Z2, 2./3.), 3./4.) * pow(A1, 3./2.) * pow(A2, 1./2.));
	G4double Q = 1. / (1. + k * g);
	return Q;
}
/////////////////////////////////////////////////////////////////////////////////////////////////////

G4int runs = 0;


G4bool LindhardDamage::ProcessHits(G4Step* aStep,G4TouchableHistory*)
{	runs = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID(); // This updates the value for each run. I feel like there is a better way of doing it, just a single time at the end of the last run.
	if (!fIsActive) {
		fSkippedWhileInactive++;
		return false;
	}
	//G4double incidentEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
	G4String incidentName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
	G4String incidentProcess = "Primary";
	if (aStep->GetTrack()->GetParentID() != 0) { // Finds the process that created the incident particle. The if statement here is needed to avoid segmentation fault when accessing process name of a primary particle
		incidentProcess = aStep->GetTrack()->GetCreatorProcess()->GetProcessName();
	}
	
	// Looking at the secondary particles themselves
	const std::vector<const G4Track*>* secondaries = aStep->GetSecondaryInCurrentStep();
	
	// Looping through all the secondaries
	for (long unsigned int i = 0; i < secondaries->size(); i++) {

		const G4Track* track = (*secondaries)[i];

		G4String particleName = track->GetParticleDefinition()->GetParticleName();
			
		G4double energy = track->GetKineticEnergy();
		
		G4String creatorProcess = track->GetCreatorProcess()->GetProcessName();
		
		
					
					
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		if (incidentProcess == "Primary") {
		
			if (energy / eV > 21) {	// Filter out everything above 21eV (that's the threshold for Silicon displacement)
					
				G4Material* material = aStep->GetTrack()->GetMaterial();
				const G4ElementVector* elements = material->GetElementVector();
				const G4double* densities = material->GetVecNbOfAtomsPerVolume(); // The units of this is atoms/mm^3, it seems. A barn is e22 mm^2
					
				G4double effectiveQ = 0;
				// Looping through all the elements in the material, calculating Q * atom density / total atom density (to use it as a sort of weight pr element) for each of them
				for (long unsigned int n = 0; n < elements->size(); n++){
					const G4Element* element = (*elements)[n];
					effectiveQ += LindhardQ(energy / eV /*makes sure unit is eV*/, track->GetParticleDefinition()->GetAtomicMass(), track->GetParticleDefinition()->GetAtomicNumber(), element->GetAtomicMassAmu() , element->GetZ()) * densities[n] / material->GetTotNbOfAtomsPerVolume();
				}
					
			G4double totalnielvalue = effectiveQ * energy; 
			AccumulateHit(aStep, totalnielvalue);
			}

		}
		/////////////////////////////////////////////////////////////////////////////////////////////////////
		
		
		
		
		
		
		
		
			
	}
	return false;

}



void LindhardDamage::UserHookForEndOfRun()
{
	for (long unsigned int n=0; n<fFirstMomentMap.size(); n++){
		fFirstMomentMap[n] = fFirstMomentMap[n] / (runs+1);
	}
	G4cout << "Damage in 0.1mm silicon: " << 0.490966 /*Niel to D*/ * fFirstMomentMap[0] / keV / (0.0233 /*rho*L*/ ) << " /95MeVmb, runs= " << runs+1 << G4endl; // Gives out D in 0.1mm Si//(NIEL for 0.1mm silicon, in units of keV cm²/g)
}

