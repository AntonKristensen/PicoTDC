// Scorer for ppka
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

#include "ppka.hh"

#include "TsVGeometryComponent.hh"
#include "TsVScorer.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"

ppka::ppka(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
                                     G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer)
    : TsVNtupleScorer(pM, mM, gM, scM, eM, scorerName, quantity, outFileName, isSubScorer),
      fEnergy(0.), fCreatorProcess(""), fParticleName(""), fRunID(0)
{
    fNtuple->RegisterColumnF(&fEnergy, "Energy", "MeV");
    fNtuple->RegisterColumnS(&fCreatorProcess, "Creator Process Name");
    fNtuple->RegisterColumnS(&fParticleName, "Particle Name");
    fNtuple->RegisterColumnI(&fRunID, "RunID");
    fNtuple->RegisterColumnS(&fIncidentParticle, "Incident Particle Name");
    fNtuple->RegisterColumnI(&Z, "Atomic number");
    fNtuple->RegisterColumnI(&A, "Atomic mass");
    fNtuple->RegisterColumnI(&fCreatorSecondary, "Created by secondary");
    fNtuple->RegisterColumnF(&fIncidentEnergy, "Incident Particle Energy", "MeV");
}

ppka::~ppka()
{;}

G4bool ppka::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
	if (!fIsActive) {
		return false;
	}
	
	// Checks if there are any secondary particles made, and if there aren't any then stop and move on
	G4int nSecondaryParticles = aStep->GetNumberOfSecondariesInCurrentStep();
	if (nSecondaryParticles == 0) { 
		return false; 
	}
	
	// Writing out how many particles got made
	//G4cout << nSecondaryParticles << G4endl;
	
	// Looking at the secondary particles themselves
	const std::vector<const G4Track*>* secondaries
			= aStep->GetSecondaryInCurrentStep();
			
	// Looping through all the secondaries
	for (long unsigned int i = 0; i < secondaries->size(); i++) {

		const G4Track* track = (*secondaries)[i];

		G4String particleName = track->GetParticleDefinition()->GetParticleName();
			
		G4double energy = track->GetKineticEnergy();
			
		G4String creatorProcess = track->GetCreatorProcess()->GetProcessName();

		//G4cout << particleName << " " << energy << "MeV, by a " << aStep->GetTrack()->GetParticleDefinition()->GetParticleName() << " via " << creatorProcess <<G4endl;
		
		// Filling in variables that get written
		fEnergy	     = energy;
		fCreatorProcess = creatorProcess;
		fParticleName = particleName;
		fRunID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
		fIncidentParticle = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
		Z = track->GetParticleDefinition()->GetAtomicNumber();
		A = track->GetParticleDefinition()->GetAtomicMass();
		fCreatorSecondary = aStep->GetTrack()->GetParentID();
		fIncidentEnergy = aStep->GetPreStepPoint()->GetKineticEnergy();
		fNtuple->Fill(); // Fill into the tuple
		}
	
	return true;
	
	//return false;
}

void ppka::AccumulateEvent()
{
    TsVNtupleScorer::AccumulateEvent();
}

//void TsScorePhaseSpace::UpdateForEndOfRun() {fEnergy	        = 0.; TsVScorer::UpdateForEndOfRun();}

void ppka::Output()
{
	std::ostringstream title;
	title << "TOPAS Primary Knocked Atoms scorer" << G4endl;
	title << "Number of Original Histories:" << G4endl << GetScoredHistories() << G4endl;
	fNtuple->fHeaderPrefix = title.str();
	fNtuple->Write();
	
	G4cout << "PPKA-scorer has now written output!" << G4endl;
}
