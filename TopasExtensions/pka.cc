// Scorer for pka
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

#include "pka.hh"

#include "TsVGeometryComponent.hh"
#include "TsVScorer.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UIcommand.hh"
#include "G4Tokenizer.hh"
#include "G4SystemOfUnits.hh"
#include "G4PSDirectionFlag.hh"

pka::pka(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
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
}

pka::~pka()
{;}	

G4bool pka::ProcessHits(G4Step* aStep, G4TouchableHistory*)
{
    if (!fIsActive) {
		return false;
	}
    if (aStep->GetTrack()->GetCurrentStepNumber() == 1) {
    		// Fill in the energy
		G4StepPoint* theStepPoint = aStep->GetPreStepPoint();
		fEnergy	     = theStepPoint->GetKineticEnergy();
		
		// Fill in the creation type
		const G4VProcess* creatorProcess = aStep->GetTrack()->GetCreatorProcess();
		if (creatorProcess){
			fCreatorProcess = creatorProcess->GetProcessName();}
		else{
			fCreatorProcess = "Primary";}
			
		// Fill in the type of particle created
		fParticleName = aStep->GetTrack()->GetParticleDefinition()->GetParticleName();
		
		// Get the ID of the run
		fRunID = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();//GetRunID();
		
		// Get the name of the incident particle which created this new particle
		fIncidentParticle = GetIncidentParticleDefinition()->GetParticleName();
		
		// Get the atomic number of the created particle
		Z = aStep->GetTrack()->GetParticleDefinition()->GetAtomicNumber();
		
		// Get the atomic mass of the created particle
		A = aStep->GetTrack()->GetParticleDefinition()->GetAtomicMass();
		
        fNtuple->Fill();
        return true;
	}
    return false;
}

void pka::AccumulateEvent()
{
    TsVNtupleScorer::AccumulateEvent();
}

//void TsScorePhaseSpace::UpdateForEndOfRun() {fEnergy	        = 0.; TsVScorer::UpdateForEndOfRun();}

void pka::Output()
{
	std::ostringstream title;
	title << "TOPAS Primary Knocked Atoms scorer" << G4endl;
	title << "Number of Original Histories:" << G4endl << GetScoredHistories() << G4endl;
	fNtuple->fHeaderPrefix = title.str();
	fNtuple->Write();
	
	G4cout << "PKA-scorer has now written output!" << G4endl;
}
