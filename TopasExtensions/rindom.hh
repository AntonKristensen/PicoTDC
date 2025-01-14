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

#ifndef rindom_hh
#define rindom_hh

#include "TsVNtupleScorer.hh"

#include <stdint.h>

class rindom : public TsVNtupleScorer
{
public:
	rindom(TsParameterManager* pM, TsMaterialManager* mM, TsGeometryManager* gM, TsScoringManager* scM, TsExtensionManager* eM,
					  G4String scorerName, G4String quantity, G4String outFileName, G4bool isSubScorer);
	~rindom();

	virtual G4bool ProcessHits(G4Step*, G4TouchableHistory*);

	void AccumulateEvent();


protected:
	void Output();
	//void Clear();

	G4float fEnergy;
	G4String fCreatorProcess;
	G4String fParticleName;
	G4int fRunID;
	G4String fIncidentParticle;
	G4int Z;
	G4int A;
	G4int fCreatorSecondary;
	G4float fIncidentEnergy;
	G4float fTOPASTime;
	G4int fNumberOfSecondaries;
	//G4bool fLastStep;
	G4float fExitEnergy;


};

#endif
