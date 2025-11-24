#ifndef EVENTACTION_HH
#define EVENTACTION_HH

#include "G4UserEventAction.hh"
#include "globals.hh"

class G4Event;
class RunAction;

class EventAction : public G4UserEventAction
{
public:
  explicit EventAction(RunAction* runAction);
  virtual ~EventAction();

  virtual void BeginOfEventAction(const G4Event* event) override;
  virtual void EndOfEventAction(const G4Event* event) override;

  // İstersen step’ten çağırıp event içi enerji biriktirebilirsin
  void AddScoringEdep(G4double edepMeV);

private:
  RunAction* fRunAction;
  G4double   fEdepEventMeV;
};

#endif
