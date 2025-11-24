#include "EventAction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4SystemOfUnits.hh"
#include "G4ios.hh"

EventAction::EventAction(RunAction* runAction)
  : G4UserEventAction(),
    fRunAction(runAction),
    fEdepEventMeV(0.0)
{
}

EventAction::~EventAction() = default;

void EventAction::BeginOfEventAction(const G4Event*)
{
  fEdepEventMeV = 0.0;
}

void EventAction::EndOfEventAction(const G4Event* event)
{
  if (!fRunAction) return;

  // Buraya istersen event bazlı istatistik aktarabilirsin.
  // Şimdilik sadece saklıyoruz.
  /*
  if (fEdepEventMeV > 0.0) {
    G4cout << "Event " << event->GetEventID()
           << " : Edep = " << fEdepEventMeV << " MeV" << G4endl;
  }
  */
}

void EventAction::AddScoringEdep(G4double edepMeV)
{
  fEdepEventMeV += edepMeV;
}
