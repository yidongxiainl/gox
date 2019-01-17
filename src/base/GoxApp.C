#include "GoxApp.h"
#include "Moose.h"
#include "AppFactory.h"
//#include "ModulesApp.h"
#include "MooseSyntax.h"

//module include
#include "ContactApp.h"
#include "HeatConductionApp.h"
#include "MiscApp.h"
#include "SolidMechanicsApp.h"
//#include "PhaseFieldApp.h"
#include "TensorMechanicsApp.h"
//#include "XFEMApp.h"

//local input for material properties
#include "PorousMediaBase.h"
#include "GoxHeatTimeDerivative.h"
#include "GasMixtureDiffusion.h"
#include "GasMixtureTimeDerivative.h"
#include "GoxHeatConduction.h"
#include "GoxHeatReaction.h"
#include "ReactionSourceSinkKernel.h"

#include "BulkDensityAux.h"


template<>
InputParameters validParams<GoxApp>()
{
  InputParameters params = validParams<MooseApp>();

  params.set<bool>("use_legacy_uo_initialization") = false;
  params.set<bool>("use_legacy_uo_aux_computation") = false;
  params.set<bool>("use_legacy_output_syntax") = false;

  return params;
}

GoxApp::GoxApp(InputParameters parameters) :
    MooseApp(parameters)
{
  Moose::registerObjects(_factory);
//  ModulesApp::registerObjects(_factory);
  ContactApp::registerObjects(_factory);
  HeatConductionApp::registerObjects(_factory);
  MiscApp::registerObjects(_factory);
  SolidMechanicsApp::registerObjects(_factory);
//  PhaseFieldApp::registerObjects(_factory);
//  XFEMApp::registerObjects(_factory);
  TensorMechanicsApp::registerObjects(_factory);

  GoxApp::registerObjects(_factory);

  Moose::associateSyntax(_syntax, _action_factory);
//  ModulesApp::associateSyntax(_syntax, _action_factory);
  SolidMechanicsApp::associateSyntax(_syntax, _action_factory);
  ContactApp::associateSyntax(_syntax, _action_factory);
  HeatConductionApp::associateSyntax(_syntax, _action_factory);
  MiscApp::associateSyntax(_syntax, _action_factory);
  TensorMechanicsApp::associateSyntax(_syntax, _action_factory);
//  XFEMApp::associateSyntax(_syntax, _action_factory);

  GoxApp::associateSyntax(_syntax, _action_factory);
}

GoxApp::~GoxApp()
{
}

// External entry point for dynamic application loading
extern "C" void GoxApp__registerApps() { GoxApp::registerApps(); }
void
GoxApp::registerApps()
{
  registerApp(GoxApp);
}

// External entry point for dynamic object registration
extern "C" void GoxApp__registerObjects(Factory & factory) { GoxApp::registerObjects(factory); }
void
GoxApp::registerObjects(Factory & factory)
{
  // Register our new material class so we can use it.
  registerMaterial(PorousMediaBase);
  registerKernel(GoxHeatTimeDerivative);
  registerKernel(GasMixtureDiffusion);
  registerKernel(GasMixtureTimeDerivative);
  registerKernel(GoxHeatConduction);
  registerKernel(GoxHeatReaction);
  registerKernel(ReactionSourceSinkKernel);

  registerAux(BulkDensityAux);
}

// External entry point for dynamic syntax association
extern "C" void GoxApp__associateSyntax(Syntax & syntax, ActionFactory & action_factory) { GoxApp::associateSyntax(syntax, action_factory); }
void
GoxApp::associateSyntax(Syntax & /*syntax*/, ActionFactory & /*action_factory*/)
{
}
