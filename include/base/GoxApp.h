#ifndef GOXAPP_H
#define GOXAPP_H

#include "MooseApp.h"

class GoxApp;

template<>
InputParameters validParams<GoxApp>();

class GoxApp : public MooseApp
{
public:
  GoxApp(InputParameters parameters);
  virtual ~GoxApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* GOXAPP_H */
