/******************************************************************************
 *       SOFA, Simulation Open-Framework Architecture, version 1.0 beta 4      *
 *                (c) 2006-2009 MGH, INRIA, USTL, UJF, CNRS                    *
 *                                                                             *
 * This program is free software; you can redistribute it and/or modify it     *
 * under the terms of the GNU General Public License as published by the Free  *
 * Software Foundation; either version 2 of the License, or (at your option)   *
 * any later version.                                                          *
 *                                                                             *
 * This program is distributed in the hope that it will be useful, but WITHOUT *
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
 * FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for    *
 * more details.                                                               *
 *                                                                             *
 * You should have received a copy of the GNU General Public License along     *
 * with this program; if not, write to the Free Software Foundation, Inc., 51  *
 * Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.                   *
 *******************************************************************************
 *                            SOFA :: Applications                             *
 *                                                                             *
 * Authors: M. Adam, J. Allard, B. Andre, P-J. Bensoussan, S. Cotin, C. Duriez,*
 * H. Delingette, F. Falipou, F. Faure, S. Fonteneau, L. Heigeas, C. Mendoza,  *
 * M. Nesme, P. Neumann, J-P. de la Plata Alcade, F. Poyer and F. Roy          *
 *                                                                             *
 * Contact information: contact@sofa-framework.org                             *
 ******************************************************************************/
#include <Cosserat/config.h>

#include <sofa/core/ObjectFactory.h>
#include <sofa/helper/system/PluginManager.h>
using sofa::helper::system::PluginManager;

#include <cstring>
#include <string>

namespace Cosserat {


#ifdef COSSERAT_USES_SOFTROBOTS
extern void registerQPSlidingConstraint(sofa::core::ObjectFactory* factory);
extern void registerCosseratActuatorConstraint(sofa::core::ObjectFactory* factory);
#endif

extern void registerCosseratNeedleSlidingConstraint(sofa::core::ObjectFactory* factory);
extern void registerCosseratSlidingConstraint(sofa::core::ObjectFactory* factory);
extern void registerPointsManager(sofa::core::ObjectFactory* factory);
extern void registerProjectionEngine(sofa::core::ObjectFactory* factory);
extern void registerBeamHookeLawForceField(sofa::core::ObjectFactory* factory);
extern void registerBeamHookeLawForceFieldRigid(sofa::core::ObjectFactory* factory);
extern void registerCosseratInternalActuation(sofa::core::ObjectFactory* factory);
extern void registerDifferenceMultiMapping(sofa::core::ObjectFactory* factory);
extern void registerDiscreteCosseratMapping(sofa::core::ObjectFactory* factory);
extern void registerDiscretDynamicCosseratMapping(sofa::core::ObjectFactory* factory);
extern void registerLegendrePolynomialsMapping(sofa::core::ObjectFactory* factory);
extern void registerRigidDistanceMapping(sofa::core::ObjectFactory* factory);

extern "C" {
SOFA_COSSERAT_API void initExternalModule();
SOFA_COSSERAT_API const char *getModuleLicense();
SOFA_COSSERAT_API const char *getModuleName();
SOFA_COSSERAT_API const char *getModuleVersion();
SOFA_COSSERAT_API const char *getModuleDescription();
SOFA_COSSERAT_API const char *getModuleComponentList();
SOFA_COSSERAT_API void registerObjects(sofa::core::ObjectFactory* factory);
}

// Here are just several convenient functions to help user to know what contains
// the plugin

void initExternalModule() {
  static bool first = true;
  if (first)
  {
    sofa::helper::system::PluginManager::getInstance().registerPlugin(MODULE_NAME);
    first = false;
  }
  // Automatically load the STLIB plugin if available.
  if (!PluginManager::getInstance().findPlugin("STLIB").empty())
  {
    PluginManager::getInstance().loadPlugin("STLIB");
  }

#ifdef SOFTROBOTS_PYTHON
  PythonEnvironment::addPythonModulePathsForPluginsByName("Cosserat");
#endif
}

void registerObjects(sofa::core::ObjectFactory* factory)
{

#ifdef COSSERAT_USES_SOFTROBOTS
  registerQPSlidingConstraint(factory);
  registerCosseratActuatorConstraint(factory);
#endif
  registerCosseratNeedleSlidingConstraint(factory);
  registerCosseratSlidingConstraint(factory);
  registerPointsManager(factory);
  registerProjectionEngine(factory);
  registerBeamHookeLawForceField(factory);
  registerBeamHookeLawForceFieldRigid(factory);
  registerCosseratInternalActuation(factory);
  registerDifferenceMultiMapping(factory);
  registerDiscreteCosseratMapping(factory);
  registerDiscretDynamicCosseratMapping(factory);
  registerLegendrePolynomialsMapping(factory);
  registerRigidDistanceMapping(factory);
}

const char *getModuleLicense() { return "LGPL"; }

const char *getModuleName() { return Cosserat::MODULE_NAME; }

const char *getModuleVersion() { return Cosserat::MODULE_VERSION; }

const char *getModuleDescription() {
  return "This plugin is used to implement slender object";
}

const char *getModuleComponentList() {
  // string containing the names of the classes provided by the plugin
  static std::string classes =
      sofa::core::ObjectFactory::getInstance()->listClassesFromTarget(
          sofa_tostring(SOFA_TARGET));
  return classes.c_str();
}

} // namespace cosserat
