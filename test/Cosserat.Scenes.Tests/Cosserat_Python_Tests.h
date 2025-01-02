/******************************************************************************
*                              Cosserat plugin                                *
*                  (c) 2024 CNRS, University of Lille, INRIA                  *
*                                                                             *
* This program is free software; you can redistribute it and/or modify it     *
* under the terms of the GNU Lesser General Public License as published by    *
* the Free Software Foundation; either version 2.1 of the License, or (at     *
* your option) any later version.                                             *
*                                                                             *
* This program is distributed in the hope that it will be useful, but WITHOUT *
* ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or       *
* FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License *
* for more details.                                                           *
*                                                                             *
* You should have received a copy of the GNU Lesser General Public License    *
* along with this program. If not, see <http://www.gnu.org/licenses/>.        *
*******************************************************************************
* Contact information: contact@sofa-framework.org                             *
******************************************************************************/
#include <SofaPython3Testing/PythonTest.h>
#include <sofa/helper/Utils.h>
#include <sofa/helper/system/FileSystem.h>
#include <SofaPython3Testing/PythonTestExtractor.h>
#include "Cosserat_Python_Tests.h"

using sofapython3::PythonTest ;
using sofapython3::PythonTestExtractor ;
using sofapython3::PrintTo ;
using std::string;

#include <sofa/core/logging/PerComponentLoggingMessageHandler.h>
#include <sofa/helper/logging/MessageDispatcher.h>
using sofa::helper::logging::MessageDispatcher;
using sofa::helper::logging::MainPerComponentLoggingMessageHandler;

namespace
{

bool init()
{
    MessageDispatcher::addHandler(&MainPerComponentLoggingMessageHandler::getInstance()) ;
    return true;
}

static int _inited_=init();

class Cosserat : public PythonTest {};

/// static build of the test list
static struct Cosserat_Scene_tests : public PythonTestExtractor
{
    Cosserat_Scene_tests()
    {
        std::string executable_directory = sofa::helper::Utils::getExecutableDirectory();
        if(sofa::helper::system::FileSystem::exists("./Cosserat.Scene.Tests/Components")){
            addTestDirectory("./Cosserat.Scene.Tests/Components", "Cosserat_Scene_");
        }else{
            addTestDirectory(executable_directory+"/Cosserat.Scene.Tests.d/Components", "Cosserat_Scene_");
        }
    }
} python_tests;

/// run test list using the custom name function getTestName.
/// this allows to do gtest_filter=*FileName*
INSTANTIATE_TEST_SUITE_P(Cosserat,
                        Cosserat,
                        ::testing::ValuesIn(python_tests.extract()),
                        Cosserat::getTestName);

TEST_P(Cosserat, all_tests) { run(GetParam()); }

}
