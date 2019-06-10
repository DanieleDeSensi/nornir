/*
 * utils.cpp
 *
 * Created on: 09/07/2015
 *
 * =========================================================================
 *  Copyright (C) 2015-, Daniele De Sensi (d.desensi.software@gmail.com)
 *
 *  This file is part of nornir.
 *
 *  nornir is free software: you can redistribute it and/or
 *  modify it under the terms of the Lesser GNU General Public
 *  License as published by the Free Software Foundation, either
 *  version 3 of the License, or (at your option) any later version.

 *  nornir is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  Lesser GNU General Public License for more details.
 *
 *  You should have received a copy of the Lesser GNU General Public
 *  License along with nornir.
 *  If not, see <http://www.gnu.org/licenses/>.
 *
 * =========================================================================
 */

#include <nornir/utils.hpp>
#include <mammut/utils.hpp>

namespace nornir{
std::vector<std::string> getXdgConfigDirs(){
    char* confHome_c = getenv("XDG_CONFIG_DIRS");
    std::vector<std::string> confHomes;
    if(!confHome_c || strcmp(confHome_c, "") == 0){
        confHomes.push_back(std::string(XDG_CONFIG_DIR_FALLBACK));
    }else{
        confHomes = mammut::utils::split(std::string(confHome_c), ':');
    }
    for(std::string& s : confHomes){
        s += "/nornir/";
    }
    return confHomes;
}

std::string getRuntimeDir(bool userSpecific){
    char* runtimeDir_c = getenv("XDG_RUNTIME_DIR");
    if(!runtimeDir_c || strcmp(runtimeDir_c, "") == 0 || !userSpecific){
        runtimeDir_c = (char*) "/tmp/";
    }
    std::string runtimeDir = std::string(runtimeDir_c) + std::string("/nornir/");

    // Create dir if it does not exist
    if(!mammut::utils::existsDirectory(runtimeDir)){
        if(system((std::string("mkdir -p ") + runtimeDir).c_str())){
            throw std::runtime_error("Impossible to create nornir runtime dir.");
        }
        if(!userSpecific){
            if(system((std::string("chmod ugo+rwx ") + runtimeDir).c_str())){
                throw std::runtime_error("Impossible to set permissions on nornir runtime dir.");
            }
        }
    }
    return runtimeDir;
}
}