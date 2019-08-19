/*
 * nornir_blackbox.cpp
 *
 * Monitors and adapt an external application by using the riff library.
 * Needs the explicit channel name. This is deprecated. Use manager-external.
 *
 * Created on: 21/06/2016
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

#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <unistd.h>
#include <sys/types.h>
#include <sys/wait.h>
#include <nornir/nornir.hpp>

using namespace nornir;
using namespace std;

static pair<std::string, char** const> getArgs(char* cmd){
    vector<string> strings;
    istringstream f(cmd);
    string s;    
    while (getline(f, s, ' ')) {
        strings.push_back(s);
    }
    char** args = (char**) malloc(sizeof(char*) * (strings.size() + 1));
    for(size_t i = 0; i < strings.size(); i++){
        args[i] = (char*) malloc(sizeof(char) * (strings[i].length() + 1));
        strcpy(args[i], strings[i].c_str());
    }
    args[strings.size()] = NULL;
    return pair<std::string, char** const>(strings[0], args);
}

int main(int argc, char * argv[]) {
    char* command = NULL;
    if(argc != 2) {
        cerr << "use: " 
             << argv[0] 
             << " command" << endl;
        return -1;
    }   
    command = argv[1];
    pid_t pid = fork();

    if(pid){
        Parameters p("parameters.xml");
        ManagerBlackBox m(pid, p);
        m.start();
        waitpid(pid, NULL, 0);
        m.join();
    }else{
        sleep(1);
        pair<string, char** const> subargs = getArgs(command);
        execv(subargs.first.c_str(), subargs.second);       
    }

    return 0;
}

