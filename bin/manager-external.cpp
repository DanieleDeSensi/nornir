/*
 * manager-external.cpp
 *
 * Created on: 21/06/2016
 *
 * This executable starts a manager which can monitor applications not written
 * with the Nornir framework. Albeit more of such applications can run concurrently,
 * no coordination between them is provided. The result may thus be inconsistent.
 * Please use this binary to control at most one application at a time.
 *
 * Since the manager needs to access some architectures priviledged operations
 * (e.g. changing frequency, reading energy, etc...), it is possible that
 * this process could be run by using 'sudo'.
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

#include <nornir/nornir.hpp>

int main(int argc, char * argv[]){
    nornir::InstrumenterServer is(true);
    is.start();
    is.join();
    return 0;
}

