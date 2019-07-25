//# fixlsmeta.cc: Program to fix the meta info of the LofarStMan
//# Copyright (C) 2009
//# ASTRON (Netherlands Institute for Radio Astronomy)
//# P.O.Box 2, 7990 AA Dwingeloo, The Netherlands
//#
//# This file is part of the LOFAR software suite.
//# The LOFAR software suite is free software: you can redistribute it and/or
//# modify it under the terms of the GNU General Public License as published
//# by the Free Software Foundation, either version 3 of the License, or
//# (at your option) any later version.
//#
//# The LOFAR software suite is distributed in the hope that it will be useful,
//# but WITHOUT ANY WARRANTY; without even the implied warranty of
//# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//# GNU General Public License for more details.
//#
//# You should have received a copy of the GNU General Public License along
//# with the LOFAR software suite. If not, see <http://www.gnu.org/licenses/>.
//#
//# $Id$

#include <lofar_config.h>
#include <LofarStMan/LofarStMan.h>
#include <Common/LofarLogger.h>
#include <Common/Exception.h>
#include <casacore/casa/IO/AipsIO.h>
#include <casacore/casa/Containers/BlockIO.h>
#include <casacore/casa/Quanta/MVTime.h>
#include <casacore/casa/Quanta/Quantum.h>
#include <casacore/casa/Inputs/Input.h>
#include <casacore/casa/iostream.h>

using namespace casacore;
using namespace LOFAR;

// Use a terminate handler that can produce a backtrace.
Exception::TerminateHandler t(Exception::terminate);

int main (int argc, char* argv[])
{
  try {
    Input inputs(1);
    // Define the input keywords
    inputs.version("201020515GvD");
    inputs.create ("ms", "",
		   "Name of input MS",
		   "string");
    inputs.create ("timedivisor", "1",
		   "Factor to divide the time interval by",
		   "float");
    inputs.create ("timemultiplier", "1",
		   "Factor to multiply the time interval with",
		   "float");
    inputs.create ("starttime", "",
                   "New start time like 2-Mar-2012/07:34:12.52"
                   "string");
    // Fill the input value from the command line.
    inputs.readArguments (argc, argv);

    // Get and check the input specification.
    String msin (inputs.getString("ms"));
    if (msin == "") {
      throw AipsError(" an input MS must be given");
    }
    double timeDivisor (inputs.getDouble("timedivisor"));
    double timeMultiplier (inputs.getDouble("timemultiplier"));
    string startStr (inputs.getString("starttime"));

    Bool asBigEndian;
    uInt alignment;
    Block<int32> ant1;
    Block<int32> ant2;
    double startTime;
    double timeIntv;
    uint32 nChan;
    uint32 nPol;
    uint32 nBytesPerNrValidSamples=2;
    double maxNrSample;
    int version;
    {
      // Open and read the meta file.
      AipsIO aio(msin + "/table.f0meta");
      version = aio.getstart ("LofarStMan");
      ASSERTSTR (version <= 3,
                 "fixlsmeta can only handle up to version 3");
      aio >> ant1 >> ant2 >> startTime >> timeIntv >> nChan
          >> nPol >> maxNrSample >> alignment >> asBigEndian;
      if (version > 1) {
        aio >> nBytesPerNrValidSamples;
      }
      aio.getend();
    }
    // Fix the interval.
    if (timeDivisor != 1  &&  timeMultiplier != 1) {
      cout << "time interval changed from " << timeIntv << " to ";
      timeIntv /= timeDivisor;
      timeIntv *= timeMultiplier;
      cout <<timeIntv << endl;
    }
    // Fix the start time.
    if (!startStr.empty()) {
      Quantity q;
      ASSERTSTR (MVTime::read(q, startStr, True),
                 startStr << " is an invalid date/time");
      cout << "start time changed from " << startTime << " to ";
      startTime = q.getValue("s");
      cout <<startTime << endl;
    }
    {
      // Create and write the meta file.
      AipsIO aio(msin + "/table.f0meta", ByteIO::New);
      aio.putstart ("LofarStMan", version);
      aio << ant1 << ant2 << startTime << timeIntv << nChan
          << nPol << maxNrSample << alignment << asBigEndian;
      if (version > 1) {
        aio << nBytesPerNrValidSamples;
      }
      aio.putend();
    }
 
  } catch (Exception& ex) {
    cerr << "Caught a LOFAR exception: " << ex << endl;
    return 1;
  } catch (AipsError& x) {
    cerr << "Caught an AIPS error: " << x.getMesg() << endl;
    return 1;
  } 
  return 0;                           // exit with success status
}
