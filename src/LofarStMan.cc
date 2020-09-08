//# LofarStMan.cc: Storage Manager for the main table of a LOFAR MS
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

#include <LofarStMan/LofarStMan.h>
#include <LofarStMan/LofarColumn.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/DataMan/DataManError.h>
#include <casacore/casa/Containers/Record.h>
#include <casacore/casa/Containers/BlockIO.h>
#include <casacore/casa/IO/AipsIO.h>
#include <casacore/casa/OS/CanonicalConversion.h>
#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/OS/DOos.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/iostream.h>

using namespace casacore;


namespace LOFAR {

LofarStMan::LofarStMan (const String& dataManName)
: DataManager    (),
  itsDataManName (dataManName),
  itsFD          (-1),
  itsRegFile     (0),
  itsSeqFile     (0)
{}

LofarStMan::LofarStMan (const String& dataManName,
                        const Record&)
: DataManager    (),
  itsDataManName (dataManName),
  itsFD          (-1),
  itsRegFile     (0),
  itsSeqFile     (0)
{}

LofarStMan::LofarStMan (const LofarStMan& that)
: DataManager    (),
  itsDataManName (that.itsDataManName),
  itsFD          (-1),
  itsRegFile     (0),
  itsSeqFile     (0)
{}

LofarStMan::~LofarStMan()
{
  for (uInt i=0; i<ncolumn(); i++) {
    delete itsColumns[i];
  }
  closeFiles();
}

DataManager* LofarStMan::clone() const
{
  return new LofarStMan (*this);
}

String LofarStMan::dataManagerType() const
{
  return "LofarStMan";
}

String LofarStMan::dataManagerName() const
{
  return itsDataManName;
}

Record LofarStMan::dataManagerSpec() const
{
  return itsSpec;
}


DataManagerColumn* LofarStMan::makeScalarColumn (const String& name,
                                                 int dtype,
                                                 const String&)
{
  LofarColumn* col;
  if (name == "TIME"  ||  name == "TIME_CENTROID") {
    col = new TimeColumn(this, dtype);
  } else if (name == "ANTENNA1") {
    col = new Ant1Column(this, dtype);
  } else if (name == "ANTENNA2") {
    col = new Ant2Column(this, dtype);
  } else if (name == "INTERVAL"  ||  name == "EXPOSURE") {
    col = new IntervalColumn(this, dtype);
  } else if (name == "FLAG_ROW") {
    col = new FalseColumn(this, dtype);
  } else {
    col = new ZeroColumn(this, dtype);
  }
  itsColumns.push_back (col);
  return col;
}

DataManagerColumn* LofarStMan::makeDirArrColumn (const String& name,
                                                 int dataType,
                                                 const String& dataTypeId)
{
  return makeIndArrColumn (name, dataType, dataTypeId);
}

DataManagerColumn* LofarStMan::makeIndArrColumn (const String& name,
                                                 int dtype,
                                                 const String&)
{
  LofarColumn* col;
  if (name == "UVW") {
    col = new UvwColumn(this, dtype);
  } else if (name == "DATA") {
    col = new DataColumn(this, dtype);
  } else if (name == "FLAG") {
    col = new FlagColumn(this, dtype);
  } else if (name == "FLAG_CATEGORY") {
    col = new FlagCatColumn(this, dtype);
  } else if (name == "WEIGHT") {
    col = new WeightColumn(this, dtype);
  } else if (name == "SIGMA") {
    col = new SigmaColumn(this, dtype);
  } else if (name == "WEIGHT_SPECTRUM") {
    col = new WSpectrumColumn(this, dtype);
  } else {
    throw DataManError (name + " is unknown column for LofarStMan");
  }
  itsColumns.push_back (col);
  return col;
}

DataManager* LofarStMan::makeObject (const String& group, const Record& spec)
{
  // This function is called when reading a table back.
  return new LofarStMan (group, spec);
}

void LofarStMan::registerClass()
{
  DataManager::registerCtor ("LofarStMan", makeObject);
}

Bool LofarStMan::isRegular() const
{
  return False;
}
Bool LofarStMan::canAddRow() const
{
  return False;
}
Bool LofarStMan::canRemoveRow() const
{
  return False;
}
Bool LofarStMan::canAddColumn() const
{
  return True;
}
Bool LofarStMan::canRemoveColumn() const
{
  return True;
}

void LofarStMan::addRow (uInt)
{
  throw DataManError ("LofarStMan cannot add rows");
}
void LofarStMan::removeRow (uInt)
{
  throw DataManError ("LofarStMan cannot remove rows");
}
void LofarStMan::addColumn (DataManagerColumn*)
{}
void LofarStMan::removeColumn (DataManagerColumn*)
{}

Bool LofarStMan::flush (AipsIO&, Bool)
{
  return False;
}

void LofarStMan::create (uInt nrows)
{
  itsNrRows = nrows;
}

void LofarStMan::open (uInt, AipsIO&)
{
  throw DataManError ("LofarStMan::open should never be called");
}
uInt LofarStMan::open1 (uInt, AipsIO&)
{
  // Read meta info.
  init();
  openFiles (table().isWritable());
  return itsNrRows;
}

void LofarStMan::prepare()
{
  for (uInt i=0; i<ncolumn(); i++) {
    itsColumns[i]->prepareCol();
  }
}

void LofarStMan::openFiles (bool writable)
{
  // Open the data file using unbuffered IO.
  // First close if needed.
  closeFiles();
  String fname (fileName() + "data");
  itsFD = FiledesIO::open (fname.c_str(), writable);
  itsRegFile = new FiledesIO (itsFD);
  // Set correct number of rows.
  itsNrRows = itsRegFile->length() / itsBlockSize * itsAnt1.size();
  // Map the file with seqnrs.
  mapSeqFile();
  // Size the buffer if needed.
  if ((long long) (itsBuffer.size()) < itsBLDataSize) {
    itsBuffer.resize (itsBLDataSize);
  }
  itsSpec.define ("useSeqnrFile", itsSeqFile!=0);
}

void LofarStMan::mapSeqFile()
{
  delete itsSeqFile;
  itsSeqFile = 0;
  try {
    itsSeqFile = new MMapIO (fileName() + "seqnr");
  } catch (...) {
    delete itsSeqFile; 
    itsSeqFile = 0;
  }
  // Check the size of the sequencenumber file, close file if it doesn't match.
  // It should contain the nr of time slots.
  if (itsSeqFile && (itsSeqFile->getFileSize() !=
		     Int64(itsNrRows / itsAnt1.size() * sizeof(uInt)))) {
    delete itsSeqFile;
    itsSeqFile = 0;
  }
}

void LofarStMan::closeFiles()
{
  if (itsFD >= 0) {
    FiledesIO::close (itsFD);
    itsFD = -1;
  }
  delete itsRegFile;
  itsRegFile = 0;
  delete itsSeqFile;
  itsSeqFile = 0;
}

void LofarStMan::resync (uInt)
{
  throw DataManError ("LofarStMan::resync should never be called");
}
uInt LofarStMan::resync1 (uInt)
{
  uInt nrows = itsRegFile->length() / itsBlockSize * itsAnt1.size();
  // Reopen file if different nr of rows.
  if (nrows != itsNrRows) {
    openFiles (table().isWritable());
  }
  return itsNrRows;
}

void LofarStMan::reopenRW()
{
  openFiles (true);
}

void LofarStMan::deleteManager()
{
  closeFiles();
  DOos::remove (fileName()+"meta", False, False);
  DOos::remove (fileName()+"data", False, False);
  DOos::remove (fileName()+"seqnr", False, False);
}

void LofarStMan::init()
{
  AipsIO aio(fileName() + "meta");
  itsVersion = aio.getstart ("LofarStMan");
  if (itsVersion > 3) {
    throw DataManError ("LofarStMan can only handle up to version 3");
  }
  Bool asBigEndian;
  uInt alignment;
  if (itsVersion == 2) {
    // In version 2 antenna1 and antenna2 were swapped.
    aio >> itsAnt2 >> itsAnt1;
  } else {
    aio >> itsAnt1 >> itsAnt2;
  }
  aio >> itsStartTime >> itsTimeIntv >> itsNChan
      >> itsNPol >> itsMaxNrSample >> alignment >> asBigEndian;
  if (itsVersion > 1) {
    aio >> itsNrBytesPerNrValidSamples;
  } else {
    itsNrBytesPerNrValidSamples = 2;
  }
  aio.getend();
  // Set start time to middle of first time slot.
  itsStartTime += itsTimeIntv*0.5;
  AlwaysAssert (itsAnt1.size() == itsAnt2.size(), AipsError);
  uInt nrbl = itsAnt1.size();
  itsDoSwap  = (asBigEndian != HostInfo::bigEndian());
  // A block contains a possibly uInt magic value, uInt seqnr, Complex data per
  // baseline,chan,pol and nsample per baseline,chan. Align it as needed.
  itsBLDataSize = itsNChan * itsNPol * 8;    // #bytes/baseline
  if (alignment <= 1) {
    itsDataStart = itsVersion==1 ? 4:8;
    itsSampStart = itsDataStart + nrbl*itsBLDataSize;
    itsBlockSize = itsSampStart + nrbl*itsNChan*itsNrBytesPerNrValidSamples;
  } else {
    itsDataStart = alignment;
    itsSampStart = itsDataStart + (nrbl*itsBLDataSize + alignment-1)
      / alignment * alignment;
    itsBlockSize = itsSampStart + (nrbl*itsNChan*itsNrBytesPerNrValidSamples +
                                   alignment-1) / alignment * alignment;
  }
  switch (itsNrBytesPerNrValidSamples) {
  case 1:
    itsNSampleBuf1.resize (itsNChan * 2);
    break;
  case 2:
    itsNSampleBuf2.resize (itsNChan * 2);
    break;
  case 4:
    itsNSampleBuf4.resize (itsNChan * 2);
    break;
  default:
    throw DataManError ("LofarStMan invalid nrbytesPerNrValidSamples");
  }
  // Fill the specification record (only used for reporting purposes).
  itsSpec.define ("version", itsVersion);
  itsSpec.define ("alignment", alignment);
  itsSpec.define ("bigEndian", asBigEndian);
  itsSpec.define ("maxNrSample", itsMaxNrSample);
  itsSpec.define ("nrBytesPerNrValidSamples", itsNrBytesPerNrValidSamples);
  itsSpec.define ("startTime", itsStartTime);
  itsSpec.define ("timeInterval", itsTimeIntv);
  itsSpec.define ("nbaseline", Int(itsAnt1.size()));
}

Double LofarStMan::time (uInt blocknr)
{
  uInt seqnr;
  const void* ptr;
  if (itsSeqFile) {
    ptr = itsSeqFile->getReadPointer(blocknr * sizeof(uInt));
  } else {
    ptr = getReadPointer (blocknr, 0, sizeof(uInt));
    // Version 2 and later have a magic value before the seqnr.
    if (itsVersion >= 2) {
      if (itsDoSwap) {
	CanonicalConversion::reverse4 (&seqnr, ptr);
      } else {
	seqnr = *static_cast<const uInt*>(ptr);
      }
      if (seqnr != 0x00000da7a) {
	throw DataManError ("Magic number mismatch in block " +
			    String::toString(blocknr) +
			    " of LofarStMan data file " + itsRegFile->fileName());
      }
      ptr = getReadPointer (blocknr, sizeof(uInt), sizeof(uInt));
    }
  }
  if (itsDoSwap) {
    CanonicalConversion::reverse4 (&seqnr, ptr);
  } else {
    seqnr = *static_cast<const uInt*>(ptr);
  }
  return itsStartTime + seqnr*itsTimeIntv;
}

void LofarStMan::getData (uInt rownr, Complex* buf)
{
  uInt blocknr = rownr / itsAnt1.size();
  uInt baseline = rownr - blocknr*itsAnt1.size();
  uInt offset  = itsDataStart + baseline * itsBLDataSize;
  const void* ptr = getReadPointer (blocknr, offset, itsBLDataSize);
  if (itsDoSwap) {
    const char* from = (const char*)ptr;
    char* to = (char*)buf;
    const char* fromend = from + itsBLDataSize;
    while (from < fromend) {
      CanonicalConversion::reverse4 (to, from);
      to += 4;
      from += 4;
    }
  } else {
    memcpy (buf, ptr, itsBLDataSize);
  }
  if (itsVersion < 3) {
    // The first RTCP versions generated conjugate data.
    for (uint i=0; i<itsBLDataSize/sizeof(Complex); ++i) {
      buf[i] = conj(buf[i]);
    }
  }
}

void LofarStMan::putData (uInt rownr, const Complex* buf)
{
  uInt blocknr = rownr / itsAnt1.size();
  uInt baseline = rownr - blocknr*itsAnt1.size();
  uInt offset  = itsDataStart + baseline * itsBLDataSize;
  void* ptr = getWritePointer (blocknr, offset, itsBLDataSize);
  // The first RTCP versions generated conjugate data.
  if (itsVersion < 3) {
    Complex val;
    const Complex* from = buf;
    char* to = (char*)ptr;
    char* toend = to + itsBLDataSize;
    while (to < toend) {
      val = conj(*from);
      if (itsDoSwap) {
        CanonicalConversion::reverse4 (to, &val);
        CanonicalConversion::reverse4 (to + sizeof(float),
                                       ((char*)(&val)) + sizeof(float));
      } else {
        memcpy (to, &val, sizeof(Complex));
      }
      to += sizeof(Complex);
      from++;
    }
  } else if (itsDoSwap) {
    const char* from = (const char*)buf;
    char* to = (char*)ptr;
    const char* fromend = from + itsBLDataSize;
    while (from < fromend) {
      CanonicalConversion::reverse4 (to, from);
      to += 4;
      from += 4;
    }
  } else {
    memcpy (ptr, buf, itsBLDataSize);
  }
  writeData (blocknr, offset, itsBLDataSize);
}

  // NOTE: Nr of samples in observations taken by Cobalt
  //       between 2013-03-28 and 2014-10-28 can erroneously
  //       be > nominal_nsamples (they were set to -1, thus max_unsigned_int).
  //       We fix that setting them to 0.

const uChar* LofarStMan::getNSample1 (uInt rownr, Bool)
{
  uInt blocknr = rownr / itsAnt1.size();
  uInt baseline = rownr - blocknr*itsAnt1.size();
  uInt offset  = itsSampStart + baseline * itsNChan*itsNrBytesPerNrValidSamples;
  const void* ptr = getReadPointer (blocknr, offset, itsNChan*itsNrBytesPerNrValidSamples);
  uChar* to = itsNSampleBuf1.storage();
  memcpy (to, ptr, itsNChan);
  for (uInt i=0; i<itsNChan; ++ i) {
    if (to[i] > itsMaxNrSample) {
      to[i] = 0;
    }
  }
  return to;
}

const uShort* LofarStMan::getNSample2 (uInt rownr, Bool)
{
  uInt blocknr = rownr / itsAnt1.size();
  uInt baseline = rownr - blocknr*itsAnt1.size();
  uInt offset  = itsSampStart + baseline * itsNChan*itsNrBytesPerNrValidSamples;
  const void* ptr = getReadPointer (blocknr, offset, itsNChan*itsNrBytesPerNrValidSamples);
  const uShort* from = (const uShort*)ptr;
  uShort* to = itsNSampleBuf2.storage();
  if (!itsDoSwap) {
    for (uInt i=0; i<itsNChan; ++i) {
      to[i] = (from[i] > itsMaxNrSample  ?  0 : from[i]);
    }
  } else {
    for (uInt i=0; i<itsNChan; ++i) {
      CanonicalConversion::reverse2 (to+i, from+i);
      if (to[i] > itsMaxNrSample) {
        to[i] = 0;
      }
    }
  }
  return to;
}

const uInt* LofarStMan::getNSample4 (uInt rownr, Bool)
{
  uInt blocknr = rownr / itsAnt1.size();
  uInt baseline = rownr - blocknr*itsAnt1.size();
  uInt offset  = itsSampStart + baseline * itsNChan*itsNrBytesPerNrValidSamples;
  const void* ptr = getReadPointer (blocknr, offset, itsNChan*itsNrBytesPerNrValidSamples);
  const uInt* from = (const uInt*)ptr;
  uInt* to = itsNSampleBuf4.storage();
  if (!itsDoSwap) {
    for (uInt i=0; i<itsNChan; ++i) {
      to[i] = (from[i] > itsMaxNrSample  ?  0 : from[i]);
    }
  } else {
    for (uInt i=0; i<itsNChan; ++i) {
      CanonicalConversion::reverse4 (to+i, from+i);
      if (to[i] > itsMaxNrSample) {
        to[i] = 0;
      }
    }
  }
  return to;
}


void* LofarStMan::readFile (uInt blocknr, uInt offset, uInt size)
{
  AlwaysAssert (size <= itsBuffer.size(), AipsError);
  itsRegFile->seek (blocknr*itsBlockSize + offset);
  itsRegFile->read (size, itsBuffer.storage());
  return itsBuffer.storage();
}

void* LofarStMan::getBuffer (uInt size)
{
  AlwaysAssert (size <= itsBuffer.size(), AipsError);
  return itsBuffer.storage();
}

void LofarStMan::writeFile (uInt blocknr, uInt offset, uInt size)
{
  AlwaysAssert (size <= itsBuffer.size(), AipsError);
  itsRegFile->seek (blocknr*itsBlockSize + offset);
  itsRegFile->write (size, itsBuffer.storage());
}

} //# end namespace
