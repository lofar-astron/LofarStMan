//# LofarStMan.h: Storage Manager for the main table of a LOFAR MS
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

#ifndef LOFAR_LOFARSTMAN_LOFARSTMAN_H
#define LOFAR_LOFARSTMAN_LOFARSTMAN_H

//# Includes
#include <casacore/tables/DataMan/DataManager.h>
#include <casacore/casa/IO/MMapIO.h>
#include <casacore/casa/IO/FiledesIO.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/Containers/Record.h>

#include <cstdint>
#include <vector>

namespace LOFAR {

//# Forward Declarations.
class LofarColumn;

// <summary>
// The Storage Manager for the main table of a raw LOFAR MS
// </summary>

// <use visibility=export>

// <reviewed reviewer="UNKNOWN" date="before2004/08/25" tests="tLofarStMan.cc">
// </reviewed>

// <prerequisite>
//# Classes you should understand before using this one.
//   <li> The Table Data Managers concept as described in module file
//        <linkto module="Tables:Data Managers">Tables.h</linkto>
// </prerequisite>

// <etymology>
// LofarStMan is the data manager which stores the data for a LOFAR MS.
// </etymology>

// <synopsis>
// LofarStMan is a specific storage manager for the main table of a LOFAR MS.
// For performance purposes the raw data from the correlator is directly
// written to a disk file. However, to be able to use the data directly as a
// MeasurementSet, this specific storage manager is created offering access to
// all mandatory columns in the main table of the MS.
//
// Similar to other storage managers, the LofarStMan files need to be part of
// the table directory. There are two files:
// <ul>
//  <li> The meta file contains the meta data describing baselines, start time,
//       integration time, etc. It needs to be written as an AipsIO file.
//       The meta info should also tell the endianness of the data file.
//  <li> The data file consists of NSEQ data blocks each containing:
//   <ul>
//    <li> 4-byte sequence number defining the time stamp.
//    <li> Complex data with shape [npol,nchan,nbasel].
//    <li> Unsigned short nr of samples used in each data point. It has shape
//         [nchan,nbasel]. It defines WEIGHT_SPECTRUM and FLAG.
//    <li> Filler bytes to align the blocks as given in the meta info.
//   </ul>
//   The sequence numbers are ascending, but there can be holes due to
//   missing time stamps.
// </ul>
// The first versions of the data file can only handle regularly shaped data
// with equal integration times. A future version might be able to deal with
// varying integration times (depending on baseline length).
//
// Most of the MS columns (like DATA_DESC_ID) are not stored in the data file;
// usually they map to the value 0. This is also true for the UVW column, so
// the UVW coordinates need to be added to the table in a separate step because
// the online system does not have the resources to do it.
//
// All columns are readonly with the exception of DATA.
// </synopsis>

// <motivation>
// The common Table storage managers are too slow for the possibly high
// output rate of the LOFAR correlator.
// </motivation>

// <example>
// The following example shows how to create a table and how to attach
// the storage manager to some columns.
// <srcblock>
//   SetupNewTable newtab("name.data", tableDesc, Table::New);
//   LofarStMan stman;                     // define storage manager
//   newtab.bindColumn ("DATA", stman);    // bind column to st.man.
//   newtab.bindColumn ("FLAG", stman);    // bind column to st.man.
//   Table tab(newtab);                    // actually create table
// </srcblock>
// </example>

//# <todo asof="$DATE:$">
//# A List of bugs, limitations, extensions or planned refinements.
//# </todo>


class LofarStMan final : public casacore::DataManager
{
public:
    // Create a Lofar storage manager with the given name.
    // If no name is used, it is set to "LofarStMan"
  explicit LofarStMan (const casacore::String& dataManagerName = "LofarStMan");

  // Create a Lofar storage manager with the given name.
  // The specifications are part of the record (as created by dataManagerSpec).
  LofarStMan (const casacore::String& dataManagerName, const casacore::Record& spec);
  
  ~LofarStMan();

  // Clone this object.
  virtual casacore::DataManager* clone() const override;
  
  // Get the type name of the data manager (i.e. LofarStMan).
  virtual casacore::String dataManagerType() const override;
  
  // Get the name given to the storage manager (in the constructor).
  virtual casacore::String dataManagerName() const override;
  
  // Record a record containing data manager specifications.
  virtual casacore::Record dataManagerSpec() const override;

  // Get the number of rows in this storage manager.
  unsigned int getNRow() const
    { return itsNrRows; }
  
  // The storage manager is not a regular one.
  bool isRegular() const override;
  
  // The storage manager cannot add rows.
  bool canAddRow() const override;
  
  // The storage manager cannot delete rows.
  bool canRemoveRow() const override;
  
  // The storage manager can add columns, which does not really do something.
  bool canAddColumn() const override;
  
  // Columns can be removed, but it does not do anything at all.
  bool canRemoveColumn() const override;
  
  // Make the object from the type name string.
  // This function gets registered in the DataManager "constructor" map.
  // The caller has to delete the object.
  static casacore::DataManager* makeObject (const casacore::String& aDataManType,
                                        const casacore::Record& spec);

  // Register the class name and the static makeObject "constructor".
  // This will make the engine known to the table system.
  static void registerClass();


  // Get data.
  // <group>
  const casacore::Block<int>& ant1() const
    { return itsAnt1; }
  const casacore::Block<int>& ant2() const
    { return itsAnt2; }
  double time (unsigned int blocknr);
  double interval() const
    { return itsTimeIntv; }
  unsigned int nchan() const
    { return itsNChan; }
  unsigned int npol() const
    { return itsNPol; }
  double maxnSample() const
    { return itsMaxNrSample; }
  void getData (unsigned int rownr, casacore::Complex* buf);
  void putData (unsigned int rownr, const casacore::Complex* buf);

  const unsigned char*  getNSample1 (unsigned int rownr, bool swapIfNeeded);
  const unsigned short* getNSample2 (unsigned int rownr, bool swapIfNeeded);
  const unsigned int*   getNSample4 (unsigned int rownr, bool swapIfNeeded);
  // </group>

  unsigned int getLofarStManVersion() const
    { return itsVersion; }

  unsigned int getNrBytesPerNrValidSamples() const
    { return itsNrBytesPerNrValidSamples; }

private:
  // Copy constructor cannot be used.
  LofarStMan (const LofarStMan& that);

  // Assignment cannot be used.
  LofarStMan& operator= (const LofarStMan& that);
  
  // Flush and optionally fsync the data.
  // It does nothing, and returns false.
  bool flush (casacore::AipsIO&, bool doFsync) override;
  
  // Let the storage manager create files as needed for a new table.
  // This allows a column with an indirect array to create its file.
  void create (unsigned int nrrow) override;
  
  // Open the storage manager file for an existing table.
  // Return the number of rows in the data file.
  // <group>
  void open (unsigned int nrrow, casacore::AipsIO&) override; //# should never be called
  unsigned int open1 (unsigned int nrrow, casacore::AipsIO&) override;
  // </group>

  // Prepare the columns (needed for UvwColumn).
  void prepare() override;

  // Resync the storage manager with the new file contents.
  // It does nothing.
  // <group>
  void resync (unsigned int nrrow) override;   //# should never be called
  unsigned int resync1 (unsigned int nrrow) override;
  // </group>
  
  // Reopen the storage manager files for read/write.
  // It does nothing.
  void reopenRW() override;
  
  // The data manager will be deleted (because all its columns are
  // requested to be deleted).
  // So clean up the things needed (e.g. delete files).
  void deleteManager() override;

  // Add rows to the storage manager.
  // It cannot do it, so throws an exception.
  void addRow (unsigned int nrrow) override;
  
  // Delete a row from all columns.
  // It cannot do it, so throws an exception.
  void removeRow (unsigned int rowNr) override;
  
  // Do the final addition of a column.
  // It won't do anything.
  void addColumn (casacore::DataManagerColumn*) override;
  
  // Remove a column from the data file.
  // It won't do anything.
  void removeColumn (casacore::DataManagerColumn*) override;
  
  // Create a column in the storage manager on behalf of a table column.
  // The caller has to delete the newly created object.
  // <group>
  // Create a scalar column.
  casacore::DataManagerColumn* makeScalarColumn (const casacore::String& aName,
					       int aDataType,
					       const casacore::String& aDataTypeID) override;
  // Create a direct array column.
  casacore::DataManagerColumn* makeDirArrColumn (const casacore::String& aName,
					       int aDataType,
					       const casacore::String& aDataTypeID) override;
  // Create an indirect array column.
  casacore::DataManagerColumn* makeIndArrColumn (const casacore::String& aName,
					       int aDataType,
					       const casacore::String& aDataTypeID) override;
  // </group>

  // Initialize by reading the header info.
  void init();

  // Open the data file and seqnr file.
  // The seqnr is always memory-mapped (it is very small).
  // The data file is only memory-mapped in 64 bit systems because the
  // address space of 32-bit systems is too small for it.
  void openFiles (bool writable);

  // Memory map the seqnr file.
  void mapSeqFile();

  // Close all files.
  void closeFiles();

  // Get a pointer to data to be read.
  const void* getReadPointer (unsigned int blocknr, unsigned int offset,
                              unsigned int size)
  {
    return readFile (blocknr, offset, size);
  }

  // Get a pointer where data can be written.
  void* getWritePointer (unsigned int /*blocknr*/, unsigned int /*offset*/,
                         unsigned int size)
  {
    return getBuffer (size);
  }

  // Write the data. It is a no-op if mmap is used.
  void writeData (unsigned int blocknr, unsigned int offset, unsigned int size)
  {
    writeFile (blocknr, offset, size);
  }

  // Read or write the data for regular files.
  void* readFile  (unsigned int blocknr, unsigned int offset, unsigned int size);
  void* getBuffer (unsigned int size);
  void  writeFile (unsigned int blocknr, unsigned int offset, unsigned int size);


  //# Declare member variables.
  // Name of data manager.
  casacore::String itsDataManName;
  // The number of rows in the columns.
  unsigned int         itsNrRows;
  // The antennae forming the baselines.
  casacore::Block<int> itsAnt1;
  casacore::Block<int> itsAnt2;
  // The start time and interval.
  double itsStartTime;
  double itsTimeIntv;
  unsigned int itsNChan;
  unsigned int itsNPol;
  unsigned int itsNrBytesPerNrValidSamples;
  // The column objects.
  std::vector<LofarColumn*> itsColumns;
  // On 32-bit systems regular IO is used.
  int    itsFD;
  casacore::FiledesIO* itsRegFile;
  casacore::Block<char> itsBuffer;   //# buffer of size itsBLDataSize for regular IO
  // The seqnr file (if present) is always memory-mapped because it is small.
  casacore::MMapIO*     itsSeqFile;
  bool   itsDoSwap;       //# true = byte-swapping is needed
  int64_t itsBlockSize;    //# size of a block containing a seqnr
  int64_t itsBLDataSize;   //# data size of a single baseline
  int64_t itsDataStart;    //# start of data in a block
  int64_t itsSampStart;    //# start of nsamples in a block
  //# Buffer to hold nsample values.
  casacore::Block<unsigned char> itsNSampleBuf1;
  casacore::Block<unsigned short> itsNSampleBuf2;
  casacore::Block<unsigned int>   itsNSampleBuf4;
  double  itsMaxNrSample; //# weight = nsample / itsMaxNrSample;
  casacore::Record itsSpec;

  unsigned int itsVersion;        //# Version of LofarStMan MeasurementSet
};


} //# end namespace

#endif
