//# tLofarStMan.cc: Test program for the LofarStMan class
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

#include "LofarStMan/LofarStMan.h"
#include <casacore/tables/Tables/TableDesc.h>
#include <casacore/tables/Tables/SetupNewTab.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableLock.h>
#include <casacore/tables/Tables/ScaColDesc.h>
#include <casacore/tables/Tables/ArrColDesc.h>
#include <casacore/tables/Tables/ScalarColumn.h>
#include <casacore/tables/Tables/ArrayColumn.h>
#include <casacore/measures/TableMeasures/TableMeasRefDesc.h>
#include <casacore/measures/TableMeasures/TableMeasValueDesc.h>
#include <casacore/measures/TableMeasures/TableMeasDesc.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Arrays/ArrayMath.h>
#include <casacore/casa/Arrays/ArrayLogical.h>
#include <casacore/casa/Containers/BlockIO.h>
#include <casacore/casa/OS/RegularFile.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/IO/RegularFileIO.h>
#include <casacore/casa/IO/RawIO.h>
#include <casacore/casa/IO/CanonicalIO.h>
#include <casacore/casa/OS/HostInfo.h>
#include <casacore/casa/Exceptions/Error.h>
#include <casacore/casa/iostream.h>
#include <casacore/casa/sstream.h>
#include <casacore/casa/version.h>

using namespace LOFAR;
using namespace casacore;

// This program tests the class LofarStMan and related classes.
// The results are written to stdout. The script executing this program,
// compares the results with the reference output file.
//
// It uses the given phase center, station positions, times, and baselines to
// get known UVW coordinates. In this way the UVW calculation can be checked.


void createField (Table& mainTable)
{
  // Build the table description.
  TableDesc td;
  td.addColumn (ArrayColumnDesc<double>("PHASE_DIR"));
  TableMeasRefDesc measRef(MDirection::J2000);
  TableMeasValueDesc measVal(td, "PHASE_DIR");
  TableMeasDesc<MDirection> measCol(measVal, measRef);
  measCol.write (td);
  // Now create a new table from the description.
  SetupNewTable newtab("tLofarStMan_tmp.data/FIELD", td, Table::New);
  Table tab(newtab, 1);
  ArrayMeasColumn<MDirection> fldcol (tab, "PHASE_DIR");
  Vector<MDirection> phaseDir(1);
  phaseDir[0] = MDirection(MVDirection(Quantity(-1.92653768, "rad"),
                                       Quantity(1.09220917, "rad")),
                           MDirection::J2000);
  fldcol.put (0, phaseDir);
  mainTable.rwKeywordSet().defineTable ("FIELD", tab);
}

void createAntenna (Table& mainTable, unsigned int nant)
{
  // Build the table description.
  TableDesc td;
  td.addColumn (ArrayColumnDesc<double>("POSITION"));
  TableMeasRefDesc measRef(MPosition::ITRF);
  TableMeasValueDesc measVal(td, "POSITION");
  TableMeasDesc<MPosition> measCol(measVal, measRef);
  measCol.write (td);
  // Now create a new table from the description.
  SetupNewTable newtab("tLofarStMan_tmp.data/ANTENNA", td, Table::New);
  Table tab(newtab, nant);
  ScalarMeasColumn<MPosition> poscol (tab, "POSITION");
  // Use the positions of WSRT RT0-3.
  Vector<double> vals(3);
  vals[0] = 3828763; vals[1] = 442449; vals[2] = 5064923;
  poscol.put (0, MPosition(Quantum<Vector<double> >(vals,"m"),
                           MPosition::ITRF));
  vals[0] = 3828746; vals[1] = 442592; vals[2] = 5064924;
  poscol.put (1, MPosition(Quantum<Vector<double> >(vals,"m"),
                           MPosition::ITRF));
  vals[0] = 3828729; vals[1] = 442735; vals[2] = 5064925;
  poscol.put (2, MPosition(Quantum<Vector<double> >(vals,"m"),
                           MPosition::ITRF));
  vals[0] = 3828713; vals[1] = 442878; vals[2] = 5064926;
  poscol.put (3, MPosition(Quantum<Vector<double> >(vals,"m"),
                           MPosition::ITRF));
  // Write the remaining columns with the same position.
  // They are not used in the UVW check.
  for (unsigned int i=4; i<nant; ++i) {
    poscol.put (i, MPosition(Quantum<Vector<double> >(vals,"m"),
                             MPosition::ITRF));
  }
  mainTable.rwKeywordSet().defineTable ("ANTENNA", tab);
}

void createTable (unsigned int nant)
{
  // Build the table description.
  // Add all mandatory columns of the MS main table.
  TableDesc td("", "1", TableDesc::Scratch);
  td.comment() = "A test of class Table";
  td.addColumn (ScalarColumnDesc<double>("TIME"));
  td.addColumn (ScalarColumnDesc<int>("ANTENNA1"));
  td.addColumn (ScalarColumnDesc<int>("ANTENNA2"));
  td.addColumn (ScalarColumnDesc<int>("FEED1"));
  td.addColumn (ScalarColumnDesc<int>("FEED2"));
  td.addColumn (ScalarColumnDesc<int>("DATA_DESC_ID"));
  td.addColumn (ScalarColumnDesc<int>("PROCESSOR_ID"));
  td.addColumn (ScalarColumnDesc<int>("FIELD_ID"));
  td.addColumn (ScalarColumnDesc<int>("ARRAY_ID"));
  td.addColumn (ScalarColumnDesc<int>("OBSERVATION_ID"));
  td.addColumn (ScalarColumnDesc<int>("STATE_ID"));
  td.addColumn (ScalarColumnDesc<int>("SCAN_NUMBER"));
  td.addColumn (ScalarColumnDesc<double>("INTERVAL"));
  td.addColumn (ScalarColumnDesc<double>("EXPOSURE"));
  td.addColumn (ScalarColumnDesc<double>("TIME_CENTROID"));
  td.addColumn (ScalarColumnDesc<bool>("FLAG_ROW"));
  td.addColumn (ArrayColumnDesc<double>("UVW",IPosition(1,3),
                                        ColumnDesc::Direct));
  td.addColumn (ArrayColumnDesc<Complex>("DATA"));
  td.addColumn (ArrayColumnDesc<float>("SIGMA"));
  td.addColumn (ArrayColumnDesc<float>("WEIGHT"));
  td.addColumn (ArrayColumnDesc<float>("WEIGHT_SPECTRUM"));
  td.addColumn (ArrayColumnDesc<bool>("FLAG"));
  td.addColumn (ArrayColumnDesc<bool>("FLAG_CATEGORY"));
  // Now create a new table from the description.
  SetupNewTable newtab("tLofarStMan_tmp.data", td, Table::New);
  // Create the storage manager and bind all columns to it.
  LofarStMan sm1;
  newtab.bindAll (sm1);
  // Finally create the table. The destructor writes it.
  Table tab(newtab);
  // Create the subtables needed to calculate UVW coordinates.
  createField (tab);
  createAntenna (tab, nant);
}

unsigned int nalign (unsigned int size, unsigned int alignment)
{
  return (size + alignment-1) / alignment * alignment - size;
}

void createData (unsigned int nseq, unsigned int nant, unsigned int nchan, unsigned int npol,
                 double startTime, double interval, const Complex& startValue,
                 unsigned int alignment, bool bigEndian, unsigned int myStManVersion,
		 unsigned int myNrBytesPerValidSamples, bool useSeqFile)
{
  AlwaysAssertExit(myStManVersion <= 3);
  // Create the baseline vectors.
  unsigned int nrbl = nant*nant;
  Block<int> ant1(nrbl);
  Block<int> ant2(nrbl);
  unsigned int inx=0;
  for (unsigned int i=0; i<nant; ++i) {
    for (unsigned int j=0; j<nant; ++j) {
      if (myStManVersion != 2) {
        // Use baselines 0,0, 0,1, 1,1 ... n,n
	ant1[inx] = j;
	ant2[inx] = i;
      } else {
        // Use baselines 0,0, 1,0, 1,1 ... n,n (will be swapped by LofarStMan)
	ant1[inx] = i;
	ant2[inx] = j;
      }
      ++inx;
    }
  }
  // Create the meta file. The destructor closes it.
  {
    // If not big-endian, use local format (which might be big-endian).
    if (!bigEndian) {
      bigEndian = HostInfo::bigEndian();
    }
    double maxNSample = 32768;
    AipsIO aio("tLofarStMan_tmp.data/table.f0meta", ByteIO::New);
    aio.putstart ("LofarStMan", myStManVersion);     // version 1, 2, or 3
    aio << ant1 << ant2 << startTime << interval << nchan
        << npol << maxNSample << alignment << bigEndian;
    if (myStManVersion > 1) {
      aio << myNrBytesPerValidSamples;
    } else {
      myNrBytesPerValidSamples = 2;    // version 1 requires 2 bytes
    }
    aio.close();
  }

  cout << "createData:" << endl;
  cout << "  nseq=" << nseq << "  nant=" << nant
       << "  nchan=" << nchan << "  npol=" << npol << endl;
  cout << "  alignment=" << alignment << "  bigEndian=" << bigEndian << endl;
  cout << "  version=" << myStManVersion
       << "  bytesPerSample=" << myNrBytesPerValidSamples
       << "  useSeqFile=" << useSeqFile << endl;

  // Now create the data file.
  auto file = std::make_shared<RegularFileIO>(
    RegularFile("tLofarStMan_tmp.data/table.f0data"), ByteIO::New);
  // Write in canonical (big endian) or local format.
  TypeIO* cfile;
#if CASACORE_MAJOR_VERSION < 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION < 6)
  if (bigEndian) {
    cfile = new CanonicalIO(file.get());
  } else {
    cfile = new RawIO(file.get());
  }
#else
  if (bigEndian) {
    cfile = new CanonicalIO(file);
  } else {
    cfile = new RawIO(file);
  }
#endif

  // Create and initialize data and nsample.
  Array<Complex> data(IPosition(2,npol,nchan));
  indgen (data, startValue, Complex(0.01, 0.01));
  Array<unsigned char>  nsample1(IPosition(1, nchan));
  Array<unsigned short> nsample2(IPosition(1, nchan));
  Array<unsigned int>   nsample4(IPosition(1, nchan));
  indgen (nsample1);
  indgen (nsample2);
  indgen (nsample4);
  nsample2(IPosition(1,nchan-1)) = 33000; // invalid nsample for Cobalt bugfix
  nsample4(IPosition(1,nchan-1)) = 33000; // invalid nsample for Cobalt bugfix

  // Allocate space for possible block alignment.
  if (alignment < 1) {
    alignment = 1;
  }
  unsigned int seqSize = (myStManVersion==1 ? 4:8);
  Block<char> align1(nalign(seqSize, alignment), 0);
  Block<char> align2(nalign(nrbl*8*data.size(), alignment), 0);

  unsigned int nsamplesSize = nrbl*myNrBytesPerValidSamples*nsample2.size();
  Block<char> align3(nalign(nsamplesSize, alignment), 0);

  // Write the data as nseq blocks.
  for (unsigned int i=0; i<nseq; ++i) {
    if (myStManVersion > 1) {
      // From version 2 on RTCP writes a magic value before the seqnr.
      unsigned int magicVal = 0x0000da7a;
      cfile->write (1, &magicVal);
    }
    cfile->write (1, &i);
    if (align1.size() > 0) {
      cfile->write (align1.size(), align1.storage());
    }

    for (unsigned int j=0; j<nrbl; ++j) {
      // The RTCP wrote the conj of the data for version 1 and 2.
      if (myStManVersion < 3) {
        Array<Complex> cdata = conj(data);
        cfile->write (data.size(), cdata.data());
      } else {
        cfile->write (data.size(), data.data());
      }
      data += Complex(0.01, 0.02);
    }
    if (align2.size() > 0) {
      cfile->write (align2.size(), align2.storage());
    }

    if (myStManVersion < 2) {
      for (unsigned int j=0; j<nrbl; ++j) {
	cfile->write (nsample2.size(), nsample2.data());
	nsample2 += static_cast<unsigned short>(1);
      }
    } else {
      for (unsigned int j=0; j<nrbl; ++j) {
	switch (myNrBytesPerValidSamples) {
	case 1:
	  {
	    cfile->write(nsample1.size(), nsample1.data());
	    nsample1 += static_cast<unsigned char>(1);
	  } break;
	case 2:
	  {
	    cfile->write(nsample2.size(), nsample2.data());
	    nsample2 += static_cast<unsigned short>(1);
	  } break;
	case 4:
	  {
	    cfile->write(nsample4.size(), nsample4.data());
	    nsample4 += 1u;
	  } break;
	}
      }
    }

    if (align3.size() > 0) {
      cfile->write (align3.size(), align3.storage());
    }
  }
  delete cfile;

  if (useSeqFile  &&  myStManVersion > 1) {
    TypeIO* sfile = 0;
    // create seperate file for sequence numbers if version > 1
    auto file = std::make_shared<RegularFileIO>(
      RegularFile("tLofarStMan_tmp.data/table.f0seqnr"), ByteIO::New);
    // Write in canonical (big endian) or local format.
#if CASACORE_MAJOR_VERSION < 3 || \
    (CASACORE_MAJOR_VERSION == 3 && CASACORE_MINOR_VERSION < 6)
    if (bigEndian) {
      sfile = new CanonicalIO(file.get());
    } else {
      sfile = new RawIO(file.get());
    }
#else
    if (bigEndian) {
      sfile = new CanonicalIO(file);
    } else {
      sfile = new RawIO(file);
    }
#endif
    for (unsigned int i=0; i<nseq; ++i) {
      sfile->write (1, &i);
    }
    delete sfile;
  } else {
    // Delete a possible existing seqfile.
    RegularFile file("tLofarStMan_tmp.data/table.f0seqnr");
    file.remove();
  }
}


void checkUVW (unsigned int row, unsigned int nant, Vector<double> uvw)
{
  // Expected outcome of UVW for antenna 0-3 and seqnr 0-1
  static double uvwVals[] = {
    0, 0, 0,
    0.423756, -127.372, 67.1947,
    0.847513, -254.744, 134.389,
    0.277918, -382.015, 201.531,
    -0.423756, 127.372, -67.1947,
    0, 0, 0,
    0.423756, -127.372, 67.1947,
    -0.145838, -254.642, 134.336,
    -0.847513, 254.744, -134.389,
    -0.423756, 127.372, -67.1947,
    0, 0, 0,
    -0.569594, -127.27, 67.1417,
    -0.277918, 382.015, -201.531,
    0.145838, 254.642, -134.336,
    0.569594, 127.27, -67.1417,
    0, 0, 0,
    0, 0, 0,
    0.738788, -127.371, 67.1942,
    1.47758, -254.742, 134.388,
    1.22276, -382.013, 201.53,
    -0.738788, 127.371, -67.1942,
    0, 0, 0,
    0.738788, -127.371, 67.1942,
    0.483976, -254.642, 134.336,
    -1.47758, 254.742, -134.388,
    -0.738788, 127.371, -67.1942,
    0, 0, 0,
    -0.254812, -127.271, 67.1421,
    -1.22276, 382.013, -201.53,
    -0.483976, 254.642, -134.336,
    0.254812, 127.271, -67.1421,
    0, 0, 0
  };
  unsigned int nrbl = nant*nant;
  unsigned int seqnr = row / nrbl;
  unsigned int bl = row % nrbl;
  unsigned int ant1 = bl % nant;
  unsigned int ant2 = bl / nant;
  // Only check first two time stamps and first four antennae.
  if (seqnr < 2  &&  ant1 < 4  &&  ant2 < 4) {
    AlwaysAssertExit (near(uvw[0],
                           uvwVals[3*(seqnr*16 + 4*ant1 + ant2)],
                           1e-5))
  }
}

// maxWeight tells maximum weight before it wraps
// (when nbytesPerSample is small).
void readTable (unsigned int nseq, unsigned int nant, unsigned int nchan, unsigned int npol,
                double startTime, double interval, const Complex& startValue,
                float maxWeight)
{
  unsigned int nbasel = nant*nant;
  // Open the table and check if #rows is as expected.
  Table tab("tLofarStMan_tmp.data");
  unsigned int nrow = tab.nrow();
  AlwaysAssertExit (nrow == nseq*nbasel);
  AlwaysAssertExit (!tab.canAddRow());
  AlwaysAssertExit (!tab.canRemoveRow());
  AlwaysAssertExit (tab.canRemoveColumn(Vector<String>(1, "DATA")));
  // Create objects for all mandatory MS columns.
  ArrayColumn<Complex> dataCol(tab, "DATA");
  ArrayColumn<float> weightCol(tab, "WEIGHT");
  ArrayColumn<float> wspecCol(tab, "WEIGHT_SPECTRUM");
  ArrayColumn<float> sigmaCol(tab, "SIGMA");
  ArrayColumn<double> uvwCol(tab, "UVW");
  ArrayColumn<bool> flagCol(tab, "FLAG");
  ArrayColumn<bool> flagcatCol(tab, "FLAG_CATEGORY");
  ScalarColumn<double> timeCol(tab, "TIME");
  ScalarColumn<double> centCol(tab, "TIME_CENTROID");
  ScalarColumn<double> intvCol(tab, "INTERVAL");
  ScalarColumn<double> expoCol(tab, "EXPOSURE");
  ScalarColumn<int> ant1Col(tab, "ANTENNA1");
  ScalarColumn<int> ant2Col(tab, "ANTENNA2");
  ScalarColumn<int> feed1Col(tab, "FEED1");
  ScalarColumn<int> feed2Col(tab, "FEED2");
  ScalarColumn<int> ddidCol(tab, "DATA_DESC_ID");
  ScalarColumn<int> pridCol(tab, "PROCESSOR_ID");
  ScalarColumn<int> fldidCol(tab, "FIELD_ID");
  ScalarColumn<int> arridCol(tab, "ARRAY_ID");
  ScalarColumn<int> obsidCol(tab, "OBSERVATION_ID");
  ScalarColumn<int> stidCol(tab, "STATE_ID");
  ScalarColumn<int> scnrCol(tab, "SCAN_NUMBER");
  ScalarColumn<bool> flagrowCol(tab, "FLAG_ROW");
  // Create and initialize expected data and weight.
  Array<Complex> dataExp(IPosition(2,npol,nchan));
  indgen (dataExp, startValue, Complex(0.01, 0.01));
  Array<float> weightExp(IPosition(2,1,nchan));
  indgen (weightExp);
  // Loop through all rows in the table and check the data.
  unsigned int row=0;
  for (unsigned int i=0; i<nseq; ++i) {

    for (unsigned int j=0; j<nant; ++j) {
      for (unsigned int k=0; k<nant; ++k) {

        // Wrap expected weight if needed.
        for (unsigned int i=0; i<weightExp.size(); ++i) {
          while (weightExp.data()[i] >= maxWeight) {
            weightExp.data()[i] -= maxWeight;
          }
        }
        if (maxWeight > 32768) {
          weightExp(IPosition(2,0,nchan-1)) = 0;   // invalid nsample
        }
        // Contents must be present except for FLAG_CATEGORY.
	AlwaysAssertExit (dataCol.isDefined (row));
	AlwaysAssertExit (weightCol.isDefined (row));
        AlwaysAssertExit (wspecCol.isDefined (row));
        AlwaysAssertExit (sigmaCol.isDefined (row));
        AlwaysAssertExit (flagCol.isDefined (row));
        AlwaysAssertExit (!flagcatCol.isDefined (row));
        // Check data, weight, sigma, weight_spectrum, flag
        AlwaysAssertExit (allNear (dataCol(row), dataExp, 1e-7));
        AlwaysAssertExit (weightCol.shape(row) == IPosition(1,npol));
        AlwaysAssertExit (allEQ (weightCol(row), float(1)));
        AlwaysAssertExit (sigmaCol.shape(row) == IPosition(1,npol));
        AlwaysAssertExit (allEQ (sigmaCol(row), float(1)));
        Array<float> weights = wspecCol(row);
        AlwaysAssertExit (weights.shape() == IPosition(2,npol,nchan));
        for (unsigned int p=0; p<npol; ++p) {

	  if (!allNear(weights(IPosition(2,p,0),IPosition(2,p,nchan-1)),
                       weightExp/float(32768), 1e-7)) {

	    std::cout << "weights: " << std::endl;
	    std::cout << weights(IPosition(2,p,0), IPosition(2,p,nchan-1)) << std::endl;

	    std::cout << "weigthExp: " << std::endl;
	    std::cout << weightExp/float(32768) << std:: endl;

	  }
        }
        Array<bool> flagExp (weights == float(0));
        AlwaysAssertExit (allEQ (flagCol(row), flagExp));
        // Check ANTENNA1 and ANTENNA2

        AlwaysAssertExit (ant1Col(row) == int(k));
        AlwaysAssertExit (ant2Col(row) == int(j));

        dataExp += Complex(0.01, 0.02);
        weightExp += float(1);
        ++row;
      }
    }
  }
  // Check values in TIME column.
  Vector<double> times = timeCol.getColumn();
  AlwaysAssertExit (times.size() == nrow);
  row=0;
  startTime += interval/2;
  for (unsigned int i=0; i<nseq; ++i) {
    for (unsigned int j=0; j<nbasel; ++j) {
      AlwaysAssertExit (near(times[row], startTime));
      ++row;
    }
    startTime += interval;
  }
  // Check the other columns.
  AlwaysAssertExit (allNear(centCol.getColumn(), times, 1e-13));
  AlwaysAssertExit (allNear(intvCol.getColumn(), interval, 1e-13));
  AlwaysAssertExit (allNear(expoCol.getColumn(), interval, 1e-13));
  AlwaysAssertExit (allEQ(feed1Col.getColumn(), 0));
  AlwaysAssertExit (allEQ(feed2Col.getColumn(), 0));
  AlwaysAssertExit (allEQ(ddidCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(pridCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(fldidCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(arridCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(obsidCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(stidCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(scnrCol.getColumn(), 0));
  AlwaysAssertExit (allEQ(flagrowCol.getColumn(), false));
  // Check the UVW coordinates.
  for (unsigned int i=0; i<nrow; ++i) {
    checkUVW (i, nant, uvwCol(i));
  }
  RefRows rownrs(0,2,1);
  Slicer slicer(IPosition(2,0,0), IPosition(2,1,1));
  cout << wspecCol(0).shape() << endl;
  Array<float> wg = wspecCol.getColumnCells (rownrs);
}

void updateTable (unsigned int nchan, unsigned int npol, const Complex& startValue)
{
  // Open the table for write.
  Table tab("tLofarStMan_tmp.data", Table::Update);
  unsigned int nrow = tab.nrow();
  // Create object for DATA column.
  ArrayColumn<Complex> dataCol(tab, "DATA");
  // Check we can write the column, but not change the shape.
  AlwaysAssertExit (tab.isColumnWritable ("DATA"));
  AlwaysAssertExit (!dataCol.canChangeShape());
  // Create and initialize data.
  Array<Complex> data(IPosition(2,npol,nchan));
  indgen (data, startValue, Complex(0.01, 0.01));
  // Loop through all rows in the table and write the data.
  for (unsigned int row=0; row<nrow; ++row) {
    dataCol.put (row, data);
    data += Complex(0.01, 0.02);
  }
}

void copyTable()
{
  Table tab("tLofarStMan_tmp.data");
  // Deep copy the table.
  tab.deepCopy ("tLofarStMan_tmp.datcp", Table::New, true);
}


int main (int argc, char* argv[])
{
  try {
    // Register LofarStMan to be able to read it back.
    LofarStMan::registerClass();
    // Get nseq, nant, nchan, npol from argv.
    unsigned int nseq=10;
    unsigned int nant=16;
    unsigned int nchan=64;
    unsigned int npol=4;
    if (argc > 1) {
      istringstream istr(argv[1]);
      istr >> nseq;
    }
    if (nseq == 0) {
      Table tab(argv[2]);
      cout << "nrow=" << tab.nrow() << endl;
      ArrayColumn<double> uvwcol(tab, "UVW");
      cout << "uvws="<< uvwcol(0) << endl;
      cout << "uvws="<< uvwcol(1) << endl;
      cout << "uvwe="<< uvwcol(tab.nrow()-1) << endl;
    } else {
      if (argc > 2) {
        istringstream istr(argv[2]);
        istr >> nant;
      }
      if (argc > 3) {
        istringstream istr(argv[3]);
        istr >> nchan;
      }
      if (argc > 4) {
        istringstream istr(argv[4]);
        istr >> npol;
      }
      // Test all possible bytes per sample.
      unsigned int nbytesPerSample[] = {0, 2,4,1};
      unsigned int maxWeight[] = {0, 256*256, 256*256*256*127, 256};
      // Test the various versions.
      for (int v=1; v<4; ++v) {
        cout << "Test version " << v << endl;
        // Create the table.
        createTable (nant);
        // Write data in big-endian and check it. Align on 512.
        // Use this start time and interval to get known UVWs.
        // They are the same as used in DPPP/tUVWFlagger.
        double interval= 30.;
        double startTime = 4472025740.0 - interval*0.5;
        createData (nseq, nant, nchan, npol, startTime, interval,
                    Complex(0.1, 0.1), 512, true, v, nbytesPerSample[v], v%2==0);
        readTable (nseq, nant, nchan, npol, startTime, interval,
                   Complex(0.1, 0.1), maxWeight[v]);
        // Update the table and check again.
        updateTable (nchan, npol, Complex(-3.52, -20.3));
        readTable (nseq, nant, nchan, npol, startTime, interval,
                   Complex(-3.52, -20.3), maxWeight[v]);
        // Write data in local format and check it. No alignment.
        createData (nseq, nant, nchan, npol, startTime, interval,
                    Complex(3.1, -5.2), 0, false, v, nbytesPerSample[v], v%2!=0);
        readTable (nseq, nant, nchan, npol, startTime, interval,
                   Complex(3.1, -5.2), maxWeight[v]);
        // Update the table and check again.
        updateTable (nchan, npol, Complex(3.52, 20.3));
        readTable (nseq, nant, nchan, npol, startTime, interval,
                   Complex(3.52, 20.3), maxWeight[v]);
        copyTable();
      }
    }
  } catch (AipsError& x) {
    cout << "Caught an exception: " << x.getMesg() << endl;
    return 1;
  }
  return 0;                           // exit with success status
}
