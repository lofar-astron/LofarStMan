//# LofarColumn.cc: A Column in the LOFAR Storage Manager
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

#include <LofarStMan/LofarColumn.h>

#include <casacore/tables/DataMan/DataManError.h>
#include <casacore/tables/Tables/Table.h>
#include <casacore/tables/Tables/TableRecord.h>
#include <casacore/measures/TableMeasures/ScalarMeasColumn.h>
#include <casacore/measures/TableMeasures/ArrayMeasColumn.h>
#include <casacore/measures/Measures/MEpoch.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MCDirection.h>
#include <casacore/measures/Measures/MPosition.h>
#include <casacore/measures/Measures/MCPosition.h>
#include <casacore/measures/Measures/MBaseline.h>
#include <casacore/measures/Measures/MCBaseline.h>
#include <casacore/measures/Measures/MeasConvert.h>
#include <casacore/measures/Measures/MeasTable.h>
#include <casacore/measures/Measures/Muvw.h>
#include <casacore/casa/Arrays/Array.h>
#include <casacore/casa/Utilities/Assert.h>
#include <casacore/casa/Exceptions/Error.h>

using namespace casacore;


namespace LOFAR {

  LofarColumn::~LofarColumn()
  {}
  bool LofarColumn::isWritable() const
  {
    return false;
  }
  void LofarColumn::setShapeColumn (const IPosition&)
  {}
  void LofarColumn::prepareCol()
  {}

  Ant1Column::~Ant1Column()
  {}
  void Ant1Column::getIntV (unsigned int rownr, int* dataPtr)
  {
    // Fill ColumnCache object.
    const Block<int>& ants = itsParent->ant1();
    columnCache().setIncrement (1);
    unsigned int strow = rownr / ants.size() * ants.size();
    columnCache().setIncrement (1);
    columnCache().set (strow, strow + ants.size() - 1, ants.storage());
    *dataPtr = ants[rownr-strow];
  }

  Ant2Column::~Ant2Column()
  {}
  void Ant2Column::getIntV (unsigned int rownr, int* dataPtr)
  {
    // Fill ColumnCache object.
    const Block<int>& ants = itsParent->ant2();
    unsigned int strow = rownr / ants.size() * ants.size();
    columnCache().setIncrement (1);
    columnCache().set (strow, strow + ants.size() - 1, ants.storage());
    *dataPtr = ants[rownr-strow];
  }

  TimeColumn::~TimeColumn()
  {}
  void TimeColumn::getdoubleV (unsigned int rownr, double* dataPtr)
  {
    // Get time of the block containing this row.
    unsigned int nrbasel = itsParent->ant1().size();
    unsigned int blnr = rownr / nrbasel;
    itsValue = itsParent->time (blnr);
    // Fill ColumnCache object.
    unsigned int strow = blnr * nrbasel;
    columnCache().setIncrement (0);
    columnCache().set (strow, strow + nrbasel - 1, &itsValue);
    *dataPtr = itsValue;
  }

  IntervalColumn::~IntervalColumn()
  {}
  void IntervalColumn::getdoubleV (unsigned int, double* dataPtr)
  {
    itsValue = itsParent->interval();
    columnCache().setIncrement (0);
    columnCache().set (0, itsParent->getNRow()-1, &itsValue);
    *dataPtr = itsValue;
  }

  ZeroColumn::~ZeroColumn()
  {}
  void ZeroColumn::getIntV (unsigned int, int* dataPtr)
  {
    itsValue = 0;
    columnCache().setIncrement (0);
    columnCache().set (0, itsParent->getNRow()-1, &itsValue);
    *dataPtr = 0;
  }

  FalseColumn::~FalseColumn()
  {}
  void FalseColumn::getBoolV (unsigned int, bool* dataPtr)
  {
    itsValue = false;
    columnCache().setIncrement (0);
    columnCache().set (0, itsParent->getNRow()-1, &itsValue);
    *dataPtr = 0;
  }

  UvwColumn::~UvwColumn()
  {}
  void UvwColumn::prepareCol()
  {
    // Read the station positions from the ANTENNA subtable
    // and convert them to a baseline in ITRF.
    const TableRecord& keyset = itsParent->table().keywordSet();
    itsCanCalc = keyset.isDefined ("ANTENNA");
    if (itsCanCalc) {
      Table anttab (keyset.asTable ("ANTENNA"));
      AlwaysAssert (anttab.nrow() > 0, AipsError);
      int nrant = anttab.nrow();
      ROScalarMeasColumn<MPosition> antcol (anttab, "POSITION");
      MPosition arrayPos;
      Vector<double> pos0;
      for (int i=0; i<nrant; ++i) {
        // Read antenna position and convert to ITRF.
        MPosition mpos = MPosition::Convert (antcol(i), MPosition::ITRF)();
        if (i == 0) {
          pos0 = mpos.getValue().getVector();
        }
        // Use position of middle station as array position.
        if (i == nrant/2) {
          arrayPos = mpos;
        }
        Vector<double> pos = mpos.getValue().getVector();
        MVPosition mvpos((pos[0] - pos0[0]),
                         (pos[1] - pos0[1]),
                         (pos[2] - pos0[2]));
        itsAntMB.push_back (MBaseline (MVBaseline(mvpos), MBaseline::ITRF));
      }
      // Read the phase reference position from the FIELD subtable.
      // Only use the first value from the PHASE_DIR array.
      Table fldtab (itsParent->table().keywordSet().asTable ("FIELD"));
      AlwaysAssert (fldtab.nrow() == 1, AipsError);
      ROArrayMeasColumn<MDirection> fldcol (fldtab, "PHASE_DIR");
      itsPhaseDir = fldcol(0).data()[0];
      // Create a reference frame. Use the middle antenna as array position.
      itsFrame.set (arrayPos);
      // Initialize the rest which is used to cache the UVW per antenna.
      // The cache is only useful if the MS is accessed in time order, but that
      // is normally the case.
      itsLastBlNr = -1;
      itsAntUvw.resize (nrant);
      itsUvwFilled.resize (nrant);
      itsUvwFilled = false;
    }
  }
  IPosition UvwColumn::shape (unsigned int)
  {
    return IPosition(1,3);
  }
  void UvwColumn::getArraydoubleV (unsigned int rownr, Array<double>* dataPtr)
  {
    if (!itsCanCalc) {
      *dataPtr = 0.;
    } else {
      // Get nr of the block containing this row.
      int nrbasel = itsParent->ant1().size();
      int blnr   = rownr / nrbasel;
      int antinx = rownr - blnr * nrbasel;
      int ant1   = itsParent->ant1()[antinx];
      int ant2   = itsParent->ant2()[antinx];
      // If a different block (i.e. time), we have to calculate the UVWs.
      if (blnr != itsLastBlNr) {
        itsLastBlNr  = blnr;
        Quantum<double> tm(itsParent->time(blnr), "s");
        itsFrame.set (MEpoch(MVEpoch(tm.get("d").getValue()), MEpoch::UTC));
        itsJ2000Dir = MDirection::Convert (itsPhaseDir,
                                           MDirection::Ref(MDirection::J2000,
                                                           itsFrame))();
        itsFrame.set (itsJ2000Dir);
        itsUvwFilled = false;
      }
      // Calculate the UVWs for this timestamp if not done yet.
      int ant = ant1;
      for (int i=0; i<2; ++i) {
        if (!itsUvwFilled[ant]) {
          MBaseline& mbl = itsAntMB[ant];
          mbl.getRefPtr()->set(itsFrame);       // attach frame
          MBaseline::Convert mcvt(mbl, MBaseline::J2000);
          MVBaseline bas = mcvt().getValue();
          MVuvw jvguvw(bas, itsJ2000Dir.getValue());
          itsAntUvw[ant] = Muvw(jvguvw, Muvw::J2000).getValue().getVector();
          itsUvwFilled[ant] = true;
        }
        ant = ant2;
      }
      // The UVW of the baseline is the difference of the antennae.
      *dataPtr = itsAntUvw[ant2] - itsAntUvw[ant1];
    }
  }

  DataColumn::~DataColumn()
  {}
  bool DataColumn::isWritable() const
  {
    return true;
  }
  IPosition DataColumn::shape (unsigned int)
  {
    return IPosition(2, itsParent->npol(), itsParent->nchan());
  }
  void DataColumn::getArrayComplexV (unsigned int rownr, Array<Complex>* dataPtr)
  {
    bool deleteIt;
    Complex* data = dataPtr->getStorage(deleteIt);
    itsParent->getData (rownr, data);
    dataPtr->putStorage (data, deleteIt);
  }
  void DataColumn::putArrayComplexV (unsigned int rownr, const Array<Complex>* dataPtr)
  {
    bool deleteIt;
    const Complex* data = dataPtr->getStorage(deleteIt);
    itsParent->putData (rownr, data);
    dataPtr->freeStorage (data, deleteIt);
  }

  FlagColumn::~FlagColumn()
  {}
  IPosition FlagColumn::shape (unsigned int)
  {
    return IPosition(2, itsParent->npol(), itsParent->nchan());
  }
  void FlagColumn::getArrayBoolV (unsigned int rownr, Array<bool>* dataPtr)
  {
    unsigned int npol = itsParent->npol();
    
    switch(itsParent->getLofarStManVersion()) {
    case 1:
    {
      const unsigned short* data = itsParent->getNSample2 (rownr, false);
      const unsigned short* dataEnd = data + itsParent->nchan();

      if (dataPtr->contiguousStorage()) {
	for (Array<bool>::contiter iter=dataPtr->cbegin(); data<dataEnd; ++data) {
	  bool flagged = (*data == 0);
	  for (unsigned int i=0; i<npol; ++i, ++iter) {
	    *iter = flagged;
	  }
	}
      } else {
	for (Array<bool>::iterator iter=dataPtr->begin();
	     data<dataEnd; ++data, ++iter) {
	  bool flagged = (*data == 0);
	  for (unsigned int i=0; i<npol; ++i, ++iter) {
	    *iter = flagged;
	  }
	}
      }
    } break;
    case 2:
    case 3:
    {
      switch (itsParent->getNrBytesPerNrValidSamples()) {

      case 1:
      {
	const unsigned char* data    = itsParent->getNSample1 (rownr, false);
	const unsigned char* dataEnd = data + itsParent->nchan();
	
	if (dataPtr->contiguousStorage()) {
	  for (Array<bool>::contiter iter=dataPtr->cbegin(); data<dataEnd; ++data) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	} else {
	  for (Array<bool>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	}
      } break;

      case 2:
      {
	const unsigned short* data = itsParent->getNSample2 (rownr, false);
	const unsigned short* dataEnd = data + itsParent->nchan();
	
	if (dataPtr->contiguousStorage()) {
	  for (Array<bool>::contiter iter=dataPtr->cbegin(); data<dataEnd; ++data) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	} else {
	  for (Array<bool>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	}
      } break;
	
      case 4:
      {
	const unsigned int* data = itsParent->getNSample4 (rownr, false);
	const unsigned int* dataEnd = data + itsParent->nchan();

	if (dataPtr->contiguousStorage()) {
	  for (Array<bool>::contiter iter=dataPtr->cbegin(); data<dataEnd; ++data) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	} else {
	  for (Array<bool>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    bool flagged = (*data == 0);
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = flagged;
	    }
	  }
	}
      } break;

      default:
	throw;
      }
    } break;
    default:
      throw;
    }
  }


  WeightColumn::~WeightColumn()
  {}
  IPosition WeightColumn::shape (unsigned int)
  {
    return IPosition(1, itsParent->npol());
  }
  void WeightColumn::getArrayfloatV (unsigned int, Array<float>* dataPtr)
  {
    *dataPtr = float(1);
  }

  SigmaColumn::~SigmaColumn()
  {}
  IPosition SigmaColumn::shape (unsigned int)
  {
    return IPosition(1, itsParent->npol());
  }
  void SigmaColumn::getArrayfloatV (unsigned int, Array<float>* dataPtr)
  {
    *dataPtr = float(1);
  }

  WSpectrumColumn::~WSpectrumColumn()
  {}
  IPosition WSpectrumColumn::shape (unsigned int)
  {
    return IPosition(2, itsParent->npol(), itsParent->nchan());
  }
  void WSpectrumColumn::getArrayfloatV (unsigned int rownr, Array<float>* dataPtr)
  {
    double maxn = itsParent->maxnSample();
    unsigned int npol = itsParent->npol();

    switch (itsParent->getLofarStManVersion()) {
    case 1:
    {
      const unsigned short* data    = itsParent->getNSample2 (rownr, true);
      const unsigned short* dataEnd = data + itsParent->nchan();

      if (dataPtr->contiguousStorage()) {
	for (Array<float>::contiter iter=dataPtr->cbegin();
	     data<dataEnd; ++data) {
	  float weight = *data / maxn;
	  for (unsigned int i=0; i<npol; ++i, ++iter) {
	    *iter = weight;
	  }
	}
      } else {
	for (Array<float>::iterator iter=dataPtr->begin();
	     data<dataEnd; ++data, ++iter) {
	  float weight = *data / maxn;
	  for (unsigned int i=0; i<npol; ++i, ++iter) {
	    *iter = weight;
	  }
	}
      }
    } break;
    case 2:
    case 3:
    {
      switch (itsParent->getNrBytesPerNrValidSamples()) {
      case 1:
      {
	const unsigned char* data    = itsParent->getNSample1(rownr, true);
	const unsigned char* dataEnd = data + itsParent->nchan();
	
	if (dataPtr->contiguousStorage()) {
	  for (Array<float>::contiter iter=dataPtr->cbegin();
	       data<dataEnd; ++data) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	} else {
	  for (Array<float>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	}
      } break;
      case 2:
      {
	const unsigned short* data    = itsParent->getNSample2(rownr, true);
	const unsigned short* dataEnd = data + itsParent->nchan();
	
	if (dataPtr->contiguousStorage()) {
	  for (Array<float>::contiter iter=dataPtr->cbegin();
	       data<dataEnd; ++data) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	} else {
	  for (Array<float>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	}
      } break;

      case 4:
      {
	const unsigned int* data    = itsParent->getNSample4(rownr, true);
	const unsigned int* dataEnd = data + itsParent->nchan();
	
	if (dataPtr->contiguousStorage()) {
	  for (Array<float>::contiter iter=dataPtr->cbegin();
	       data<dataEnd; ++data) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	} else {
	  for (Array<float>::iterator iter=dataPtr->begin();
	       data<dataEnd; ++data, ++iter) {
	    float weight = *data / maxn;
	    for (unsigned int i=0; i<npol; ++i, ++iter) {
	      *iter = weight;
	    }
	  }
	}
      } break;
      default:
      	throw;
      }
    } break;
    default:
      throw;
    }
  }

  FlagCatColumn::~FlagCatColumn()
  {}
  bool FlagCatColumn::isShapeDefined (unsigned int)
  {
    return false;
  }
  IPosition FlagCatColumn::shape (unsigned int)
  {
    throw DataManError ("LofarStMan: no data in column FLAG_CATEGORY");
  }

} //# end namespace
