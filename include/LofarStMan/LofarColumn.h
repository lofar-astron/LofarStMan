//# LofarColumn.h: A Column in the LOFAR Storage Manager
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

#ifndef LOFAR_LOFARSTMAN_LOFARCOLUMN_H
#define LOFAR_LOFARSTMAN_LOFARCOLUMN_H


//# Includes
#include "LofarStMan/LofarStMan.h"

#include <casacore/tables/DataMan/StManColumn.h>
#include <casacore/measures/Measures/MeasFrame.h>
#include <casacore/measures/Measures/MDirection.h>
#include <casacore/measures/Measures/MBaseline.h>
#include <casacore/casa/Arrays/Vector.h>
#include <casacore/casa/Arrays/IPosition.h>
#include <casacore/casa/Containers/Block.h>
#include <casacore/casa/OS/Conversion.h>

#include <vector>

namespace LOFAR {

// <summary>
// A column in the LOFAR Storage Manager.
// </summary>

// <use visibility=local>

// <reviewed reviewer="UNKNOWN" date="before2004/08/25" tests="tLofarStMan.cc">
// </reviewed>

// <prerequisite>
//# Classes you should understand before using this one.
//   <li> <linkto class=LofarStMan>LofarStMan</linkto>
// </prerequisite>

// <synopsis>
// For each column a specific Column class exists.
// </synopsis>

class LofarColumn : public casacore::StManColumn
{
public:
  explicit LofarColumn (LofarStMan* parent, int dtype)
    : StManColumn (dtype),
      itsParent   (parent)
  {}
  ~LofarColumn() override;
  // Most columns are not writable (only DATA is writable).
  bool isWritable() const override;
  // Set column shape of fixed shape columns; it does nothing.
  void setShapeColumn (const casacore::IPosition& shape) final override;
  // Prepare the column. By default it does nothing.
  virtual void prepareCol();
protected:
  LofarStMan* itsParent;
};

// <summary>ANTENNA1 column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class Ant1Column final : public LofarColumn
{
public:
  explicit Ant1Column (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  ~Ant1Column() override;
  void getIntV (unsigned int rowNr, int* dataPtr) override;
};

// <summary>ANTENNA2 column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class Ant2Column final : public LofarColumn
{
public:
  explicit Ant2Column (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  ~Ant2Column() override;
  void getIntV (unsigned int rowNr, int* dataPtr) override;
};

// <summary>TIME and TIME_CENTROID column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class TimeColumn final : public LofarColumn
{
public:
  explicit TimeColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  ~TimeColumn() override;
  void getdoubleV (unsigned int rowNr, double* dataPtr) override;
private:
  double itsValue;
};

// <summary>INTERVAL and EXPOSURE column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class IntervalColumn final : public LofarColumn
{
public:
  explicit IntervalColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  ~IntervalColumn() override;
  void getdoubleV (unsigned int rowNr, double* dataPtr) override;
private:
  double itsValue;
};

// <summary>All columns in the LOFAR Storage Manager with value 0.</summary>
// <use visibility=local>
class ZeroColumn : public LofarColumn
{
public:
  explicit ZeroColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~ZeroColumn();
  virtual void getIntV (unsigned int rowNr, int* dataPtr);
private:
  int itsValue;
};

// <summary>All columns in the LOFAR Storage Manager with value false.</summary>
// <use visibility=local>
class FalseColumn : public LofarColumn
{
public:
  explicit FalseColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~FalseColumn();
  virtual void getBoolV (unsigned int rowNr, bool* dataPtr);
private:
  bool itsValue;
};

// <summary>UVW column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class UvwColumn : public LofarColumn
{
public:
  explicit UvwColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~UvwColumn();
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArraydoubleV (unsigned int rowNr,
                                casacore::Array<double>* dataPtr);
  virtual void prepareCol();
private:
  casacore::MDirection              itsPhaseDir;    //# could be SUN, etc.
  casacore::MDirection              itsJ2000Dir;    //# Phase dir in J2000
  casacore::MeasFrame               itsFrame;
  std::vector<casacore::MBaseline>       itsAntMB;
  std::vector<casacore::Vector<double> > itsAntUvw;
  casacore::Block<bool>             itsUvwFilled;
  int                           itsLastBlNr;
  bool                          itsCanCalc;     //# false = UVW cannot be calc.
};

// <summary>DATA column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class DataColumn : public LofarColumn
{
public:
  explicit DataColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~DataColumn();
  virtual bool isWritable() const;
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArrayComplexV (unsigned int rowNr,
                                 casacore::Array<casacore::Complex>* dataPtr);
  virtual void putArrayComplexV (unsigned int rowNr,
                                 const casacore::Array<casacore::Complex>* dataPtr);
};

// <summary>FLAG column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class FlagColumn : public LofarColumn
{
public:
  explicit FlagColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~FlagColumn();
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArrayBoolV (unsigned int rowNr,
                              casacore::Array<bool>* dataPtr);
};

// <summary>WEIGHT column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class WeightColumn : public LofarColumn
{
public:
  explicit WeightColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~WeightColumn();
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArrayfloatV (unsigned int rowNr,
                               casacore::Array<float>* dataPtr);
};

// <summary>SIGMA column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class SigmaColumn : public LofarColumn
{
public:
  explicit SigmaColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~SigmaColumn();
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArrayfloatV (unsigned int rowNr,
                               casacore::Array<float>* dataPtr);
};

// <summary>WEIGHT_SPECTRUM column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class WSpectrumColumn : public LofarColumn
{
public:
  explicit WSpectrumColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~WSpectrumColumn();
  virtual casacore::IPosition shape (unsigned int rownr);
  virtual void getArrayfloatV (unsigned int rowNr,
                               casacore::Array<float>* dataPtr);
};

// <summary>FLAG_CATEGORY column in the LOFAR Storage Manager.</summary>
// <use visibility=local>
class FlagCatColumn : public LofarColumn
{
public:
  explicit FlagCatColumn (LofarStMan* parent, int dtype)
    : LofarColumn(parent, dtype) {}
  virtual ~FlagCatColumn();
  virtual bool isShapeDefined (unsigned int rownr);
  virtual casacore::IPosition shape (unsigned int rownr);
};


} //# end namespace

#endif
