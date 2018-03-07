/*
  LFA Lab - Library to simplify local Fourier analysis.
  Copyright (C) 2018  Hannah Rittich

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
*/

#include "SystemSymbolProperties.h"
#include "MathUtil.h"

namespace lfa {

SystemSymbolProperties::SystemSymbolProperties(
  int rows,
  int cols,
  FoProperties element_properties)
  : m_rows(rows),
    m_cols(cols),
    m_element_properties(element_properties)
{

}

SystemSymbolProperties
  SystemSymbolProperties::operator* (const SystemSymbolProperties& other)
  const
{
  if (m_cols != other.m_rows)
    throw logic_error("Columns and rows of the systems mismatch.");

  return SystemSymbolProperties(m_rows,
                                other.m_cols,
                                m_element_properties
                                  * other.m_element_properties);
}

SystemSymbolProperties
  SystemSymbolProperties::operator+ (const SystemSymbolProperties& other)
  const
{
  if (m_cols != other.m_cols || m_rows != other.m_rows)
    throw logic_error("Rows or columns mismatch.");

  return SystemSymbolProperties(m_rows,
                                m_cols,
                                m_element_properties
                                  + other.m_element_properties);
}

SystemSymbolProperties
  operator* (double scalar, const SystemSymbolProperties& other)
{
  return other;
}

SystemSymbolProperties SystemSymbolProperties::inverse() const
{
  if (m_rows != m_cols)
    throw logic_error("Cannot invert non-square system.");

  if (m_element_properties.input() != m_element_properties.output())
    throw logic_error("Cannot invert system with non-square entries.");

  return SystemSymbolProperties(m_cols,
                                m_rows,
                                m_element_properties.inverse());
}

SystemSymbolProperties SystemSymbolProperties::adjoint() const
{
  return SystemSymbolProperties(m_cols,
                                m_rows,
                                m_element_properties.adjoint());
}

SystemSymbolProperties properties_of_symbol_system(
  MatrixContainer<FoProperties> properties)
{
  if (properties.cols() <= 0 || properties.rows() <= 0) {
    throw std::logic_error("System has zero rows or columns, which is not allowed.");
  }

  // ensure that all operators have the same outputs and same inputs
  for (int i=0; i < properties.rows(); ++i) {
    for (int j=0; j < properties.cols(); ++j) {
      if (properties(0,0).input() != properties(i,j).input()) {
        throw logic_error("Inputs are not homogeneous.");
      }
      if (properties(0,0).output() != properties(i,j).output()) {
        throw logic_error("Outputs are not homogeneous.");
      }
    }
  }

  return SystemSymbolProperties(properties.rows(),
                                properties.cols(),
                                properties(0,0));
}

}

