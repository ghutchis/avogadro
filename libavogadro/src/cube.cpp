/**********************************************************************
  Cube - Primitive class to encapsulate volumetric data

  Copyright (C) 2008 Marcus D. Hanwell

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.sourceforge.net/>

  Avogadro is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  Avogadro is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program; if not, write to the Free Software
  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
  02110-1301, USA.
 **********************************************************************/

#include "cube.h"

#include <QDebug>

namespace Avogadro {

  Cube::Cube() : m_data(0), m_min(0.0, 0.0, 0.0), m_max(0.0, 0.0, 0.0),
    m_spacing(0.0, 0.0, 0.0), m_points(0, 0, 0)
  {
  }

  Cube::~Cube()
  {
  }

  bool Cube::setLimits(const Eigen::Vector3d &min, const Eigen::Vector3d &max,
                       const Eigen::Vector3i &points)
  {
    // We can calculate all necessary properties and initialise our data
    Eigen::Vector3d delta = max - min;
    m_spacing = Eigen::Vector3d(delta.x() / points.x(),
                                delta.y() / points.y(),
                                delta.z() / points.z());
    m_min = min;
    m_max = max;
    m_points = points;
    m_data.resize(m_points.x() * m_points.y() * m_points.z());
    return true;
  }

  bool Cube::setLimits(const Eigen::Vector3d &min, const Eigen::Vector3d &max,
                       double spacing)
  {
    Eigen::Vector3i points;
    Eigen::Vector3d delta = max - min;
    delta = delta / spacing;
    points = Eigen::Vector3i(delta.x(), delta.y(), delta.z());
    return setLimits(min, max, points);
  }

  bool Cube::setLimits(const Molecule &mol, double spacing, double padding)
  {
    QList<Atom *> atoms = mol.atoms();
    Eigen::Vector3d min(0.0, 0.0, 0.0);
    Eigen::Vector3d max(0.0, 0.0, 0.0);
    foreach (Atom *atom, atoms) {
      if (atom->pos().x() < min.x())
        min.x() = atom->pos().x();
      else if (atom->pos().x() > max.x())
        max.x() = atom->pos().x();
      if (atom->pos().y() < min.y())
        min.y() = atom->pos().y();
      else if (atom->pos().y() > max.y())
        max.y() = atom->pos().y();
      if (atom->pos().z() < min.z())
        min.z() = atom->pos().z();
      else if (atom->pos().z() > max.z())
        max.z() = atom->pos().z();
    }
    min = min - Eigen::Vector3d(padding, padding, padding);
    max = max + Eigen::Vector3d(padding, padding, padding);
    return setLimits(min, max, spacing);
  }

  std::vector<double> Cube::data()
  {
    return m_data;
  }

  bool Cube::setData(const std::vector<double> &values)
  {
    if (values.size() == m_points.x() * m_points.y() * m_points.z()) {
      m_data = values;
      qDebug() << "Loaded in cube data" << m_data.size();
      return true;
    }
    else
      return false;
  }

  int Cube::index(const Eigen::Vector3d &pos)
  {
    int x, y, z;
    // Calculate how many steps each coordinate is along its axis
    x = (pos.x() - m_min.x()) / m_spacing.x();
    y = (pos.y() - m_min.y()) / m_spacing.y();
    z = (pos.z() - m_min.z()) / m_spacing.z();
    return x*m_points.y()*m_points.z() + y*m_points.z() + z;
  }

  Eigen::Vector3d Cube::position(int index)
  {
    int x, y, z;
    x = static_cast<int>(index / (m_points.y()*m_points.z()));
    y = static_cast<int>((index - (x*m_points.y()*m_points.z())) / m_points.z());
    z = index % m_points.z();
    return Eigen::Vector3d(x * m_spacing.x() + m_min.x(),
                           y * m_spacing.y() + m_min.y(),
                           z * m_spacing.z() + m_min.z());
  }

  double Cube::value(int i, int j, int k) const
  {
    int index = i*m_points.y()*m_points.z() + j*m_points.z() + k;
    if (index < m_data.size())
      return m_data.at(index);
    else
      return 0.0;
  }

  double Cube::value(const Eigen::Vector3d &) const
  {
    // This is a really expensive operation and so should be avoided
    // Interpolate the value at the supplied vector
  }

  bool Cube::setValue(int i, int j, int k, double value)
  {
    int index = i*m_points.y()*m_points.z() + j*m_points.z() + k;
    if (index < m_data.size()) {
      m_data[index] = value;
      return true;
    }
    else
      return false;
  }

} // End namespace Avogadro

#include "cube.moc"
