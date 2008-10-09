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

 #ifndef CUBE_H
 #define CUBE_H

#include <avogadro/primitive.h>

#include <vector>

namespace Avogadro {

  class A_EXPORT Cube : public Primitive
  {
  Q_OBJECT

  public:
    Cube();
    ~Cube();

   /**
    * @return The minimum point in the cube.
    */
    Eigen::Vector3d min() { return m_min; }

   /**
    * @return The maximum point in the cube.
    */
    Eigen::Vector3d max() { return m_max; }

   /**
    * @return The spacing of the grid.
    */
    Eigen::Vector3d spacing() { return m_spacing; }

   /**
    * @return The x, y and z dimensions of the cube.
    */
    Eigen::Vector3i dimensions() { return m_points; }

    /**
     * Set the limits of the cube.
     * @param min The minimum point in the cube.
     * @param max The maximum point in the cube.
     * @param points The number of (integer) points in the cube.
     */
    bool setLimits(const Eigen::Vector3d &min, const Eigen::Vector3d &max,
                   const Eigen::Vector3i &points);
    bool setLimits(const Eigen::Vector3d &min, const Eigen::Vector3d &max,
                   double spacing);
    bool setLimits(const Molecule &mol, double spacing, double padding);

    /**
     * @return Vector containing all the data in a one-dimensional array.
     */
    std::vector<double> data();

    /**
     * Set the values in the cube to those passed in the vector.
     */
    bool setData(const std::vector<double> &values);

    /**
     * @return Index of the point closest to the position supplied.
     * @param pos Position to get closest index for.
     */
    int index(const Eigen::Vector3d &pos);

    /**
     * @param index Index to be translated to a position.
     * @return Position of the given index.
     */
    Eigen::Vector3d position(int index);

    double value(int i, int j, int k) const;
    bool setValue(int i, int j, int k, double value);

    inline void setName(QString name) { m_name = name; }
    inline QString name() { return m_name; }

  private:
    std::vector<double> m_data;
    Eigen::Vector3d m_min, m_max, m_spacing;
    Eigen::Vector3i m_points;
    QString m_name;
  };
} // End namespace Avogadro

 #endif
 