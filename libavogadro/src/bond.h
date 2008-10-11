/**********************************************************************
  Bond - Bond class derived from the base Primitive class

  Copyright (C) 2007 Donald Ephraim Curtis
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

#ifndef BOND_H
#define BOND_H

#include <avogadro/primitive.h>
#include <avogadro/atom.h>

namespace OpenBabel {
  class OBBond;
}

namespace Avogadro {

  /**
   * @class Bond primitive.h <avogadro/primitive.h>
   * @brief Bond Class
   * @author Donald Ephraim Curtis
   *
   * The Bond class is a Primitive subclass that provides a wrapper around
   * OpenBabel::OBBond.  This class is provided to give more control of
   * the OpenBabel::OBBond class through slots/signals provided by the
   * Primitive superclass.
   */
  class Atom;
  class BondPrivate;
  class A_EXPORT Bond : public Primitive
  {
    Q_OBJECT

    public:
      /**
       * Constructor
       *
       * @param parent the object parent.
       */
      Bond(QObject *parent=0);

      inline unsigned int index() const // replaces GetIdx
      {
        return m_index;
      }
      inline void setIndex(unsigned int index)
      {
        m_index = index;
      }
      /// FIXME: More functions that need fixing up!
      inline unsigned long int beginAtomId() const { return m_beginAtomId; }
      void setBegin(Atom* atom) { m_beginAtomId = atom->id(); }
      inline unsigned long int endAtomId() const { return m_endAtomId; }
      void setEnd(Atom* atom) { m_endAtomId = atom->id(); }
      /// The order of the bond - 1 = single, 2 = double etc
      inline int order() const { return m_order; }
      inline void setOrder(int order) { m_order = order; }
      inline double length() const { return 0.0; }

      bool setOBBond(OpenBabel::OBBond *obbond);

      Bond& operator=(const Bond& other);

    private:
      unsigned int m_index;
      unsigned long int m_beginAtomId, m_endAtomId;
      int m_order;
      /* shared d_ptr with Primitive */
      Q_DECLARE_PRIVATE(Bond)
  };

} // End namespace Avogadro

#endif
