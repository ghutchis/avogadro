/**********************************************************************
  Fragment - Fragment class derived from the base Primitive class

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

#ifndef FRAGMENT_H
#define FRAGMENT_H

#include <avogadro/primitive.h>

#include <QList>

namespace Avogadro {

  /**
   * @class Fragment fragment.h <avogadro/fragment.h>
   * @brief Fragment Class
   * @author Marcus D. Hanwell
   *
   * The Fragment class is a Primitive subclass that provides a generic way
   * of addressing fragments. This is intended to be suitable for rings,
   * residues, molecule fragments etc. That is anything that needs to address
   * a subset of atoms/bonds in a Molecule.
   */
  class FragmentPrivate;
  class A_EXPORT Fragment : public Primitive
  {
    Q_OBJECT

    public:
      /**
       * Constructor
       *
       * @param parent the object parent.
       */
      Fragment(QObject *parent=0);

      explicit Fragment(Type type, QObject *parent=0);

      ~Fragment();

      /**
       * @return the name of the fragment.
       * @note Replaces GetName().
       */
      inline QString name() { return m_name; }

      /**
       * Set the name of the fragment.
       */
      inline void setName(QString name) { m_name = name; }

      void addAtom(unsigned long int id);
      void removeAtom(unsigned long int id);
      QList<unsigned long int> atoms();
      void addBond(unsigned long int id);
      void removeBond(unsigned long int id);
      QList<unsigned long int> bonds();

    protected:
      QString m_name;
      QList<unsigned long int> m_atoms;
      QList<unsigned long int> m_bonds;

    private:
      Q_DECLARE_PRIVATE(Fragment)
  };

} // End namespace Avoagdro

#endif
