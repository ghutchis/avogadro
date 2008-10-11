/**********************************************************************
  Molecule - Molecule class derived from the base Primitive class

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

#ifndef MOLECULE_H
#define MOLECULE_H

#include <avogadro/primitive.h>

namespace OpenBabel {
  class OBAtom;
  class OBBond;
  class OBMol;
}

namespace Avogadro {

  // Declare new classes
  class Atom;
  class Bond;
  class Residue;
  class Cube;

  /**
   * @class Molecule primitive.h <avogadro/primitive.h>
   * @brief Molecule Class
   * @author Donald Ephraim Curtis
   *
   * The Molecule class implements the OpenBabel::OBMol virtual functions
   * in order to not only use our primitive objects but also to provide signals
   * based on internal OpenBabel actions.  In terms of a Model-View architecture,
   * this is our model class and is used by our various views to hold
   * all required data.
   */
  class MoleculePrivate;
  class A_EXPORT Molecule : public Primitive
  {
    Q_OBJECT

    public:
      /**
       * Constructor
       *
       * @param parent the object parent.
       */
      Molecule(QObject *parent=0);
      Molecule(const Molecule &other);
      virtual ~Molecule();
      void update();

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Creates a new Atom object.
       *
       * @return pointer to a newly allocated Atom object
       */
//      Atom *CreateAtom();
//      void createAtom(Atom *) { ; }

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Creates a new Bond object.
       *
       * @return pointer to a newly allocated Bond object
       */
//      Bond *CreateBond();

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Creates a new Residue object.
       *
       * @return pointer to a newly allocated Residue object
       */
//      Residue *CreateResidue();

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Deletes an Atom object.
       *
       * @param atom the atom to delete
       */
//      void DestroyAtom(Atom* atom);

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Deletes an Bond object.
       *
       * @param bond the bond to delete
       */
//      void DestroyBond(Bond* bond);

      /**
       * Virtual function inherited from OpenBabel::OBMol.
       * Deletes an Residue object.
       *
       * @param residue the residue to delete
       */
//      void DestroyResidue(Residue* residue);

      /**
       * These are new functions to replace OBMol::NewAtom
       * to use our new unique identifier functionality.
       */
      Atom *newAtom();
      Atom *newAtom(unsigned long id);
      void deleteAtom(Atom *atom);
      void deleteAtom(unsigned long int id);

      /**
       * These are new functions to replace OBMol::NewBond
       * to use our new unique identifier functionality.
       */
      Bond *newBond();
      Bond *newBond(unsigned long id);
      void deleteBond(Bond *bond);
      void deleteBond(unsigned long int id);

      Atom *getAtomById(unsigned long id) const;

      Bond *getBondById(unsigned long id) const;

      void addHydrogens(Atom *atom = 0);
      void deleteHydrogens(Atom *atom);
      void deleteHydrogens();

      Cube *newCube();
      void deleteCube(Cube *cube);

      /// FIXME These are new functions that need fleshing out to work properly
      unsigned int numAtoms() const;
      unsigned int numBonds() const;
      unsigned int numResidues() const;

      Atom* atom(int index); // Replaces GetAtom
      Bond* bond(int index); // Replaces GetBond

      Bond* bond(unsigned long int id1, unsigned long int id2);
      Bond* bond(const Atom*, const Atom*);

      QList<Atom *> atoms() const;
      QList<Bond *> bonds() const;
      QList<Cube *> cubes() const;

      /** FIXME Implement me!
       * Delete all elements of the molecule
       */
       void clear();

      /**
       * Get access to an OpenBabel atom, this is a copy of the internal data
       * structure in OpenBabel form, you must call setOBAtom in order to save
       * any changes you make to this object.
       */
      OpenBabel::OBMol OBMol();
      bool setOBMol(OpenBabel::OBMol *obmol);

      const Eigen::Vector3d & center() const;
      const Eigen::Vector3d & normalVector() const;
      const double & radius() const;
      const Atom *farthestAtom() const;
      void translate(const Eigen::Vector3d&) { ; }

      Molecule& operator=(const Molecule& other);

      Molecule& operator+=(const Molecule& other);

    protected:
      MoleculePrivate * const d_ptr;

    private:
      /* shared d_ptr with Primitive */
      Q_DECLARE_PRIVATE(Molecule)

      void computeGeomInfo() const;

    private Q_SLOTS:
      /**
       * Function which handles when a child primitive has been
       * updated.  The response is to find the sender object
       * and then emit a signal passing the sender as a parameter.
       *
       * @sa primitiveAdded
       * @sa primitiveUpdated
       * @sa primitiveRemoved
       */
      void updatePrimitive();

    Q_SIGNALS:
      /**
       * Emitted when a child primitive is added.
       *
       * @param primitive pointer to the primitive that was added
       */
      void primitiveAdded(Primitive *primitive);
      /**
       * Emitted when a child primitive is updated.
       *
       * @param primitive pointer to the primitive that was updated
       */
      void primitiveUpdated(Primitive *primitive);
      /**
       * Emitted when a child primitive is deleted.
       *
       * @param primitive pointer to the primitive that was updated before it is free'd
       */
      void primitiveRemoved(Primitive *primitive);
  };

} // End namespace Avogadro

#endif
