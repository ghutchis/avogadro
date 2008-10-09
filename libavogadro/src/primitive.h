/**********************************************************************
  Primitive - Wrapper class around the OpenBabel classes

  Copyright (C) 2007 Donald Ephraim Curtis

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

#ifndef PRIMITIVE_H
#define PRIMITIVE_H

#include <avogadro/global.h>

#include <QAbstractItemModel>

#include <Eigen/Core>

// Forward declarations
class QReadWriteLock;

namespace OpenBabel {
  class OBAtom;
  class OBBond;
  class OBMol;
}

namespace Avogadro {

  // Declare new classes
  class Atom;
  class Bond;
  class Molecule;
  class Cube;

  /**
   * @class Primitive primitive.h <avogadro/primitive.h>
   * @brief Base class for all primitives (Molecule, Atom, Bond, Residue, ...).
   */

  class PrimitivePrivate;
  class A_EXPORT Primitive : public QObject
  {
    Q_OBJECT
    Q_PROPERTY(Type type READ type)
    Q_ENUMS(Type)

    public:
      /**
       * This enum allows us to iterate through the various types
       * of primitives.
       */
      enum Type {
        /// Untyped Primitive
        OtherType=0,
        /// Molecule Primitive
        MoleculeType,
        /// Atom Primitive
        AtomType,
        /// Bond Primitive
        BondType,
        /// Residue Primitive
        ResidueType,
        /// Chain Primitive (i.e., a set of residues)
        ChainType,
        /// Surface Primitive
        SurfaceType,
        /// Plane Primitive
        PlaneType,
        /// Grid Primitive
        GridType,
        /// Points (i.e., non-atoms)
        PointType,
        /// Vectors (i.e., arrows, dipole moments)
        VectorType,
        /// Non-bonded interactions (i.e., non-bond connections)
        NonbondedType,
        /// Text annoations
        TextType,
        /// End Placeholder
        LastType,
        /// First Placeholder
        FirstType=OtherType
      };


      /**
       * Default constructor.
       * @param parent the object parent
       */
      Primitive(QObject *parent = 0);
      /**
       * Constructor
       * @param type the primitive type
       * @param parent the object parent
       */
      explicit Primitive(Type type, QObject *parent=0);
      /**
       * Deconstructor
       */
      virtual ~Primitive();

      /**
       * Function used to push changes to a primitive to
       * the rest of the system.  At this time there is no
       * way (other than this) to generate a signal when
       * properties of a primitive change.
       *
       * In the case of the Atom primitive, this should be called
       * when changes to coordinates have been made.
       */
      void update();

      /**
       * @property Type
       * Holds the primitive type
       */

      /**
       * @return the primitive type (one of Primitive::Type)
       */
      Type type() const;

      QReadWriteLock *lock();

      void setId(unsigned long m_id);
      unsigned long id() const;


    Q_SIGNALS:
      /**
       * Emitted when the primitive has been updated.
       */
      void updated();

    protected:
      PrimitivePrivate * const d_ptr;
      Primitive(PrimitivePrivate &dd, QObject *parent = 0);
      Primitive(PrimitivePrivate &dd, Type type, QObject *parent=0);

    private:
      Q_DECLARE_PRIVATE(Primitive)

  };

  /**
   * @class Atom primitive.h <avogadro/primitive.h>
   * @brief Atom Class
   * @author Donald Ephraim Curtis
   *
   * The Atom class is a Primitive subclass that provides a wrapper around
   * OpenBabel::OBAtom.  This class is provided to give more control of
   * the OpenBabel::OBAtom class through slots/signals provided by the
   * Primitive superclass.
   */
  class AtomPrivate;
  class A_EXPORT Atom : public Primitive
  {
    Q_OBJECT

    public:
      /**
       * Constructor
       *
       * @param parent the object parent.
       */
      Atom(QObject *parent=0);

      /** Returns the position of the atom, as a Eigen::Vector3d. This is similar to
        * the OBAtom::GetVector() method, which returns the position as a OpenBabel::vector3.
        *
        * Rationale for inlining: this method only does a cast on the return value of OBAtom::GetVector().
        * The memory layouts of the types between which it casts are not likely to change: both
        * types represent 3D vectors of doubles, and there's only one sane way to represent them:
        * struct{ double x,y,z; }.
        *
        * @return OBAtom::GetVector() but reinterpret_casted as a const Eigen::Vector3d &
        */
      inline const Eigen::Vector3d &pos () const
      {
        return m_pos;
      }

      /** Sets the position of the atom, from a Eigen::Vector3d. This is similar to
        * the OBAtom::SetVector() method, which sets the position from a OpenBabel::vector3.
        *
        * Rationale for inlining: this method only does a cast on the argument of OBAtom::SetVector().
        * The memory layouts of the types between which it casts are not likely to change: both
        * types represent 3D vectors of doubles, and there's only one sane way to represent them:
        * struct{ double x,y,z; }.
        */
      inline void setPos(const Eigen::Vector3d &vec)
      {
        m_pos = vec;
      }

      /// FIXME: I should return a real ID!
      inline unsigned int index() const  // Replaces OpenBabel::OBAtom::GetIdx
      {
        return m_index;
      }
      inline void setIndex(unsigned int index)
      {
        m_index = index;
      }
      /// FIXME These are new functions that need fleshing out to work properly
      inline int atomicNumber() const // Replaces GetAtomicNum
      { return m_atomicNum; }

      inline void setAtomicNumber(int num)
      { m_atomicNum = num; }

      void addBond(Bond* bond);
      void deleteBond(Bond* bond);
      QList<unsigned long int> bonds() { return m_bonds; }

      QList<unsigned long int> neighbors() { return m_neighbors; }

      double valence() { return 0.0; }

      bool isHydrogen() const { return m_atomicNum == 1; }

      /// FIXME New function that needs implenenting
      inline double partialCharge() const
      { return 0.0; }

      /// Our OpenBabel conversion functions
      OpenBabel::OBAtom OBAtom();
      bool setOBAtom(OpenBabel::OBAtom *obatom);

      Atom& operator=(const Atom& other);

    private:
      /* shared d_ptr with Primitive */
      unsigned int m_index;
      Eigen::Vector3d m_pos;
      int m_atomicNum;
      QList<unsigned long int> m_bonds, m_neighbors;
      Q_DECLARE_PRIVATE(Atom)
  };

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

  /**
   * @class Residue primitive.h <avogadro/primitive.h>
   * @brief Residue Class
   * @author Donald Ephraim Curtis
   *
   * The Residue class is a Primitive subclass that provides a wrapper around
   * OpenBabel::OBResidue.  This class is provided to give more control of
   * the OpenBabel::OBResidue class through slots/signals provided by the
   * Primitive superclass.
   */
  class A_EXPORT Residue : public Primitive
  {
    Q_OBJECT

    public:
      /**
       * Constructor
       *
       * @param parent the object parent.
       */
      Residue(QObject *parent=0): Primitive(ResidueType, parent) { }

      /// FIXME More functions that need to be fleshed out
      QString name() // Replaces GetName
      { return "FIXME"; }

      QString numString() // Replaces GetName
      { return "FIXME"; }
  };

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

} // namespace Avogadro

Q_DECLARE_METATYPE(Avogadro::Primitive*)

#endif // PRIMITIVE_H
