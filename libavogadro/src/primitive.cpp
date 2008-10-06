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

#include <config.h>

#include <avogadro/primitive.h>
#include <Eigen/Regression>
#include <Eigen/Geometry>

#include <QReadWriteLock>

#include <QDebug>

#include <openbabel/mol.h>

namespace Avogadro {

  class PrimitivePrivate {
    public:
      PrimitivePrivate() : type(Primitive::OtherType), id(-1) {};

      enum Primitive::Type type;
      QReadWriteLock lock;

      unsigned long id;
  };

  Primitive::Primitive(QObject *parent) : QObject(parent), d_ptr(new PrimitivePrivate) {}

  Primitive::Primitive(enum Type type, QObject *parent) : QObject(parent), d_ptr(new PrimitivePrivate)
  {
    Q_D(Primitive);
    d->type = type;
  }

  Primitive::Primitive(PrimitivePrivate &dd, QObject *parent) : QObject(parent), d_ptr(&dd) {}

  Primitive::Primitive(PrimitivePrivate &dd, enum Type type, QObject *parent) : QObject(parent), d_ptr(&dd)
  {
    Q_D(Primitive);
    d->type = type;
  }

  Primitive::~Primitive()
  {
    delete d_ptr;
  }

  enum Primitive::Type Primitive::type() const
  {
    Q_D(const Primitive);
    return d->type;
  }

  QReadWriteLock *Primitive::lock()
  {
    Q_D(Primitive);
    return &d->lock;
  }

  void Primitive::update()
  {
    emit updated();
  }

  void Primitive::setId(unsigned long m_id)
  {
    Q_D(Primitive);
    d->id = m_id;
  }

  unsigned long Primitive::id() const
  {
    Q_D(const Primitive);
    return d->id;
  }


  class AtomPrivate : public PrimitivePrivate {
    public:
      AtomPrivate() : PrimitivePrivate() {}

  };

  Atom::Atom(QObject *parent) : Primitive(*new AtomPrivate, AtomType, parent),
    m_index(0), m_pos(0., 0., 0.), m_atomicNum(0)
  {
  }

  void Atom::addBond(Bond* bond)
  {
    m_bonds.push_back(bond->id());
    // Update the neighbors list
    if (bond->beginAtomId() == id())
      m_neighbors.push_back(bond->endAtomId());
    else
      m_neighbors.push_back(bond->beginAtomId());
  }

  void Atom::deleteBond(Bond* bond)
  {
    int index = m_bonds.indexOf(bond->id());
    if (index >= 0)
      m_bonds.removeAt(index);

    // Update the neighbors list too
    if (bond->beginAtomId() == id())
      m_neighbors.removeAt(m_neighbors.indexOf(bond->endAtomId()));
    else
      m_neighbors.removeAt(m_neighbors.indexOf(bond->beginAtomId()));
  }

  OpenBabel::OBAtom Atom::OBAtom()
  {
    // Need to copy all relevant data over to the OBAtom
    OpenBabel::OBAtom obatom;
    obatom.SetVector(m_pos.x(), m_pos.y(), m_pos.z());
    obatom.SetAtomicNum(m_atomicNum);

    return obatom;
  }

  bool Atom::setOBAtom(OpenBabel::OBAtom *obatom)
  {
    // Copy all needed OBAtom data to our atom
    m_pos = Eigen::Vector3d(obatom->x(), obatom->y(), obatom->z());
    m_atomicNum = obatom->GetAtomicNum();
    update();
    return true;
  }

  Atom& Atom::operator=(const Atom& other)
  {
    // Virtually everything here is invariant apart from the index and possibly id
    m_pos = other.m_pos;
    m_atomicNum = other.m_atomicNum;
    m_bonds = other.m_bonds;
    return *this;
  }

  class BondPrivate : public PrimitivePrivate {
    public:
      BondPrivate() : PrimitivePrivate() {}
  };

  Bond::Bond(QObject *parent) : Primitive(*new BondPrivate, BondType, parent),
    m_index(0), m_beginAtomId(0), m_endAtomId(0), m_order(1)
  {
  }

  bool Bond::setOBBond(OpenBabel::OBBond *obbond)
  {
    m_order = obbond->GetBondOrder();
    return true;
  }

  Bond& Bond::operator=(const Bond& other)
  {
    m_beginAtomId = other.m_beginAtomId;
    m_endAtomId = other.m_endAtomId;
    m_order = other.m_order;
    return *this;
  }

  class MoleculePrivate : public PrimitivePrivate {
    public:
      MoleculePrivate() : PrimitivePrivate(), farthestAtom(0), invalidGeomInfo(true), autoId(true), obmol(0) {}
      mutable Eigen::Vector3d       center;
      mutable Eigen::Vector3d       normalVector;
      mutable double                radius;
      mutable Atom *                farthestAtom;
      mutable bool                  invalidGeomInfo;
      bool                          autoId;

      // std::vector used over QVector due to index issues, QVector uses ints
      std::vector<Atom *>           atoms;
      std::vector<Bond *>           bonds;

      // Used to store the index based list (not unique ids)
      QList<Atom *>                 atomList;
      QList<Bond *>                 bondList;

      // Our OpenBabel OBMol object
      OpenBabel::OBMol *            obmol;
  };

  Molecule::Molecule(QObject *parent) : Primitive(*new MoleculePrivate, MoleculeType, parent)
  {
    connect(this, SIGNAL(updated()), this, SLOT(updatePrimitive()));
  }

  Molecule::Molecule(const Molecule &other) : Primitive(*new MoleculePrivate, MoleculeType, other.parent())
  {
    *this = other;
    connect(this, SIGNAL(updated()), this, SLOT(updatePrimitive()));
  }

  Molecule::~Molecule()
  {
    // Need to iterate through all atoms/bonds and destroy them
    Q_D(Molecule);
    foreach (Atom *atom, d->atomList)
      atom->deleteLater();
    foreach (Bond *bond, d->bondList)
      bond->deleteLater();
  }
/*
  Atom * Molecule::CreateAtom()
  {
    Q_D(Molecule);

    d->lock.lockForWrite();
    Atom *atom = new Atom(this);
    connect(atom, SIGNAL(updated()), this, SLOT(updatePrimitive()));

    if(!d->autoId) {
      d->lock.unlock();
      return(atom);
    }

    atom->setId(d->atoms.size());
    d->atoms.push_back(atom);
    d->lock.unlock();

    emit primitiveAdded(atom);
    return(atom);
  }

  Bond * Molecule::CreateBond()
  {
    Q_D(Molecule);

    d->lock.lockForWrite();
    Bond *bond = new Bond(this);
    connect(bond, SIGNAL(updated()), this, SLOT(updatePrimitive()));

    if(!d->autoId) {
      d->lock.unlock();
      return(bond);
    }

    bond->setId(d->bonds.size());
    d->bonds.push_back(bond);
    d->lock.unlock();

    emit primitiveAdded(bond);
    return(bond);
  }

  Residue * Molecule::CreateResidue()
  {
    Residue *residue = new Residue(this);
    connect(residue, SIGNAL(updated()), this, SLOT(updatePrimitive()));
    emit primitiveAdded(residue);
    return(residue);
  }
*/
  Atom *Molecule::newAtom()
  {
    Q_D(Molecule);
    Atom *atom = new Atom;
    d->atoms.push_back(atom);
    d->atomList.push_back(atom);
    atom->setId(d->atoms.size()-1);
    atom->setIndex(d->atomList.size()-1);
    connect(atom, SIGNAL(updated()), this, SLOT(updatePrimitive()));
    emit primitiveAdded(atom);
    return atom;
  }

  // do some fancy footwork when we add an atom previously created
  Atom *Molecule::newAtom(unsigned long id)
  {
    Q_D(Molecule);

    // we have to bypass the emit given by CreateAtom()
    d->autoId = false;
    Atom *atom = new Atom;
    d->autoId = true;

    if(id >= d->atoms.size())
      d->atoms.resize(id+1,0);
    atom->setId(id);
    d->atoms[id] = atom;

    // Does this still want to have the same index as before somehow?
    d->atomList.push_back(atom);
    atom->setIndex(d->atomList.size()-1);

    // now that the id is correct, emit the signal
    connect(atom, SIGNAL(updated()), this, SLOT(updatePrimitive()));
    emit primitiveAdded(atom);
    return(atom);
  }

  void Molecule::deleteAtom(Atom *atom)
  {
    Q_D(Molecule);
    if(atom) {
      // When deleting an atom this also implicitly deletes any bonds to the atom
      QList<unsigned long int> bonds = atom->bonds();
      foreach (unsigned long int bond, bonds)
        deleteBond(getBondById(bond));

      d->atoms[atom->id()] = 0;
      // 1 based arrays stored/shown to user
      int index = atom->index();
      d->atomList.removeAt(index);
      for (int i = index; i < d->atomList.size(); ++i)
        d->atomList[i]->setIndex(i);
      atom->deleteLater();
      disconnect(atom, SIGNAL(updated()), this, SLOT(updatePrimitive()));
      emit primitiveRemoved(atom);
      qDebug() << "Atom" << atom->id() << atom->index() << "deleted";
    }
  }

  void Molecule::deleteAtom(unsigned long int id)
  {
    Q_D(Molecule);
    if (id < d->atoms.size())
      deleteAtom(d->atoms[id]);
  }

  Atom *Molecule::getAtomById(unsigned long id) const
  {
    Q_D(const Molecule);
    if(id < d->atoms.size())
      return d->atoms[id];
    else
      return 0;
  }

  Atom *Molecule::atom(int index)
  {
    Q_D(Molecule);
    if (index >= 0 && index < d->atomList.size())
      return d->atomList[index];
    else
      return 0;
  }

  Bond *Molecule::newBond()
  {
    Q_D(Molecule);
    Bond *bond = new Bond;
    d->bonds.push_back(bond);
    d->bondList.push_back(bond);
    bond->setId(d->bonds.size()-1);
    bond->setIndex(d->bondList.size()-1);
    connect(bond, SIGNAL(updated()), this, SLOT(updatePrimitive()));
    emit primitiveAdded(bond);
    return bond;
  }

  Bond *Molecule::newBond(unsigned long id)
  {
    Q_D(Molecule);

    d->autoId = false;
    Bond *bond = new Bond;
    d->autoId = true;

    if(id >= d->bonds.size())
      d->bonds.resize(id+1,0);
    bond->setId(id);
    d->bonds[id] = bond;

    d->bondList.push_back(bond);
    bond->setIndex(d->bondList.size()-1);

    // now that the id is correct, emit the signal
    connect(bond, SIGNAL(updated()), this, SLOT(updatePrimitive()));
    emit primitiveAdded(bond);
    return(bond);
  }

  void Molecule::deleteBond(Bond *bond)
  {
    Q_D(Molecule);
    if(bond) {
      d->bonds[bond->id()] = 0;
      // 1 based arrays stored/shown to user
      int index = bond->index();
      d->bondList.removeAt(index);
      for (int i = index; i < d->bondList.size(); ++i)
        d->bondList[i]->setIndex(i);
      // Also delete the bond from the attached atoms
      (getAtomById(bond->beginAtomId()))->deleteBond(bond);
      (getAtomById(bond->endAtomId()))->deleteBond(bond);

      bond->deleteLater();
      disconnect(bond, SIGNAL(updated()), this, SLOT(updatePrimitive()));
      emit primitiveRemoved(bond);
      qDebug() << "Bond" << bond->id() << bond->index() << "deleted";
    }
  }

  void Molecule::deleteBond(unsigned long int id)
  {
    Q_D(Molecule);
    if (id < d->bonds.size())
      deleteBond(d->bonds[id]);
  }

  Bond *Molecule::getBondById(unsigned long id) const
  {
    Q_D(const Molecule);
    if(id < d->bonds.size())
    {
      return d->bonds[id];
    }
    return 0;
  }

  Bond *Molecule::bond(int index)
  {
    Q_D(Molecule);
    if (index >= 0 && index < d->bondList.size())
      return d->bondList[index];
    else
      return 0;
  }

  void Molecule::addHydrogens(Atom *atom)
  {
    // Construct an OBMol, call AddHydrogens and translate the changes
    OpenBabel::OBMol obmol = OBMol();
    if (atom)
      obmol.AddHydrogens(obmol.GetAtom(atom->index()+1));
    else
      obmol.AddHydrogens();
    // All new atoms in the OBMol must be the additional hydrogens
    for (unsigned int i = numAtoms()+1; i <= obmol.NumAtoms(); ++i) {
      if (obmol.GetAtom(i)->IsHydrogen()) {
        OpenBabel::OBAtom *obatom = obmol.GetAtom(i);
        Atom *atom = newAtom();
        atom->setOBAtom(obatom);
        // Get the neighbor atom
        OpenBabel::OBBondIterator iter;
        OpenBabel::OBAtom *next = obatom->BeginNbrAtom(iter);
        Bond *bond = newBond();
        bond->setEnd(Molecule::atom(atom->index()));
        bond->setBegin(Molecule::atom(next->GetIdx()-1));
        Molecule::atom(next->GetIdx()-1)->addBond(bond);
        atom->addBond(bond);
      }
    }
  }
  
  void Molecule::deleteHydrogens(Atom *atom)
  {
    // Delete any connected hydrogen atoms
    QList<unsigned long int> neighbors = atom->neighbors();
    foreach (unsigned long int a, neighbors)
      if (getAtomById(a)->isHydrogen())
        deleteAtom(a);
  }

  void Molecule::deleteHydrogens()
  {
    Q_D(Molecule);
    foreach (Atom *atom, d->atomList)
      if (atom->isHydrogen())
        deleteAtom(atom);
  }

  unsigned int Molecule::numAtoms() const
  {
    Q_D(const Molecule);
    return d->atomList.size();
  }

  unsigned int Molecule::numBonds() const
  {
    Q_D(const Molecule);
    return d->bondList.size();
  }

  unsigned int Molecule::numResidues() const
  {
    // FIXME Need to add residues
    return 0;
  }

  void Molecule::updatePrimitive()
  {
    Q_D(Molecule);
    Primitive *primitive = qobject_cast<Primitive *>(sender());
    d->invalidGeomInfo = true;
    emit primitiveUpdated(primitive);
  }

  void Molecule::update()
  {
    Q_D(Molecule);
    d->invalidGeomInfo = true;
    emit updated();
  }

  Bond* Molecule::bond(unsigned long int id1, unsigned long int id2)
  {
    // Take two atom IDs and see if we have a bond between the two
    Q_D(Molecule);
    foreach (Bond *bond, d->bondList) {
      if (bond->beginAtomId() == id1 && bond->endAtomId() == id2)
        return bond;
      if (bond->beginAtomId() == id2 && bond->endAtomId() == id1)
        return bond;
    }
    return 0;
  }

  Bond* Molecule::bond(const Atom *a1, const Atom *a2)
  {
    if (a1 && a2)
      return bond(a1->id(), a2->id());
    else
      return 0;
  }

  QList<Atom *> Molecule::atoms()
  {
    Q_D(Molecule);
    return d->atomList;
  }

  QList<Bond *> Molecule::bonds()
  {
    // Make a QList containing all current atoms
    Q_D(Molecule);
    return d->bondList;
  }

  OpenBabel::OBMol Molecule::OBMol()
  {
    Q_D(Molecule);
    // Right now we make an OBMol each time
    OpenBabel::OBMol obmol;
    obmol.BeginModify();
    foreach (Atom *atom, d->atomList) {
      OpenBabel::OBAtom *a = obmol.NewAtom();
      OpenBabel::OBAtom obatom = atom->OBAtom();
      *a = obatom;
//      qDebug() << "Atoms" << obmol.NumAtoms();
    }
    foreach (Bond *bond, d->bondList) {
      obmol.AddBond(getAtomById(bond->beginAtomId())->index() + 1,
                    getAtomById(bond->endAtomId())->index() + 1, bond->order());
//      qDebug() << "Bonds" << obmol.NumBonds();
    }
    obmol.EndModify();

    qDebug() << "OBMol() run" << obmol.NumAtoms() << obmol.NumBonds();

    return obmol;
  }

  bool Molecule::setOBMol(OpenBabel::OBMol *obmol)
  {
    // Take an OBMol, copy everything we need and store this object
    Q_D(Molecule);
    d->obmol = obmol;
    // Copy all the parts of the OBMol to our Molecule
    std::vector<OpenBabel::OBNodeBase*>::iterator i;
    for (OpenBabel::OBAtom *obatom = static_cast<OpenBabel::OBAtom *>(obmol->BeginAtom(i));
          obatom; obatom = static_cast<OpenBabel::OBAtom *>(obmol->NextAtom(i))) {
      Atom *atom = newAtom();
      atom->setOBAtom(obatom);
      qDebug() << "Old atom:" << obatom->GetIdx() << obatom->GetX() << obatom->GetY() << obatom->GetZ() << obatom->GetAtomicNum();

      qDebug() << "New atom:" << atom->index() << atom->pos().x() << atom->pos().y() << atom->pos().z() << atom->atomicNumber();
    }
    // Now bonds, we use the indices of the atoms to get the bonding right
    std::vector<OpenBabel::OBEdgeBase*>::iterator j;
    for (OpenBabel::OBBond *obbond = static_cast<OpenBabel::OBBond*>(obmol->BeginBond(j));
         obbond; obbond = static_cast<OpenBabel::OBBond*>(obmol->NextBond(j))) {
      Bond *bond = newBond();
      bond->setOBBond(obbond);
      // Get the begin and end atoms - we use a 0 based index, OB uses 1 based
      bond->setBegin(atom(obbond->GetBeginAtom()->GetIdx()-1));
      bond->setEnd(atom(obbond->GetEndAtom()->GetIdx()-1));
      // Set the bond to the atoms too, remember the 0 based and 1 based arrays
      atom(obbond->GetBeginAtom()->GetIdx()-1)->addBond(bond);
      atom(obbond->GetEndAtom()->GetIdx()-1)->addBond(bond);
    }
    return true;
  }

  const Eigen::Vector3d & Molecule::center() const
  {
    Q_D(const Molecule);
    if( d->invalidGeomInfo ) computeGeomInfo();
    return d->center;
  }

  const Eigen::Vector3d & Molecule::normalVector() const
  {
    Q_D(const Molecule);
    if( d->invalidGeomInfo ) computeGeomInfo();
    return d->normalVector;
  }

  const double & Molecule::radius() const
  {
    Q_D(const Molecule);
    if( d->invalidGeomInfo ) computeGeomInfo();
    return d->radius;
  }

  const Atom * Molecule::farthestAtom() const
  {
    Q_D(const Molecule);
    if( d->invalidGeomInfo ) computeGeomInfo();
    return d->farthestAtom;
  }

  void Molecule::clear()
  {
    Q_D(Molecule);
    d->atoms.resize(0);
    d->bonds.resize(0);
    foreach (Atom *atom, d->atomList) {
      atom->deleteLater();
      emit primitiveRemoved(atom);
    }
    d->atomList.clear();
    foreach (Bond *bond, d->bondList) {
      bond->deleteLater();
      emit primitiveRemoved(bond);
    }
    d->bondList.clear();
  }

  Molecule &Molecule::operator=(const Molecule& other)
  {
    Q_D(Molecule);
    clear();
    d->autoId = false;
    const MoleculePrivate *e = other.d_func();
    d->atoms.resize(e->atoms.size(), 0);
    d->bonds.resize(e->bonds.size(), 0);

    // Copy the atoms and bonds over
    for (unsigned int i = 0; i < e->atoms.size(); ++i) {
      if (e->atoms.at(i) > 0) {
        Atom *atom = new Atom;
        *atom = *(e->atoms[i]);
        atom->setId(e->atoms[i]->id());
        atom->setIndex(e->atoms[i]->index());
        d->atoms[i] = atom;
        d->atomList.push_back(atom);
        emit primitiveAdded(atom);
      }
    }

    for (unsigned int i = 0; i < e->bonds.size(); ++i) {
      if (e->bonds.at(i)) {
        Bond *bond = new Bond;
        *bond = *(e->bonds[i]);
        bond->setId(e->bonds[i]->id());
        bond->setIndex(e->bonds[i]->index());
        d->bonds[i] = bond;
        d->bondList.push_back(bond);
        emit primitiveAdded(bond);
      }
    }

    d->autoId = true;

    return *this;
  }

  Molecule &Molecule::operator+=(const Molecule& other)
  {
    const MoleculePrivate *e = other.d_func();
    // Create a temporary map from the old indices to the new for bonding
    QList<int> map;
    foreach (Atom *a, e->atomList) {
      Atom *atom = newAtom();
      *atom = *a;
      map.push_back(atom->id());
      emit primitiveAdded(atom);
    }
    foreach (Bond *b, e->bondList) {
      Bond *bond = newBond();
      *bond = *b;
      bond->setBegin(getAtomById(map.at(other.getAtomById(b->beginAtomId())->index())));
      bond->setEnd(getAtomById(map.at(other.getAtomById(b->endAtomId())->index())));
      emit primitiveAdded(bond);
    }
    return *this;
  }

  void Molecule::computeGeomInfo() const
  {
    Q_D(const Molecule);
    d->invalidGeomInfo = true;
    d->farthestAtom = 0;
    d->center.setZero();
    d->normalVector.setZero();
    d->radius = 0.0;
    if(numAtoms() != 0)
    {
      // compute center
      foreach (Atom *atom, d->atomList)
        d->center += atom->pos();

      d->center /= numAtoms();

      // compute the normal vector to the molecule's best-fitting plane
      int i = 0;
      Eigen::Vector3d ** atomPositions = new Eigen::Vector3d*[numAtoms()];
      foreach (Atom *atom, d->atomList)
        atomPositions[i++] = const_cast<Eigen::Vector3d*>(&atom->pos());

      Eigen::Hyperplane<double, 3> planeCoeffs;
      Eigen::fitHyperplane(numAtoms(), atomPositions, &planeCoeffs);
      delete[] atomPositions;
      d->normalVector = planeCoeffs.normal();

      // compute radius and the farthest atom
      d->radius = -1.0; // so that ( squaredDistanceToCenter > d->radius ) is true for at least one atom.
      foreach (Atom *atom, d->atomList) {
        double distanceToCenter = (atom->pos() - d->center).norm();
        if(distanceToCenter > d->radius) {
          d->radius = distanceToCenter;
          d->farthestAtom = atom;
        }
      }
    }
    d->invalidGeomInfo = false;
  }

}

#include "primitive.moc"
