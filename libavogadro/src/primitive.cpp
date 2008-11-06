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

namespace Avogadro {

  class PrimitivePrivate {
    public:
      PrimitivePrivate() {}
  };

  Primitive::Primitive(QObject *parent) : QObject(parent),
    d_ptr(new PrimitivePrivate), m_type(Primitive::OtherType), m_id(-1),
    m_lock(new QReadWriteLock) {}

  Primitive::Primitive(enum Type type, QObject *parent) : QObject(parent),
    d_ptr(new PrimitivePrivate), m_type(type), m_id(-1), m_lock(new QReadWriteLock)
  {}

  Primitive::Primitive(PrimitivePrivate &dd, QObject *parent) : QObject(parent), d_ptr(&dd) {}

  Primitive::Primitive(PrimitivePrivate &dd, enum Type type, QObject *parent) : QObject(parent), d_ptr(&dd), m_type(type), m_id(-1), m_lock(new QReadWriteLock) {}

  Primitive::~Primitive()
  {
    delete d_ptr;
    delete m_lock;
  }

  enum Primitive::Type Primitive::type() const
  {
    return m_type;
  }

  QReadWriteLock *Primitive::lock()
  {
    return m_lock;
  }

  void Primitive::update()
  {
    emit updated();
  }

  void Primitive::setId(unsigned long id)
  {
    m_id = id;
  }

  unsigned long Primitive::id() const
  {
    return m_id;
  }

  void Primitive::setIndex(unsigned long index)
  {
    m_index = index;
  }

  unsigned long Primitive::index() const
  {
    return m_index;
  }

}

#include "primitive.moc"
