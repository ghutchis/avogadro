/**********************************************************************
  HBondEngine - Hydrogen Bond Engine

  Copyright (C) 2007 by Marcus D. Hanwell
  Copyright (C) 2007 by Tim Vandermeersch

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

#include "hbondengine.h"
#include <config.h>

#include <avogadro/primitive.h>
#include <avogadro/atom.h>
#include <avogadro/molecule.h>
#include <avogadro/color.h>
#include <avogadro/glwidget.h>

#include <openbabel/mol.h>
#include <openbabel/obiter.h>

#include <QMessageBox>
#include <QString>
#include <QDebug>

using namespace std;
using namespace OpenBabel;
using namespace Eigen;

namespace Avogadro {

  HBondEngine::HBondEngine(QObject *parent) : Engine(parent), m_settingsWidget(0),
                                              m_width(2), m_radius(2.0), m_angle(120)
  {
    setDescription(tr("Renders hydrogen bonds"));
  }

  Engine *HBondEngine::clone() const
  {
    HBondEngine *engine = new HBondEngine(parent());
    engine->setAlias(alias());
    engine->setWidth(m_width);
    engine->setRadius(m_radius);
    engine->setAngle(m_angle);
    engine->setEnabled(isEnabled());

    return engine;
  }

  HBondEngine::~HBondEngine()
  {
  }

  bool HBondEngine::renderOpaque(PainterDevice *pd)
  {
    Molecule *molecule = const_cast<Molecule *>(pd->molecule());
    OBMol mol = molecule->OBMol();

    pd->painter()->setColor(1.0, 1.0, 1.0);
    int stipple = 0xF0F0; // pattern for lines

    FOR_PAIRS_OF_MOL (p, mol) {
      OBAtom *a = mol.GetAtom((*p)[0]);
      OBAtom *b = mol.GetAtom((*p)[1]);

      if (a->GetDistance(b) > m_radius)
        continue;

      if (a->GetDistance(b) < 0.7) // too close
        continue;

      if (a->IsHbondDonorH() && b->IsHbondAcceptor()) {
        double angle = 180.0; // default, if no neighbours on H
        FOR_NBORS_OF_ATOM (c, a)
          angle = c->GetAngle(a, b);
        if (angle < m_angle)
          continue;

        OBAtom *atom1 = mol.GetAtom((*p)[0]);
        OBAtom *atom2 = mol.GetAtom((*p)[1]);
        Vector3d v1(atom1->GetX(), atom1->GetY(), atom1->GetZ());
        Vector3d v2(atom2->GetX(), atom2->GetY(), atom2->GetZ());
        pd->painter()->drawMultiLine(v1, v2, m_width, 1, stipple);
      } else if (b->IsHbondDonorH() && a->IsHbondAcceptor()) {
        double angle = 180.0; // default, if no neighbours on H
      	FOR_NBORS_OF_ATOM (c, b)
          angle = c->GetAngle(b, a);
        if (angle < m_angle)
          continue;

        OBAtom *atom1 = mol.GetAtom((*p)[0]);
        OBAtom *atom2 = mol.GetAtom((*p)[1]);
        Vector3d v1(atom1->GetX(), atom1->GetY(), atom1->GetZ());
        Vector3d v2(atom2->GetX(), atom2->GetY(), atom2->GetZ());
        pd->painter()->drawMultiLine(v1, v2, m_width, 1, stipple);
      }
    }

    return true;
  }

  double HBondEngine::radius(const PainterDevice *, const Primitive *) const
  {
    return 0.0;
  }

  QWidget* HBondEngine::settingsWidget()
  {
    if(!m_settingsWidget)
      {
        m_settingsWidget = new HBondSettingsWidget();
        connect(m_settingsWidget->widthSlider, SIGNAL(valueChanged(int)), this, SLOT(setWidth(int)));
        connect(m_settingsWidget->radiusSpin, SIGNAL(valueChanged(double)), this, SLOT(setRadius(double)));
        connect(m_settingsWidget->angleSpin, SIGNAL(valueChanged(double)), this, SLOT(setAngle(double)));
        connect(m_settingsWidget, SIGNAL(destroyed()), this, SLOT(settingsWidgetDestroyed()));
        m_settingsWidget->widthSlider->setValue(m_width);
        m_settingsWidget->radiusSpin->setValue(m_radius);
        m_settingsWidget->angleSpin->setValue(m_angle);
      }
    return m_settingsWidget;
  }

  void HBondEngine::setWidth(int value)
  {
    m_width = (double) value;
    emit changed();
  }

  void HBondEngine::setRadius(double value)
  {
    m_radius = value;
    emit changed();
  }

  void HBondEngine::setAngle(double value)
  {
    m_angle = value;
    emit changed();
  }


  void HBondEngine::settingsWidgetDestroyed()
  {
    qDebug() << "Destroyed Settings Widget";
    m_settingsWidget = 0;
  }

  void HBondEngine::writeSettings(QSettings &settings) const
  {
    Engine::writeSettings(settings);
    settings.setValue("width", m_width);
    settings.setValue("radius", m_radius);
    settings.setValue("angle", m_angle);
  }

  void HBondEngine::readSettings(QSettings &settings)
  {
    Engine::readSettings(settings);
    setWidth(settings.value("width", 2.0).toDouble());
    setRadius(settings.value("radius", 2.0).toDouble());
    setAngle(settings.value("angle", 120.0).toDouble());
    if (m_settingsWidget) {
      m_settingsWidget->widthSlider->setValue(m_width);
      m_settingsWidget->radiusSpin->setValue(m_radius);
      m_settingsWidget->angleSpin->setValue(m_angle);
    }
  }

}

#include "hbondengine.moc"

Q_EXPORT_PLUGIN2(hbondengine, Avogadro::HBondEngineFactory)
