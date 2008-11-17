/**********************************************************************
  POVRayExtension - Extension for generating POV-Ray rendered images

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

#ifndef POVRAYEXTENSION_H
#define POVRAYEXTENSION_H

#include "povraydialog.h"

#include <avogadro/glwidget.h>
#include <avogadro/extension.h>

class QProcess;

namespace Avogadro
{
  class POVRayExtension : public Extension
  {
  Q_OBJECT

  public:
    POVRayExtension(QObject* parent = 0);
    virtual ~POVRayExtension();

    virtual QString name() const { return QObject::tr("POV-Ray"); }
    virtual QString description() const
    {
      return QObject::tr("Export images rendered using POV-Ray");
    }

    virtual QList<QAction *> actions() const;

    virtual QString menuPath(QAction* action) const;

    virtual QUndoCommand* performAction(QAction *action, GLWidget *widget);

    void setMolecule(Molecule *molecule);

  private:
    GLWidget* m_glwidget;
    POVRayDialog* m_POVRayDialog;
    QList<QAction *> m_actions;
    Molecule *m_molecule;
    QProcess *m_process;

  private Q_SLOTS:
    void render();

  };

  class POVRayExtensionFactory : public QObject, public PluginFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::PluginFactory)
    AVOGADRO_EXTENSION_FACTORY(POVRayExtension,
        tr("POV-Ray Extension"),
        tr("Extension for creating POV-Ray files and"
          " rendering them using the command line POV-Ray program."))

  };

} // End namespace Avogadro

#endif
