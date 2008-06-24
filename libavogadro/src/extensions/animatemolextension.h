/**********************************************************************
 Template - Extension Template

  Copyright (C) 2008 by Author

  This file is part of the Avogadro molecular editor project.
  For more information, see <http://avogadro.sourceforge.net/>

  Some code is based on Open Babel
  For more information, see <http://openbabel.sourceforge.net/>

  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation version 2 of the License.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
 ***********************************************************************/

#ifndef ANIMATEMOLEXTENSION_H
#define ANIMATEMOLEXTENSION_H

#include <avogadro/extension.h>
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

#include <QTimeLine>
#include "animatemoldialog.h"

namespace Avogadro {

  class AnimateMolExtension : public Extension
  {
    Q_OBJECT

    public:
      //! constructor
      AnimateMolExtension(QObject *parent=0);
      //! Deconstructor
      virtual ~AnimateMolExtension();

      virtual QString name() const { return QObject::tr("Animate somethinig Mol"); }
      virtual QString description() const { return QObject::tr("Animate Mol Extension.  This is an abstract class."); };

      virtual QList<QAction *> actions() const;
      virtual QString menuPath(QAction *action) const;

      virtual QDockWidget * dockWidget();
      virtual QUndoCommand* performAction(QAction *action, GLWidget *widget);

      virtual void setMolecule(Molecule *molecule);
      


  protected:
      Molecule *m_molecule;
      

      QList<QAction *> m_actions;
      AnimateMolDialog *m_animateMolDialog;
      QTimeLine *m_timeLine;

      int m_frameCount;
    
    protected Q_SLOTS:
      virtual void loadFile(QString file)=0;
      void writeTrajectoryFile(QString file); 
      void setDuration(int i);
      void setLoop(int state);
      void setFrame(int i);
      void stop();

  };

} // end namespace Avogadro

#endif
