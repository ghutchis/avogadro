/**********************************************************************
  Animation - Basic animation

  Copyright (C) 2008 by Tim Vandermeersch

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

#include "animationextension.h"
#include "trajvideomaker.h"
#include <avogadro/primitive.h>
#include <avogadro/color.h>
#include <avogadro/glwidget.h>

#include <openbabel/obconversion.h>

#include <QMessageBox>

using namespace OpenBabel;

namespace Avogadro {

  AnimationExtension::AnimationExtension(QObject *parent) : Extension(parent), m_molecule(0),
							    m_animationDialog(0), m_timeLine(0), m_frameCount(0), m_widget(0)
  {
    QAction *action = new QAction(this);
    action->setText(tr("Animation..."));
    m_actions.append(action);
    
    action = new QAction( this );
    action->setSeparator(true);
    m_actions.append(action);
  }

  AnimationExtension::~AnimationExtension()
  {
    if (m_animationDialog)
    {
      delete m_animationDialog;
      m_animationDialog = 0;
    }
 
    if (m_timeLine)
    {
      delete m_timeLine;
      m_timeLine = 0;
    }
  }

  QList<QAction *> AnimationExtension::actions() const
  {
    return m_actions;
  }

  QString AnimationExtension::menuPath(QAction *) const
  {
    return tr("&Extensions");
  }

  void AnimationExtension::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QUndoCommand* AnimationExtension::performAction(QAction *, GLWidget* widget)
  {
    m_widget = widget;

    if (!m_animationDialog)
    {
      m_timeLine = new QTimeLine;
      m_animationDialog = new AnimationDialog;

      connect(m_animationDialog, SIGNAL(fileName(QString)), this, SLOT(loadFile(QString)));
      connect(m_animationDialog, SIGNAL(sliderChanged(int)), this, SLOT(setFrame(int)));
      connect(m_animationDialog, SIGNAL(fpsChanged(int)), this, SLOT(setDuration(int)));
      connect(m_animationDialog, SIGNAL(loopChanged(int)), this, SLOT(setLoop(int)));

      connect(m_timeLine, SIGNAL(frameChanged(int)), this, SLOT(setFrame(int)));
      connect(m_animationDialog, SIGNAL(play()), m_timeLine, SLOT(start()));
      connect(m_animationDialog, SIGNAL(pause()), m_timeLine, SLOT(stop()));
      connect(m_animationDialog, SIGNAL(stop()), this, SLOT(stop()));
      connect(m_animationDialog, SIGNAL(videoFileInfo(QString)), this, SLOT(saveVideo(QString)));
    } 

    m_animationDialog->show();

    return 0;
  }

  void AnimationExtension::loadFile(QString file)
  {
    //qDebug() << "AnimationExtension::loadFile()" << endl;

    if (file.isEmpty())
      return;

    OBConversion conv;
    OBFormat *inFormat = conv.FormatFromExt(( file.toAscii() ).data() );
    if ( !inFormat || !conv.SetInFormat( inFormat ) ) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Cannot read file format of file %1." )
          .arg( file ) );
      return;
    }

    if (!conv.ReadFile(m_molecule, file.toStdString())) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Read trajectory file %1 failed." )
          .arg( file ) );
      return;
    }

    m_frameCount = m_molecule->NumConformers();
    m_animationDialog->setFrameCount(m_frameCount);
    m_animationDialog->setFrame(1);
    m_timeLine->setFrameRange(1, m_frameCount);
    setDuration(m_animationDialog->fps());
  }

  void AnimationExtension::setDuration(int i)
  {
    int interval = 1000 / i;
    m_timeLine->setUpdateInterval(interval);
    int duration = interval * m_frameCount;
    m_timeLine->setDuration(duration);
  }

  void AnimationExtension::setLoop(int state)
  {
    if (state == Qt::Checked) {
      m_timeLine->setLoopCount(0);
    } else {
      m_timeLine->setLoopCount(1);
    }
  }

  void AnimationExtension::setFrame(int i)
  {
    // if (m_timeLine->state() != QTimeLine::Running)
    //  m_timeLine->setCurrentTime(m_timeLine->updateInterval() * i);

    m_animationDialog->setFrame(i);
    m_molecule->SetConformer(i - 1);
    m_molecule->update();
  }

  void AnimationExtension::stop()
  {
    m_timeLine->stop();
    m_timeLine->setCurrentTime(0);
    setFrame(1);
  }

  void AnimationExtension::saveVideo(QString videoFileName) 
  {       
    if (videoFileName.isEmpty()) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "Must specify a valid .avi file name" ));
      return;
    }
    
    if (!videoFileName.endsWith(".avi")){
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "Must specify a valid .avi file name" ));
      return;
    }
    
    if (!m_widget) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "GL widget was not correctly initialized in order to save video" ));
      return;
    }
    
    //first, split out the directory and filenames
    QString dir, fileName, prefix;

    int slashPos = videoFileName.lastIndexOf("/");
    
    if (slashPos < 0) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "Invalid video filename.  Must include full directory path" ));
      return;
    }

    dir = videoFileName.left(slashPos) + "/"; 
    fileName = videoFileName.right(videoFileName.length() - (slashPos+1));
    if (fileName.isEmpty()) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "Invalid video filename.  Must include full directory path and name, ending with .avi" ));
      return;
    }
    
    //if (fileName.endsWith(".avi")) {
    prefix = fileName.left(fileName.length() - 4);
    

    //Make the directory where the snapshots will be saved
    QString snapshotsDir = dir + prefix + "/";
    QString mkdirCommand = "mkdir " + snapshotsDir;
    system(mkdirCommand.toStdString().c_str());

    TrajVideoMaker::makeVideo(m_widget, snapshotsDir, videoFileName);

  }


} // end namespace Avogadro

#include "animationextension.moc"
Q_EXPORT_PLUGIN2(animationextension, Avogadro::AnimationExtensionFactory)
