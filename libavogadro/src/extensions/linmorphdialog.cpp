/**********************************************************************
  AnimationDialog - Dialog for animation extension

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

#include "linmorphdialog.h"

#include <QPushButton>
#include <QButtonGroup>
#include <QDebug>

#include <QFileDialog>
#include <QFile>

#include <QMessageBox>

#include <openbabel/plugin.h>

using namespace OpenBabel;

namespace Avogadro {

  LinMorphDialog::LinMorphDialog( QWidget *parent, Qt::WindowFlags f ) : QDialog( parent, f )
  {
    //  qDebug() << "LinMorphDialog::LinMorphDialog()" << endl;
    ui.setupUi(this);
    connect(ui.loadButton, SIGNAL(clicked()), this, SLOT(loadFile()));
    connect(ui.frameSlider, SIGNAL(valueChanged(int)), this, SIGNAL(sliderChanged(int)));
    connect(ui.fpsSpin, SIGNAL(valueChanged(int)), this, SIGNAL(fpsChanged(int)));
    connect(ui.loopBox, SIGNAL(stateChanged(int)), this, SIGNAL(loopChanged(int)));
    connect(ui.numFramesSpin, SIGNAL(valueChanged(int)), this, SIGNAL(frameCountChanged(int)));    
    connect(ui.saveTrajectoryButton, SIGNAL(clicked()), this, SLOT(saveTrajFile()));
    connect(ui.playButton, SIGNAL(clicked()), this, SIGNAL(play()));
    connect(ui.pauseButton, SIGNAL(clicked()), this, SIGNAL(pause()));
    connect(ui.stopButton, SIGNAL(clicked()), this, SIGNAL(stop()));


  }

  LinMorphDialog::~LinMorphDialog()
  {
    //  qDebug() << "LinMorphDialog::~LinMorphDialog()" << endl;
  }

  void LinMorphDialog::loadFile()
  {
    // Load a file
    QString file = QFileDialog::getOpenFileName(this,
      tr("Open mol file"), ui.fileEdit->text(),
						tr("Any files (*.*)"));
						//tr("Trajectory files (*.xtc)"));
    ui.fileEdit->setText(file);
    
    emit fileName(file);
  }


  void LinMorphDialog::saveTrajFile()
  {
    // Load a file
    QString file = QFileDialog::getSaveFileName(this,
      tr("Traj file to save to"), ui.fileEdit->text(),
						tr("Trajectory files (*.xtc)"));
    //ui.fileEdit->setText(file);
    
    emit trajFileName(file);
  }

  void LinMorphDialog::setFrame(int i)
  {
    QString str = tr("%1/%2").arg( i ).arg( m_frameCount );

    ui.frameEdit->setText(str);
    ui.frameSlider->setValue(i);
  }

  void LinMorphDialog::setFrameCount(int i)
  {
    m_frameCount = i;
    ui.frameSlider->setMaximum(i);
    setFrame(1);
  }

  int LinMorphDialog::fps()
  {
    return ui.fpsSpin->value();
  }

}

#include "linmorphdialog.moc"
