/**********************************************************************
  MovieDialog - dialog for making a movie from snapshots

  Copyright (C) 2008 by Naomi Fox

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

#include "moviedialog.h"

#include <QPushButton>
#include <QButtonGroup>
#include <QDebug>

#include <QFileDialog>
#include <QFile>

#include <QMessageBox>
#include <QInputDialog>

#include <openbabel/plugin.h>

using namespace OpenBabel;

namespace Avogadro {

  MovieDialog::MovieDialog( QWidget *parent, Qt::WindowFlags f ) : QDialog( parent, f )
  {
    ui.setupUi(this);
    connect(ui.movieSaveButton, SIGNAL(clicked()), this, SLOT(saveMovie()));
  }

  MovieDialog::~MovieDialog()
  {
    //  qDebug() << "MovieDialog::~MovieDialog()" << endl;
  }

  void MovieDialog::saveMovie()
  {
    QString sMovieFileName = QFileDialog::getSaveFileName(this, tr("Save Movie File"),
                                "/home",
                                tr("Movie (*.avi)"));

    emit movieFileInfo(sMovieFileName);
  }
  
}

#include "moviedialog.moc"
