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


#ifndef MOVIEDIALOG_H
#define MOVIEDIALOG_H

#include <QDialog>
#include <QButtonGroup>
#include <QModelIndex>

#include "ui_moviedialog.h"

namespace Avogadro
{
  class MovieDialog : public QDialog
  {
      Q_OBJECT

    public:
      //! Constructor
      explicit MovieDialog( QWidget *parent = 0, Qt::WindowFlags f = 0 );
      //! Deconstructor
      ~MovieDialog();

    private:
      Ui::MovieDialog ui;

    public Q_SLOTS:
      void saveMovie();

    Q_SIGNALS:
      // signal to emit info needed by animation extension to get movie params
      // for now, only has filename
      void movieFileInfo(QString filename);
  };
}

#endif
