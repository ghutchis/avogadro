/**********************************************************************
  AnimateMolDialog - Dialog for animation extension

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

#ifndef ANIMATEMOLDIALOG_H
#define ANIMATEMOLDIALOG_H

#include <QDialog>
#include <QButtonGroup>
#include <QModelIndex>

#include "ui_animatemoldialog.h"

namespace Avogadro
{
  class AnimateMolDialog : public QDialog
  {
      Q_OBJECT

    public:
      //! Constructor
      explicit AnimateMolDialog( QWidget *parent = 0, Qt::WindowFlags f = 0 );
      //! Desconstructor
      ~AnimateMolDialog();

      void setFrame(int i);
      void setFrameCount(int i);

      int fps();

    private:
      Ui::AnimateMolDialog ui;
      int m_frameCount;

    public Q_SLOTS:
      void loadFile();

    Q_SIGNALS:
      void fileName(QString filename);
      void sliderChanged(int i);
      void fpsChanged(int i);
      bool loopChanged(int state);
      void back();
      void play();
      void pause();
      void stop();
      void forward();
  };
}

#endif
