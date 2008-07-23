/**********************************************************************
  ManipulateTool - Manipulation Tool for Avogadro

  Copyright (C) 2007 by Marcus D. Hanwell
  Copyright (C) 2007 by Geoffrey R. Hutchison
  Copyright (C) 2007 by Benoit Jacob

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

#ifndef APPLYFORCETOOL_H
#define APPLYFORCETOOL_H

//#include <eigen/vector.h>

#include <avogadro/glwidget.h>
#include <avogadro/tool.h>

#include <openbabel/mol.h>


#include <QGLWidget>
#include <QObject>
#include <QStringList>
#include <QImage>
#include <QAction>
#include <QUndoCommand>



namespace Avogadro {

  /**
   * @class ManipulateTool
   * @brief Manipulate the position of atoms
   * @author Marcus D. Hanwell
   *
   * This tool enables the manipulation of the position of
   * the selected atoms.
   */
  class ApplyForceTool : public Tool
  {
    Q_OBJECT

    public:
      //! Constructor
      ApplyForceTool(QObject *parent = 0);
      //! Deconstructor
      virtual ~ApplyForceTool();

      //! \name Description methods
      //@{
      //! Tool Name (ie Draw)
      virtual QString name() const { return(tr("Apply Force")); }
      //! Tool Description (ie. Draws atoms and bonds)
      virtual QString description() const { return(tr("Apply Force Tool")); }
      //@}

      //! \name Tool Methods
      //@{
      //! \brief Callback methods for ui.actions on the canvas.
      /*!
      */
      virtual QUndoCommand* mousePress(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* mouseRelease(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* mouseMove(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* wheel(GLWidget *widget, const QWheelEvent *event);

      virtual int usefulness() const;

      virtual bool paint(GLWidget *widget);
      
 
  protected:
      Atom *              m_clickedAtom;
      bool                m_leftButtonPressed;  // rotation
      bool                m_midButtonPressed;   // scale / zoom
      bool                m_rightButtonPressed; // translation
      
      QPoint              m_lastDraggingPosition;
      
      
      //nkf
      //the current force vector
      //I can't get this vector3d working
      Eigen::Vector3d       m_currentForce;

      double m_forceX;
      double m_forceY;
      double m_forceZ;

      void updateForce(GLWidget *widget, const Eigen::Vector3d &what, const QPoint &from, const QPoint &to);
  };


  class ApplyForceToolFactory : public QObject, public ToolFactory
    {
      Q_OBJECT
      Q_INTERFACES(Avogadro::ToolFactory)

      public:
        Tool *createInstance(QObject *parent = 0) { return new ApplyForceTool(parent); }
    };

} // end namespace Avogadro

#endif
