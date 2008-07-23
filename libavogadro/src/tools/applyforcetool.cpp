/**********************************************************************
  
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

#include "applyforcetool.h"

#include <avogadro/navigate.h>
#include <avogadro/primitive.h>
#include <avogadro/color.h>
#include <avogadro/glwidget.h>
#include <avogadro/camera.h>

#include <avogadro/cylinder.h>

#include <openbabel/obiter.h>
#include <openbabel/mol.h>

#include <QtPlugin>
#include <QDebug>
#include <QString>

using namespace std;
using namespace OpenBabel;
using namespace Eigen;



#define TESS_LEVEL 32
#define RIBBON_WIDTH 0.05
#define RIBBON_LENGTH 0.6
#define RIBBON_ARROW_WIDTH 0.15
#define RIBBON_ARROW_LENGTH 0.25
#define RIBBON_APERTURE 0.07
#define MINIMUM_APPARENT_SIZE 0.04
#define MAXIMUM_APPARENT_SIZE 0.25
#define SIZE_FACTOR_WHEN_NOTHING_CLICKED 0.25
#define ZOOM_SIZE_FACTOR 0.3
#define ATOM_SIZE_FACTOR 1.1

namespace Avogadro {

  ApplyForceTool::ApplyForceTool(QObject *parent) : Tool(parent), m_clickedAtom(0), m_leftButtonPressed(false), m_midButtonPressed(false), m_rightButtonPressed(false), m_forceX(0.), m_forceY(0.), m_forceZ(0.) 
  {
  
    QAction *action = activateAction();
    action->setIcon(QIcon(QString::fromUtf8(":/manipulate/manipulate.png")));
    action->setToolTip(tr(" Apply Force Tool \n\n"
          "Left Mouse:   Click and drag to move atoms\n"
			  "Middle Mouse: Click and drag to move atoms further away or closer\n")); 
    //action->setShortcut(Qt::Key_F10);
  }

  ApplyForceTool::~ApplyForceTool()
  {

  }

  /// apply a force... i.e. draw a vector showing magnitude and direction 
  /// for now.
  void ApplyForceTool::updateForce(GLWidget *widget, const Eigen::Vector3d &what,
                                 const QPoint &from, const QPoint &to) 
  {
    // Set the cursor - this needs to be reset to Qt::ArrowCursor after
    // Currently, there's a Qt/Mac bug -- SizeAllCursor looks like a spreadsheet cursor
#ifdef Q_WS_MAC
    widget->setCursor(Qt::CrossCursor);
#else 
    widget->setCursor(Qt::SizeAllCursor);
#endif

    // Translate the selected atoms in the x and y sense of the view
    Eigen::Vector3d fromPos = widget->camera()->unProject(from, what);
    Eigen::Vector3d toPos = widget->camera()->unProject(to, what);

    MatrixP3d atomTranslation;
    atomTranslation.loadTranslation(toPos - fromPos);
 
    m_currentForce = toPos - what;

    QString qS = "Force to be applied at atom :" ;
    qS += " " + QString::number(m_currentForce.x()); 
    qS += " " + QString::number(m_currentForce.y()); 
    qS += " " + QString::number(m_currentForce.z()); 
    qDebug() << qS;


    /*
    widget->molecule()->BeginModify();
    if (widget->selectedPrimitives().size())
      foreach(Primitive *p, widget->selectedPrimitives())
        if (p->type() == Primitive::AtomType)
          static_cast<Atom *>(p)->setPos(atomTranslation * static_cast<Atom *>(p)->pos());
    if (m_clickedAtom && !widget->isSelected(m_clickedAtom))
      m_clickedAtom->setPos(atomTranslation * m_clickedAtom->pos());
    widget->molecule()->EndModify();
    widget->molecule()->update();
    */
  }

  QUndoCommand* ApplyForceTool::mousePress(GLWidget *widget, const QMouseEvent *event)
  {
    m_lastDraggingPosition = event->pos();
    m_currentForce.loadZero();

    // Make sure there aren't modifier keys clicked with the left button
    // If the user has a Mac and only a one-button mouse, everything
    // looks like a left button
    if (event->buttons() & Qt::LeftButton &&
        event->modifiers() == Qt::NoModifier)
    {
      m_leftButtonPressed = true;
      // Set the cursor - this needs to be reset to Qt::ArrowCursor after
      // Currently, there's a Qt/Mac bug -- SizeAllCursor looks like a spreadsheet cursor
#ifdef Q_WS_MAC
      widget->setCursor(Qt::CrossCursor);
#else 
      widget->setCursor(Qt::SizeAllCursor);
#endif
    }

    // On a Mac, click and hold the Shift key
    if (event->buttons() & Qt::MidButton ||
        (event->buttons() & Qt::LeftButton &&
         event->modifiers() == Qt::ShiftModifier))
    {
      m_midButtonPressed = true;
      // Set the cursor - this needs to be reset to Qt::ArrowCursor after
      widget->setCursor(Qt::SizeVerCursor);
    }

    // On a Mac, click and hold either the Command or Control Keys
    // (Control or Meta in Qt-speak)
    if (event->buttons() & Qt::RightButton ||
        (event->buttons() & Qt::LeftButton &&
        (event->modifiers() == Qt::ControlModifier
         || event->modifiers() == Qt::MetaModifier)))
    {
      m_rightButtonPressed = true;
      // Set the cursor - this needs to be reset to Qt::ArrowCursor after
      widget->setCursor(Qt::ClosedHandCursor);
    }

    m_clickedAtom = widget->computeClickedAtom(event->pos());

    widget->update();

    //QUndoCommand* undo = new MoveAtomCommand(widget->molecule());
    //return undo;
    return NULL;
  }

  QUndoCommand* ApplyForceTool::mouseRelease(GLWidget *widget, const QMouseEvent*)
  {
    m_leftButtonPressed = false;
    m_midButtonPressed = false;
    m_rightButtonPressed = false;
    m_clickedAtom = 0;

    // Set the cursor back to the default cursor
    widget->setCursor(Qt::ArrowCursor);
    widget->update();
    //QUndoCommand* undo = new MoveAtomCommand(widget->molecule());
    //return undo;
    return NULL;
  }

  QUndoCommand* ApplyForceTool::mouseMove(GLWidget *widget, const QMouseEvent *event)
  {
    if(!widget->molecule())
      return 0;

    QPoint deltaDragging = event->pos() - m_lastDraggingPosition;

    // Manipulation can be performed in two ways - centred on an individual atom

    if (m_clickedAtom)
    {
      if (m_leftButtonPressed)
      {
        // translate the molecule following mouse movement
        updateForce(widget, m_clickedAtom->pos(), m_lastDraggingPosition, event->pos());
      }
    }

    m_lastDraggingPosition = event->pos();
    widget->update();

    

    return 0;
  }

  QUndoCommand* ApplyForceTool::wheel(GLWidget*widget, const QWheelEvent*event)
  {
    return NULL;
  }

  bool ApplyForceTool::paint(GLWidget *widget)
  {
    if (m_clickedAtom) {
    //draw the force vector
      
    Vector3d center = m_clickedAtom->pos();
    double size = qMax(widget->radius(m_clickedAtom) * ATOM_SIZE_FACTOR,
             MINIMUM_APPARENT_SIZE * widget->camera()->distance(center));
    double shift = widget->radius(m_clickedAtom);

    // Set up the axes and some vectors to work with
    Vector3d xAxis = widget->camera()->backTransformedXAxis();
    Vector3d yAxis = widget->camera()->backTransformedYAxis();
    Vector3d zAxis = widget->camera()->backTransformedZAxis();
    Vector3d v;

    // Horizontal arrow, pointing left
    glDisable(GL_LIGHTING);
    Color arrowColor(1.0, 1.0, 0.3, 0.7);
    arrowColor.apply();
    
    /// nkf - need to make all this stuff prettier using lin alg
    v = center + shift * zAxis;
    glBegin(GL_QUAD_STRIP);
    glVertex3dv((v + RIBBON_WIDTH*size*yAxis).array());
    glVertex3dv((v - RIBBON_WIDTH*size*yAxis).array());
    v += m_currentForce;
    glVertex3dv((v - RIBBON_WIDTH*size*yAxis).array());
    glVertex3dv((v + RIBBON_WIDTH*size*yAxis).array());
    glEnd();
    glBegin(GL_TRIANGLES);
    glVertex3dv((v + RIBBON_ARROW_WIDTH*size*yAxis).array());
    glVertex3dv((v - RIBBON_ARROW_WIDTH*size*yAxis).array());
    glVertex3dv((v + RIBBON_ARROW_LENGTH*size*xAxis).array());
    glEnd();
    glEnable(GL_LIGHTING);
    }
    // this will draw the force vector with the center of  atom selected 
    // as the origin
   
    /*
    int selectedSize = widget->selectedPrimitives().size();
    if(m_clickedAtom)
    {
      if(m_leftButtonPressed)
      {
        m_eyecandy->drawTranslation(widget, m_clickedAtom, m_clickedAtom->pos());
      }
      else if(m_midButtonPressed)
      {
        m_eyecandy->drawZoom(widget, m_clickedAtom, m_clickedAtom->pos());
      }
      else if(m_rightButtonPressed && selectedSize)
      {
        m_eyecandy->drawRotation(widget, m_clickedAtom,
            m_xAngleEyecandy, m_yAngleEyecandy, m_clickedAtom->pos());
      }
    }
    else if(selectedSize)
    {
      if(m_leftButtonPressed)
      {
        m_eyecandy->drawTranslation(widget, m_selectedPrimitivesCenter, 1.5, 0.);
      }
      else if(m_midButtonPressed)
      {
        m_eyecandy->drawZoom(widget, m_selectedPrimitivesCenter, 1.5);
      }
      else if(m_rightButtonPressed)
      {
        m_eyecandy->drawRotation(widget, m_selectedPrimitivesCenter, 3.,
            m_xAngleEyecandy, m_yAngleEyecandy);
      }
    }
    */
    return true;
  }

  int ApplyForceTool::usefulness() const
  {
    return 70000;
  }
}

#include "applyforcetool.moc"

Q_EXPORT_PLUGIN2(applyforcetool, Avogadro::ApplyForceToolFactory)
