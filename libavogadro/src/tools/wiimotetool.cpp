/**********************************************************************
  WiiMoteTool - Manipulation Tool using WiiMote for Avogadro

  Copyright (C) 2007 by Shahzad Ali

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

#include "wiimotetool.h"
#include "quaternion.h"
#include "navigate.h"

#ifdef WIN32
#include <float.h>
#include <math.h>
#define isnan(x) _isnan(x)
#endif

#include <avogadro/primitive.h>
#include <avogadro/color.h>
#include <avogadro/glwidget.h>
#include <avogadro/camera.h>
#include <avogadro/toolgroup.h>

#include <openbabel/obiter.h>
#include <openbabel/mol.h>

#include <QtPlugin>
#include <QString>

#include <QMessageBox>


#include <fcntl.h>
#include <sys/ioctl.h>
#include <sys/types.h>
#include <unistd.h>

#include <linux/input.h>
#include <linux/uinput.h>


#include <bluetooth/bluetooth.h>
#include "cwiid.h"

/**
#include "conf.h"
#include "util.h"
#include "wmplugin.h"
#include "c_plugin.h"
#include "py_plugin.h"
**/

using namespace std;
using namespace OpenBabel;
using namespace Avogadro;
using namespace Eigen;

// ############################ WiiMoteTool ################################
//Globals
struct stick {
        char valid;
        uint8_t x;
        uint8_t y;
        uint8_t max;
};

//cwiid_wiimote_t *wiimote = NULL;
//cwiid_mesg_callback_t cwiid_callback;

bdaddr_t bdaddr;
struct acc_cal wm_cal, nc_cal;
struct cwiid_ir_mesg ir_data;
struct stick nc_stick;
struct stick cc_l_stick, cc_r_stick;
//struct conf conf;

// ##########  Constructor  ##########

WiiMoteTool::WiiMoteTool(QObject *parent) : Tool(parent),
                                                    m_settingsWidget(NULL),
                                                    m_clickedAtom(NULL),
                                                    m_clickedBond(NULL),
                                                    m_selectedBond(NULL),
                                                    m_skeleton(NULL),
                                                    m_referencePoint(NULL),
                                                    m_currentReference(NULL),
                                                    m_snapped(false),
                                                    m_toolGroup(NULL),
                                                    m_leftButtonPressed(false),
                                                    m_midButtonPressed(false),
                                                    m_rightButtonPressed(false),
                                                   
m_movedSinceButtonPressed(false),
                                                    m_showAngles(true),
                                                    m_snapToEnabled(true),
                                                    m_snapToAngle(10),
                                                    m_detached( false )
{
  QAction *action = activateAction();
  action->setIcon(QIcon(QString::fromUtf8(":/wiimote/wiimote.png")));
  action->setToolTip(tr("Manipulation using WiiMote Tool\n\n"
        "Left Mouse:   Click and drag to rotate the view\n"
        "Middle Mouse: Click and drag to zoom in or out\n"
        "Right Mouse:  Click and drag to move the view\n\n"
        "Left Click & drag on a Bond to set the Manipulation Plane:\n"
        "- Left Click & Drag one of the Atoms in the Bond to change the angle\n"
        "- Right Click & Drag one of the Atoms in the Bond to change the length"));
  //action->setShortcut(Qt::Key_F9);
  //&Avogadro::WiiMoteTool::cwiid_callback_ = &Avogadro::WiiMoteTool::cwiid_callback;
  //wiimote = NULL;
}

// ##########  Desctructor  ##########

WiiMoteTool::~WiiMoteTool()
{
  delete m_referencePoint;
  m_referencePoint = NULL;
  delete m_currentReference;
  m_currentReference = NULL;

  if ( !m_detached ) {
    if ( m_uinputThread->isRunning() ) {
      m_uinputThread->stop();
      m_uinputThread->wait();
    }
    delete m_uinputThread;
  }
  
  if (m_settingsWidget)
  {
    m_snapToAngleLabel->deleteLater();
    m_spacer->deleteLater();
    m_showAnglesBox->deleteLater();
    m_snapToCheckBox->deleteLater();
    m_snapToAngleBox->deleteLater();
    m_layout->deleteLater();

    m_settingsWidget->deleteLater();
  }
}

// ##########  clearData  ##########

void WiiMoteTool::clearData()
{
  m_clickedAtom = NULL;
  m_clickedBond = NULL;
  m_selectedBond = NULL;
  delete m_referencePoint;
  m_referencePoint = NULL;
  delete m_currentReference;
  m_currentReference = NULL;
  m_toolGroup = NULL;
  m_leftButtonPressed = false;
  m_midButtonPressed = false;
  m_rightButtonPressed = false;
  m_movedSinceButtonPressed = false;
  m_snapped = false;
}

// ##########  moleculeChanged  ##########

void WiiMoteTool::moleculeChanged(Molecule* previous, Molecule* next)
{
  if (previous) {
    disconnect(previous, 0 , this, 0);
  }

  if (next) {
    connect((Primitive*)next, SIGNAL(primitiveRemoved(Primitive*)), this,
            SLOT(primitiveRemoved(Primitive*)));
  }

  clearData();
}

// ##########  primitiveRemoved  ##########

void WiiMoteTool::primitiveRemoved(Primitive *primitive)
{
  if (primitive == m_clickedAtom || primitive == m_clickedBond ||
      primitive == m_selectedBond) {
    clearData();
  }
}

// ##########  connectToolGroup  ##########

void WiiMoteTool::connectToolGroup(GLWidget *widget, ToolGroup *toolGroup)
{
  if(widget->toolGroup() != toolGroup && widget->toolGroup())
  {
    disconnect(widget->toolGroup(), 0, this, 0);
    connect(widget->toolGroup(), SIGNAL(toolActivated(Tool*)),
          this, SLOT(toolChanged(Tool*)));
    toolGroup = widget->toolGroup();
  }
}

// ##########  toolChanged  ##########

void WiiMoteTool::toolChanged(Tool* tool)
{
  if(tool != this && m_glwidget)
  {
    m_glwidget->update();
    clearData();
  }
}

// ##########  usefulness  ##########

int WiiMoteTool::usefulness() const
{
  return 2000000;
}

// ##########  mousePress  ##########

QUndoCommand* WiiMoteTool::mousePress(GLWidget *widget, const QMouseEvent
*event)
{
  if(m_glwidget != widget)
  {
    disconnect(widget, 0 , this, 0);
    connect(widget, SIGNAL(moleculeChanged(Molecule*, Molecule*)), this,
            SLOT(moleculeChanged(Molecule*, Molecule*)));
    m_glwidget = widget;
    moleculeChanged(NULL, m_glwidget->molecule());
    connectToolGroup(widget, m_toolGroup);
  }

  m_undo = 0;

  m_lastDraggingPosition = event->pos();
  m_movedSinceButtonPressed = false;

#ifdef Q_WS_MAC
  m_leftButtonPressed = (event->buttons() & Qt::LeftButton
                         && event->modifiers() == Qt::NoModifier);
  // On the Mac, either use a three-button mouse
  // or hold down the Option key (AltModifier in Qt notation)
  m_midButtonPressed = ((event->buttons() & Qt::MidButton) ||
                        (event->buttons() & Qt::LeftButton && event->modifiers()
                        & Qt::AltModifier));
  // Hold down the Command key (ControlModifier in Qt notation) for right button
  m_rightButtonPressed = ((event->buttons() & Qt::RightButton) ||
                          (event->buttons() & Qt::LeftButton &&
event->modifiers()
                          & Qt::ControlModifier));
#else
  m_leftButtonPressed = (event->buttons() & Qt::LeftButton);
  m_midButtonPressed = (event->buttons() & Qt::MidButton);
  m_rightButtonPressed = (event->buttons() & Qt::RightButton);
#endif

  m_clickedAtom = NULL;
  m_clickedBond = NULL;

  int oldName = m_selectedBond ? m_selectedBond->GetIdx() : -1;

  // Check if the mouse clicked on any Atoms or Bonds.
  Primitive *clickedPrim = m_glwidget->computeClickedPrimitive(event->pos());

  if (clickedPrim && clickedPrim->type() == Primitive::AtomType)
  {
    //Create an undo instance for this manipulation
    m_undo = new WiiMoteMoveCommand(m_glwidget->molecule());
    // Atom clicked on.
    m_clickedAtom = (Atom*)clickedPrim;

    if ((m_rightButtonPressed || m_leftButtonPressed) &&
isAtomInBond(m_clickedAtom, m_selectedBond))
    {
      m_skeleton = new SkeletonTree();
      m_skeleton->populate(m_clickedAtom, m_selectedBond,
m_glwidget->molecule());
    }
  }
  else if (clickedPrim && clickedPrim->type() == Primitive::BondType)
  {
    // Bond clicked on.
    m_clickedBond = (Bond*)clickedPrim;

    // If the Bond was clicked on with the left mouse button, set it as the
    // currently selected bond and reset the reference point (if the Bond has
    // changed).
    if (m_leftButtonPressed)
    {
      m_selectedBond = m_clickedBond;

      if ((int)m_selectedBond->GetIdx() != oldName)
      {
        delete m_referencePoint;
        m_referencePoint = NULL;

        delete m_currentReference;
        m_currentReference = NULL;

        m_snapped = false;

        Atom *leftAtom = static_cast<Atom*>(m_selectedBond->GetBeginAtom());
        Atom *rightAtom = static_cast<Atom*>(m_selectedBond->GetEndAtom());

        Vector3d left = leftAtom->pos();
        Vector3d right = rightAtom->pos();
        Vector3d leftToRight = right - left;

        Vector3d x = Vector3d(1, 0, 0);
        Vector3d y = Vector3d(0, 1, 0);

        Vector3d A = leftToRight.cross(x);
        Vector3d B = leftToRight.cross(y);

        m_referencePoint = A.norm() >= B.norm() ? new Vector3d(A) : new
Vector3d(B);
        *m_referencePoint = m_referencePoint->normalized();

        Vector3d *reference = calculateSnapTo(widget, m_selectedBond,
                                              m_referencePoint, m_snapToAngle);

        if (reference && m_snapToEnabled)
        {
          m_snapped = true;
          m_currentReference = reference;
          *m_currentReference = m_currentReference->normalized();
        }
        else {
          m_currentReference = new Vector3d(*m_referencePoint);
        }
      }
    }
    /*else if (m_rightButtonPressed)
    {
      m_selectedBond = m_clickedBond;

      delete m_referencePoint;
      m_referencePoint = NULL;

      delete m_currentReference;
      m_currentReference = NULL;
    }*/
  }

  m_glwidget->update();
  return 0;
}

// ##########  mouseRelease  ##########

QUndoCommand* WiiMoteTool::mouseRelease(GLWidget *widget, const QMouseEvent*)
{
  if (!m_clickedAtom && !m_clickedBond && !m_movedSinceButtonPressed)
  {
    delete m_referencePoint;
    m_referencePoint = NULL;
    delete m_currentReference;
    m_currentReference = NULL;
    m_snapped = false;
    m_selectedBond = NULL;
  }
//  else if (!m_clickedAtom && m_clickedBond && !m_movedSinceButtonPressed) {
  else if (!m_movedSinceButtonPressed) {
    m_undo = 0;
  }
  if (m_skeleton)
  {
    delete m_skeleton;
    m_skeleton = NULL;
  }

  m_glwidget = widget;
  m_leftButtonPressed = false;
  m_midButtonPressed = false;
  m_rightButtonPressed = false;
  m_clickedAtom = NULL;
  m_clickedBond = NULL;

  m_glwidget->update();
  return m_undo;
}

// ##########  mouseMove  ##########

QUndoCommand* WiiMoteTool::mouseMove(GLWidget *widget, const QMouseEvent *event)
{
  m_glwidget = widget;
  if (!m_glwidget->molecule()) {
    return 0;
  }

  QPoint deltaDragging = event->pos() - m_lastDraggingPosition;

  if ((event->pos() - m_lastDraggingPosition).manhattanLength() > 2) {
    m_movedSinceButtonPressed = true;
  }

  // Mouse navigation has two modes - atom centred when an atom is clicked
  // and scene if no atom has been clicked.

#ifdef Q_WS_MAC
  if (event->buttons() & Qt::LeftButton && event->modifiers() == Qt::NoModifier)
#else
  if (event->buttons() & Qt::LeftButton)
#endif
  {
    if (m_clickedBond && m_selectedBond && m_referencePoint)
    {
      Atom *beginAtom = static_cast<Atom*>(m_selectedBond->GetBeginAtom());
      Atom *endAtom = static_cast<Atom*>(m_selectedBond->GetEndAtom());

      Vector3d rotationVector = beginAtom->pos() - endAtom->pos();
      rotationVector = rotationVector / rotationVector.norm();

      Vector3d begin = widget->camera()->project(beginAtom->pos());
      Vector3d end = widget->camera()->project(endAtom->pos());

      Vector3d zAxis = Vector3d(0, 0, 1);
      Vector3d beginToEnd = end - begin;
      beginToEnd -= Vector3d(0, 0, beginToEnd.z());

      Vector3d direction = zAxis.cross(beginToEnd);
      direction = direction / direction.norm();

      Vector3d mouseMoved = Vector3d( - deltaDragging.x(), - deltaDragging.y(),
0);

      double magnitude = mouseMoved.dot(direction) / direction.norm();

      *m_referencePoint = performRotation(magnitude * (M_PI / 180.0),
                                          rotationVector, Vector3d(0, 0, 0),
                                          *m_referencePoint);

      Eigen::Vector3d *reference = calculateSnapTo(widget, m_selectedBond,
                                          m_referencePoint, m_snapToAngle);
      if (reference && m_snapToEnabled)
      {
        m_snapped = true;
        delete m_currentReference;
        m_currentReference = reference;
        *m_currentReference = m_currentReference->normalized();
      }
      else
      {
        m_snapped = false;
        delete m_currentReference;
        m_currentReference = new Vector3d(*m_referencePoint);
      }
    }
    else if (isAtomInBond(m_clickedAtom, m_selectedBond))
    {
      //Do atom rotation.
      Atom *otherAtom;

      if (m_clickedAtom == static_cast<Atom*>(m_selectedBond->GetBeginAtom()))
        otherAtom = static_cast<Atom*>(m_selectedBond->GetEndAtom());
      else
        otherAtom = static_cast<Atom*>(m_selectedBond->GetBeginAtom());

      Vector3d center = otherAtom->pos();
      Vector3d clicked = m_clickedAtom->pos();

      Vector3d centerProj = widget->camera()->project(center);
      centerProj -= Vector3d(0,0,centerProj.z());
      Vector3d clickedProj = widget->camera()->project(clicked);
      clickedProj -= Vector3d(0,0,clickedProj.z());
      Vector3d referenceProj = widget->camera()->project(*m_currentReference +
center);
      referenceProj -= Vector3d(0,0,referenceProj.z());

      Vector3d referenceVector = referenceProj - centerProj;
      referenceVector = referenceVector.normalized();

      Vector3d directionVector = clickedProj - centerProj;
      directionVector = directionVector.normalized();

      Vector3d rotationVector = referenceVector.cross(directionVector);
      rotationVector = rotationVector.normalized();

      Vector3d currMouseVector = Vector3d(event->pos().x(), event->pos().y(), 0)
                                  - centerProj;
      if(currMouseVector.norm() > 5)
      {
        currMouseVector = currMouseVector.normalized();
        double mouseAngle = acos(directionVector.dot(currMouseVector) /
                            currMouseVector.norm2());

        if(mouseAngle > 0)
        {
          Vector3d tester;

          tester = performRotation(mouseAngle, rotationVector, Vector3d(0, 0,
0),
                                   directionVector);
          double testAngle1 = acos(tester.dot(currMouseVector) /
                                   currMouseVector.norm2());

          tester = performRotation(-mouseAngle, rotationVector, Vector3d(0, 0,
0),
                                   directionVector);
          double testAngle2 = acos(tester.dot(currMouseVector) /
                                   currMouseVector.norm2());

          if(testAngle1 > testAngle2 || isnan(testAngle2)) {
            mouseAngle = -mouseAngle;
          }

          Vector3d direction = clicked - center;
          if (m_skeleton)
          {
            Vector3d currCrossDir =
m_currentReference->cross(direction).normalized();

            m_skeleton->skeletonRotate(mouseAngle, currCrossDir, center);
            *m_referencePoint = performRotation(mouseAngle, currCrossDir,
                                Vector3d(0, 0, 0), *m_referencePoint);
            *m_currentReference = performRotation(mouseAngle, currCrossDir,
                                  Vector3d(0, 0, 0), *m_currentReference);
          }
        }
      }
    }
    else {
      // rotation around the center of the molecule
      Navigate::rotate(m_glwidget, m_glwidget->center(), deltaDragging.x(),
deltaDragging.y());
    }
  }
#ifdef Q_WS_MAC
  // On the Mac, either use a three-button mouse
  // or hold down the Option key (AltModifier in Qt notation)
  else if ((event->buttons() & Qt::MidButton) || (event->buttons() &
           Qt::LeftButton && event->modifiers() & Qt::AltModifier))
#else
  else if (event->buttons() & Qt::MidButton)
#endif
  {
    if (m_clickedAtom)
    {
      // Perform the rotation
      Navigate::tilt(m_glwidget, m_clickedAtom->pos(), deltaDragging.x());

      // Perform the zoom toward the center of a clicked atom
      Navigate::zoom(m_glwidget, m_clickedAtom->pos(), deltaDragging.y());
    }
    else if (m_clickedBond)
    {
      Atom *begin = static_cast<Atom *>(m_clickedBond->GetBeginAtom());
      Atom *end = static_cast<Atom *>(m_clickedBond->GetEndAtom());

      Vector3d btoe = end->pos() - begin->pos();
      double newLen = btoe.norm() / 2;
      btoe = btoe / btoe.norm();

      Vector3d mid = begin->pos() + btoe * newLen;

      // Perform the rotation
      Navigate::tilt(m_glwidget, mid, deltaDragging.x());

      // Perform the zoom toward the centre of a clicked bond
      Navigate::zoom(m_glwidget, mid, deltaDragging.y());
    }
    else
    {
      // Perform the rotation
      Navigate::tilt(m_glwidget, m_glwidget->center(), deltaDragging.x());

      // Perform the zoom toward molecule center
      Navigate::zoom(m_glwidget, m_glwidget->center(), deltaDragging.y());
    }
  }
#ifdef Q_WS_MAC
  // On the Mac, either use a three-button mouse
  // or hold down the Command key (ControlModifier in Qt notation)
  else if ((event->buttons() & Qt::RightButton) || (event->buttons() &
           Qt::LeftButton && event->modifiers() & Qt::ControlModifier))
#else
  else if (event->buttons() & Qt::RightButton)
#endif
  {
    if (isAtomInBond(m_clickedAtom, m_selectedBond))
    {
      // Adjust the length of the bond following the mouse movement.

      Atom *otherAtom;

      if (m_clickedAtom == static_cast<Atom*>(m_selectedBond->GetBeginAtom()))
        otherAtom = static_cast<Atom*>(m_selectedBond->GetEndAtom());
      else
        otherAtom = static_cast<Atom*>(m_selectedBond->GetBeginAtom());

      Vector3d clicked = m_clickedAtom->pos();
      Vector3d other = otherAtom->pos();
      Vector3d direction = clicked - other;

      Vector3d mouseLast = widget->camera()->unProject(m_lastDraggingPosition);
      Vector3d mouseCurr = widget->camera()->unProject(event->pos());
      Vector3d mouseDir = mouseCurr - mouseLast;

      Vector3d component = mouseDir.dot(direction) / direction.norm2() *
direction;

      if (m_skeleton) {
        m_skeleton->skeletonTranslate(component.x(), component.y(),
component.z());
      }
    }
    else {
      // Translate the molecule following mouse movement.
      Navigate::translate(m_glwidget, m_glwidget->center(),
m_lastDraggingPosition, event->pos());
    }
  }

  m_lastDraggingPosition = event->pos();
  m_glwidget->update();

  return 0;
}

// ##########  wheel  ##########

QUndoCommand* WiiMoteTool::wheel(GLWidget *widget, const QWheelEvent *event)
{
  m_glwidget = widget;

  m_clickedAtom = NULL;
  m_clickedBond = NULL;

  Primitive *clickedPrim = m_glwidget->computeClickedPrimitive(event->pos());

  if (clickedPrim && clickedPrim->type() == Primitive::AtomType)
  {
    Atom *clickedAtom = (Atom*)clickedPrim;
    // Perform the zoom toward clicked atom
    Navigate::zoom(m_glwidget, clickedAtom->pos(), - MOUSE_WHEEL_SPEED *
event->delta());
  }
  else if (clickedPrim && clickedPrim->type() == Primitive::BondType)
  {
    Bond *clickedBond = (Bond*)clickedPrim;

    Atom *begin = static_cast<Atom *>(clickedBond->GetBeginAtom());
    Atom *end = static_cast<Atom *>(clickedBond->GetEndAtom());

    Vector3d btoe = end->pos() - begin->pos();
    double newLen = btoe.norm() / 2;
    btoe = btoe / btoe.norm();

    Vector3d mid = begin->pos() + btoe * newLen;

    // Perform the zoom toward the centre of a clicked bond
    Navigate::zoom(m_glwidget, mid, - MOUSE_WHEEL_SPEED * event->delta());
  }
  else {
    // Perform the zoom toward molecule center
    Navigate::zoom(m_glwidget, m_glwidget->center(), - MOUSE_WHEEL_SPEED *
event->delta());
  }

  m_glwidget->update();

  return 0;
}

// ##########  paint  ##########

bool WiiMoteTool::paint(GLWidget *widget)
{
  if(widget->toolGroup()->activeTool() != this) {
    clearData();
  }

  if ((m_leftButtonPressed && !m_clickedBond && !isAtomInBond(m_clickedAtom,
m_selectedBond))
       || (m_midButtonPressed && !m_clickedBond && !m_clickedAtom)
       || (m_rightButtonPressed && !isAtomInBond(m_clickedAtom,
m_selectedBond))) {
    drawSphere(widget, widget->center(), 0.10, 1.0);
  }

  if (m_leftButtonPressed && m_clickedAtom && (!m_selectedBond ||
      !isAtomInBond(m_clickedAtom, m_selectedBond))) {
    drawAtomAngles(widget, m_clickedAtom);
  }

  if (m_selectedBond)
  {
    Atom *begin = static_cast<Atom*>(m_selectedBond->GetBeginAtom());
    Atom *end = static_cast<Atom*>(m_selectedBond->GetEndAtom());

    if (m_currentReference)
    {
      // Draw bond length text.
      QString length = tr("Bond Length:  ") +
                      QString::number(m_selectedBond->GetLength(), 10, 1) +
                      QString::fromUtf8(" Å (Angstrom)");

      glColor4f(1.0, 1.0, 1.0, 1.0);
      widget->painter()->setColor(1.0, 1.0, 1.0, 1.0);
      widget->painter()->drawText(QPoint(5, widget->height() - 25), length);

      if (m_rightButtonPressed && (m_clickedAtom == begin || m_clickedAtom ==
end)) {
        drawSkeletonAngles(widget, m_skeleton);
      }
      else
      {
        if (m_showAngles)
        {
          // Draw the angles around the two atoms.
          if (!m_clickedAtom || m_rightButtonPressed || m_midButtonPressed ||
              (m_leftButtonPressed && begin != m_clickedAtom)) {
            drawAngles(widget, begin, m_selectedBond);
          }

          if (!m_clickedAtom || m_rightButtonPressed || m_midButtonPressed ||
              (m_leftButtonPressed && end != m_clickedAtom)) {
            drawAngles(widget, end, m_selectedBond);
          }
        }
        else
        {
          // Draw the angles around the two atoms.
          if (m_leftButtonPressed && end == m_clickedAtom) {
            drawAngles(widget, begin, m_selectedBond);
          }

          if (m_leftButtonPressed && begin == m_clickedAtom) {
            drawAngles(widget, end, m_selectedBond);
          }
        }

        if (m_clickedAtom && m_leftButtonPressed &&
            isAtomInBond(m_clickedAtom, m_selectedBond)) {
          drawSkeletonAngles(widget, m_skeleton);
        }
      }

      // Draw the manipulation rectangle.
      if (m_snapped && m_snapToEnabled)
      {
        double rgb[3] = {1.0, 1.0, 0.2};
        drawManipulationRectangle(widget, m_selectedBond, m_currentReference,
rgb);
      }
      else
      {
        double rgb[3] = {0.0, 0.2, 0.8};
        drawManipulationRectangle(widget, m_selectedBond, m_currentReference,
rgb);
      }
    }
    /*else
    {
      OBBondIterator bondIter = begin->EndBonds();
      Atom *v = (Atom*)begin->BeginNbrAtom(bondIter);

      if (v != NULL)
      {
        do
        {
          if (v == end) {
            continue;
          }
          else
          {
            Vector3d tmp1 = end->pos() - begin->pos();
            Vector3d tmp2 = v->pos() - begin->pos();
            double tAngle = acos(tmp1.dot(tmp2) / (tmp2.norm() * tmp1.norm()))
                            * 180.0 / M_PI;

            if(!(tAngle > 1 && tAngle < 179)) {
              continue;
            }
          }

          Vector3d vVec = v->pos() - begin->pos();
          vVec = vVec.normalized();

          Vector3d tmp = end->pos() - begin->pos();
          vVec = begin->pos() + vVec * tmp.norm();

          drawAngleSector(widget, begin->pos(), end->pos(), vVec);
        }
        while ((v = (Atom*)begin->NextNbrAtom(bondIter)) != NULL);
      }

      bondIter = end->EndBonds();
      v = (Atom*)end->BeginNbrAtom(bondIter);

      if (v != NULL)
      {
        do
        {
          if (v == begin) {
            continue;
          }
          else
          {
            Eigen::Vector3d tmp1 = begin->pos() - end->pos();
            Eigen::Vector3d tmp2 = v->pos() - end->pos();
            double tAngle = acos(tmp1.dot(tmp2) / (tmp2.norm() * tmp1.norm()))
                  * 180.0 / M_PI;

            if(!(tAngle > 1 && tAngle < 179)) {
              continue;
            }
          }

          Vector3d vVec = v->pos() - end->pos();
          vVec = vVec.normalized();

          Vector3d tmp = begin->pos() - end->pos();
          vVec = end->pos() + vVec * tmp.norm();

          drawAngleSector(widget, end->pos(), begin->pos(), vVec);
        }
        while ((v = (Atom*)end->NextNbrAtom(bondIter)) != NULL);
      }
    }*/
  }

  return true;
}

// ##########  isAtomInBond  ##########

bool WiiMoteTool::isAtomInBond(Atom *atom, Bond *bond)
{
  if (!atom || !bond) {
    return false;
  }

  if (atom == static_cast<Atom*>(bond->GetBeginAtom())) {
    return true;
  }

  return atom == static_cast<Atom*>(bond->GetEndAtom());
}

// ##########  drawAtomAngles  ##########

void WiiMoteTool::drawAtomAngles(GLWidget *widget, Atom *atom)
{
  if (!atom || !widget) {
    return;
  }

  OBBondIterator bondIter = atom->EndBonds();

  Atom *u = (Atom*)atom->BeginNbrAtom(bondIter);
  Atom *v = NULL;

  if (u != NULL)
  {
    do
    {
      OBBondIterator tmpIter = bondIter;

      while ((v = (Atom*)atom->NextNbrAtom(tmpIter)) != NULL) {
        drawAngleSector(widget, atom->pos(), u->pos(), v->pos());
      }
    }
    while((u = (Atom*)atom->NextNbrAtom(bondIter)) != NULL);
  }
}

// ##########  drawSkeletonAngles  ##########

void WiiMoteTool::drawSkeletonAngles(GLWidget *widget, SkeletonTree *skeleton)
{
  if (!skeleton || !widget) {
    return;
  }

  Atom *atom = skeleton->rootAtom();
  Bond *bond = skeleton->rootBond();

  Atom *ref = NULL;
  if (atom == static_cast<Atom*>(bond->GetBeginAtom())) {
    ref = static_cast<Atom*>(bond->GetEndAtom());
  }
  else if (atom == static_cast<Atom*>(bond->GetEndAtom())) {
    ref = static_cast<Atom*>(bond->GetBeginAtom());
  }
  else {
    return;
  }

  OBBondIterator bondIter = atom->EndBonds();
  Atom *v = (Atom*)atom->BeginNbrAtom(bondIter);

  if (v != NULL)
  {
    do
    {
      if (v == ref) {
        continue;
      }

      if (!skeleton->containsAtom(v)) {
        drawAngleSector(widget, atom->pos(), ref->pos(), v->pos());
      }
    }
    while ((v = (Atom*)atom->NextNbrAtom(bondIter)) != NULL);
  }
}

// ##########  drawAngles  ##########

void WiiMoteTool::drawAngles(GLWidget *widget, Atom *atom, Bond *bond)
{
  if (!atom || !bond || !widget) {
    return;
  }

  assert(isAtomInBond(atom, bond));

  Atom *ref = NULL;
  if (atom == static_cast<Atom*>(bond->GetBeginAtom())) {
    ref = static_cast<Atom*>(bond->GetEndAtom());
  }
  else if (atom == static_cast<Atom*>(bond->GetEndAtom())) {
    ref = static_cast<Atom*>(bond->GetBeginAtom());
  }
  else {
    return;
  }

  OBBondIterator bondIter = atom->EndBonds();
  Atom *v = (Atom*)atom->BeginNbrAtom(bondIter);

  if (v != NULL)
  {
    do
    {
      if (v == ref) {
        continue;
      }

      drawAngleSector(widget, atom->pos(), ref->pos(), v->pos());
    }
    while ((v = (Atom*)atom->NextNbrAtom(bondIter)) != NULL);
  }
}

// ##########  drawAngleSector  ##########

void WiiMoteTool::drawAngleSector(GLWidget *widget, Eigen::Vector3d origin,
                                      Eigen::Vector3d direction1,
Eigen::Vector3d direction2)
{
  // Get vectors representing the lines from centre to left and centre to right.
  Eigen::Vector3d u = direction1 - origin;
  Eigen::Vector3d v = direction2 - origin;

  // Calculate the length of the vectors (half the length of the shortest vector.)
  double radius = qMin(u.norm(), v.norm()) * 0.5;
  double lineWidth = 1.5;

  // Adjust the length of u and v to the length calculated above.
  u = (u / u.norm()) * radius;
  v = (v / v.norm()) * radius;

  // Angle between u and v.
  double uvAngle = acos(u.dot(v) / v.norm2()) * 180.0 / M_PI;

  // If angle is less than 1 (will be approximated to 0), attempting to draw
  // will crash, so return.
  if (abs(uvAngle) <= 1) {
    return;
  }

  // Vector perpindicular to both u and v.
  Eigen::Vector3d n = u.cross(v);

  Eigen::Vector3d x = Vector3d(1, 0, 0);
  Eigen::Vector3d y = Vector3d(0, 1, 0);

  if (n.norm() < 1e-16)
  {
    Eigen::Vector3d A = u.cross(x);
    Eigen::Vector3d B = u.cross(y);

    n = A.norm() >= B.norm() ? A : B;
  }

  n = n / n.norm();

  Vector3d point = performRotation((uvAngle / 2 * (M_PI / 180.0)), n,
                                    Vector3d(0, 0, 0), u);

  QString angle = QString::number(uvAngle, 10, 1) + QString::fromUtf8("°");
  glColor4f(1.0, 1.0, 1.0, 1.0);
  widget->painter()->setColor(1.0, 1.0, 1.0, 1.0);
  widget->painter()->drawText(point + origin, angle);

  glEnable(GL_BLEND);
  widget->painter()->setColor(0, 0.5, 0, 0.4);
  glDepthMask(GL_FALSE);
  widget->painter()->drawShadedSector(origin, direction1, direction2, radius);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);

  widget->painter()->setColor(1.0, 1.0, 1.0, 1.0);
  widget->painter()->drawArc(origin, direction1, direction2, radius, lineWidth);
}

// ##########  calcualteSnapTo  ##########

Eigen::Vector3d* WiiMoteTool::calculateSnapTo(GLWidget *widget, Bond *bond,
    Eigen::Vector3d *referencePoint, double maximumAngle)
{
  if(!referencePoint || !bond || !widget) {
    return NULL;
  }

  double angle = -1;
  Eigen::Vector3d *smallestRef = NULL;
  Atom *b = static_cast<Atom*>(bond->GetBeginAtom());
  Atom *e = static_cast<Atom*>(bond->GetEndAtom());

  OBBondIterator bondIter = b->EndBonds();
  Atom *t = (Atom*)b->BeginNbrAtom(bondIter);

  Eigen::Vector3d begin = b->pos();
  Eigen::Vector3d end = e->pos();
  Eigen::Vector3d target;

  if (t != NULL)
  {
    do
    {
      if (t == e) {
        continue;
      }

      target = t->pos();

      Eigen::Vector3d u = end - begin;
      Eigen::Vector3d v = target - begin;
      double tAngle = acos(u.dot(v) / (v.norm() * u.norm())) * 180.0 / M_PI;

      if(!(tAngle > 1 && tAngle < 179)) {
        continue;
      }

      Eigen::Vector3d orth1 = u.cross(v);
      Eigen::Vector3d orth2 = referencePoint->cross(u);

      tAngle = acos(orth1.dot(orth2) / (orth1.norm() * orth2.norm())) * 180.0 /
M_PI;
      tAngle = tAngle > 90 ? 180 - tAngle : tAngle;

      if(angle < 0)
      {
        angle = tAngle;
        smallestRef = new Vector3d(v);
      }
      else if(tAngle < angle)
      {
        angle = tAngle;
        delete smallestRef;
        smallestRef = new Vector3d(v);
      }
    }
    while ((t = (Atom*)b->NextNbrAtom(bondIter)) != NULL);
  }

  bondIter = e->EndBonds();
  t = (Atom*)e->BeginNbrAtom(bondIter);

  if (t != NULL)
  {
    do
    {
      if (t == b) {
        continue;
      }

      target = t->pos();

      Eigen::Vector3d u = begin - end;
      Eigen::Vector3d v = target - end;
      double tAngle = acos(u.dot(v) / (v.norm() * u.norm())) * 180.0 / M_PI;

      if(!(tAngle > 1 && tAngle < 179)) {
        continue;
      }

      Eigen::Vector3d orth1 = u.cross(v);
      Eigen::Vector3d orth2 = referencePoint->cross(u);

      tAngle = acos(orth1.dot(orth2) / (orth1.norm() * orth2.norm())) * 180.0 /
M_PI;
      tAngle = tAngle > 90 ? 180 - tAngle : tAngle;

      if(angle < 0)
      {
        angle = tAngle;
        smallestRef = new Vector3d(v);
      }
      else if(tAngle < angle)
      {
        angle = tAngle;
        delete smallestRef;
        smallestRef = new Vector3d(v);
      }
    }
    while ((t = (Atom*)e->NextNbrAtom(bondIter)) != NULL);
  }

  if (angle > maximumAngle)
  {
    if (smallestRef) {
      delete smallestRef;
    }

    return NULL;
  }

  return smallestRef;
}

// ##########  drawManipulationRectangle  ##########

void WiiMoteTool::drawManipulationRectangle(GLWidget *widget, Bond *bond,
    Eigen::Vector3d *referencePoint, double rgb[3])
{
  if (!bond || !widget || !referencePoint) {
    return;
  }

  Atom *leftAtom = static_cast<Atom*>(bond->GetBeginAtom());
  Atom *rightAtom = static_cast<Atom*>(bond->GetEndAtom());

  Eigen::Vector3d left = leftAtom->pos();
  Eigen::Vector3d right = rightAtom->pos();

  Eigen::Vector3d leftToRight = right - left;

  Eigen::Vector3d vec = leftToRight.cross(*referencePoint);
  Eigen::Vector3d planeVec = vec.cross(leftToRight);

  double length = 1;

  planeVec = length * (planeVec / planeVec.norm());

  Eigen::Vector3d topLeft = widget->camera()->modelview() * (left + planeVec);
  Eigen::Vector3d topRight = widget->camera()->modelview() * (right + planeVec);
  Eigen::Vector3d botRight = widget->camera()->modelview() * (right - planeVec);
  Eigen::Vector3d botLeft = widget->camera()->modelview() * (left - planeVec);

  float alpha = 0.4;
  double lineWidth = 1.5;

  glEnable(GL_BLEND);
  widget->painter()->setColor(rgb[0], rgb[1], rgb[2], alpha);
  glDepthMask(GL_FALSE);
  widget->painter()->drawShadedQuadrilateral(topLeft, topRight, botRight,
botLeft);
  glDepthMask(GL_TRUE);
  glDisable(GL_BLEND);
  widget->painter()->setColor(1.0, 1.0, 1.0, 1.0);
  widget->painter()->drawQuadrilateral(topLeft, topRight, botRight, botLeft,
lineWidth);
}

// ##########  drawSphere  ##########

void WiiMoteTool::drawSphere(GLWidget *widget,  const Eigen::Vector3d &position,
                                 double radius, float alpha )
{
  glEnable(GL_BLEND);
  widget->painter()->setColor(1.0, 1.0, 0.3, alpha);
  widget->painter()->drawSphere(position, radius);
  glDisable(GL_BLEND);
}

// ##########  performRotation  ##########

Eigen::Vector3d WiiMoteTool::performRotation(double angle,
    Eigen::Vector3d rotationVector, Eigen::Vector3d centerVector,
    Eigen::Vector3d positionVector)
{
  Quaternion qLeft = Quaternion::createRotationLeftHalf(angle, rotationVector);
  Quaternion qRight = qLeft.multiplicitiveInverse();

  return Quaternion::performRotationMultiplication(qLeft, positionVector -
                     centerVector, qRight) + centerVector;
}

// ##########  showAnglesChanged  ##########

void WiiMoteTool::showAnglesChanged(int state)
{
  m_showAngles = state == Qt::Checked ? true : false;

  if (m_glwidget) {
    m_glwidget->update();
  }
}

// ##########  snapToCheckBoxChanged  ##########

void WiiMoteTool::snapToCheckBoxChanged(int state)
{
  m_snapToEnabled = state == Qt::Checked ? true : false;
  m_snapToAngleBox->setEnabled(m_snapToEnabled);

  if(!m_selectedBond) {
    return;
  }

  Eigen::Vector3d *reference = calculateSnapTo(m_glwidget, m_selectedBond,
                                               m_referencePoint, m_snapToAngle);
  if (reference && m_snapToEnabled)
  {
    m_snapped = true;
    delete m_currentReference;
    m_currentReference = reference;
    *m_currentReference = m_currentReference->normalized();
  }
  else
  {
    m_snapped = false;
    delete m_currentReference;
    m_currentReference = new Vector3d(*m_referencePoint);
  }

  if (m_glwidget) {
    m_glwidget->update();
  }
}

// ##########  snapToAngleChanged  ##########

void WiiMoteTool::snapToAngleChanged(int newAngle)
{
  m_snapToAngle = newAngle;

  if(!m_selectedBond) {
    return;
  }

  Eigen::Vector3d *reference = calculateSnapTo(m_glwidget, m_selectedBond,
                                               m_referencePoint, m_snapToAngle);
  if (reference && m_snapToEnabled)
  {
    m_snapped = true;
    delete m_currentReference;
    m_currentReference = reference;
    *m_currentReference = m_currentReference->normalized();
  }
  else
  {
    m_snapped = false;
    delete m_currentReference;
    m_currentReference = new Vector3d(*m_referencePoint);
  }

  if (m_glwidget) {
    m_glwidget->update();
  }
}

void WiiMoteTool::cwiid_callback(cwiid_wiimote_t *wiimote, int mesg_count,
                    union cwiid_mesg mesg[], struct timespec *timestamp)
{
  int i;

  for (i=0; i < mesg_count; i++) {
    switch (mesg[i].type) {
      case CWIID_MESG_BTN:
        send_btn_event((struct cwiid_btn_mesg *) &mesg[i]);
        //cout << "Msg Btn" << endl;
        break;
      case CWIID_MESG_ACC:
        //cout << "ACC: " << &mesg[i].acc_mesg << endl;
       // cwiid_acc(&mesg_array[i].acc_mesg);
        break;
      case CWIID_MESG_IR:
        //cout << "IR: " << &mesg[i].ir_mesg << endl;
       // cwiid_ir(&mesg_array[i].ir_mesg);
        break;
      case CWIID_MESG_NUNCHUK:
        //process_nunchuk_mesg((struct cwiid_nunchuk_mesg *) &mesg[i]);
        cout << "Nun Chuck" << endl;
        break;
      case CWIID_MESG_CLASSIC:
        //process_classic_mesg((struct cwiid_classic_mesg *) &mesg[i]);
        cout << "Classic" << endl;
        break;
      case CWIID_MESG_ERROR:
        cout << "Received Error, Disconnecting" << endl;
        //Send disconnect event
        break;
      default:
        break;
    }
  }
  /**for (i=0; (i < CONF_MAX_PLUGINS) && conf.plugins[i].name; i++) {
    process_plugin(&conf.plugins[i], mesg_count, mesg);
  }**/
  send_event(m_conf, EV_SYN, SYN_REPORT, 0);
}

void WiiMoteTool::send_btn_event(struct cwiid_btn_mesg *mesg)
{
  static uint16_t prev_buttons = 0;
  uint16_t pressed, released;
  __s32 axis_value;
  int i;

  /* Wiimote Button/Key Events */
  pressed = mesg->buttons & ~prev_buttons;
  released = ~mesg->buttons & prev_buttons;
  for (i=0; i < CONF_WM_BTN_COUNT; i++) {
    if (m_conf.wiimote_bmap[i].active) {
      if (pressed & m_conf.wiimote_bmap[i].mask) {
        send_event(m_conf, EV_KEY, m_conf.wiimote_bmap[i].action, 1);
      }
      else if (released & m_conf.wiimote_bmap[i].mask) {
        send_event(m_conf, EV_KEY, m_conf.wiimote_bmap[i].action, 0);
      }
    }
  }
  prev_buttons = mesg->buttons;

  /* Wiimote.Dpad.X 
  if (conf.amap[CONF_WM_AXIS_DPAD_X].active) {
    axis_value = 0;
    if (mesg->buttons & CWIID_BTN_LEFT) {
      axis_value = -1;
    }
    else if (mesg->buttons & CWIID_BTN_RIGHT) {
      axis_value = 1;
    }
    if (conf.amap[CONF_WM_AXIS_DPAD_X].flags & CONF_INVERT) {
      axis_value *= -1;
    }
    send_event(&conf, conf.amap[CONF_WM_AXIS_DPAD_X].axis_type,
                conf.amap[CONF_WM_AXIS_DPAD_X].action, axis_value);
  }

  /* Wiimote.Dpad.Y 
  if (conf.amap[CONF_WM_AXIS_DPAD_Y].active) {
    axis_value = 0;
    if (mesg->buttons & CWIID_BTN_DOWN) {
      axis_value = -1;
    }
    else if (mesg->buttons & CWIID_BTN_UP) {
      axis_value = 1;
    }
    if (conf.amap[CONF_WM_AXIS_DPAD_Y].flags & CONF_INVERT) {
      axis_value *= -1;
    }
    send_event(&conf, conf.amap[CONF_WM_AXIS_DPAD_Y].axis_type,
                conf.amap[CONF_WM_AXIS_DPAD_Y].action, axis_value);
  }**/
}


int WiiMoteTool::send_event(struct conf_st conf, __u16 type, __u16 code, __s32 value)
{
  struct input_event event;

  memset(&event, 0, sizeof(event));
  event.type = type;
  event.code = code;
  event.value = value;

  int size = write(conf.fd, &event, sizeof(event));
  bool error = size != sizeof(event);
  //cout << type << " " << code << " " << value << " #" << size << " X" << conf.fd << endl;
  if (error) {
    cout << "Error on send_event" << endl;
    return -1;
  }
  return 0;
}


void* ptO2Object;
void WiiMoteTool::cwiid_callback_wrapper(cwiid_wiimote_t *wiimote, int mesg_count,
                                 union cwiid_mesg mesg[], struct timespec *timestamp)
{
  WiiMoteTool* thisTool = (WiiMoteTool*) ptO2Object;
  thisTool->cwiid_callback(wiimote, mesg_count, mesg,timestamp);
}

void WiiMoteTool::connectClicked()
{
  cout << "Connect WiiMote Clicked" << endl;

  if ( QMessageBox::information(NULL, "Wii-Remote -- Ready?",
       "Put Wiimote in discoverable mode (Press 1+2) and press Ok.",
       tr("&Ok"), tr("&Cancel"), QString::null, 0 , 1) )
  {
  } else 
  {
    ptO2Object = this;
    char reset_bdaddr = 0;
    if (bacmp(&bdaddr, BDADDR_ANY) == 0) {
      reset_bdaddr = 1;
    }

    if ((wiimote = cwiid_open(&bdaddr, CWIID_FLAG_MESG_IFC)) == NULL) {
      QMessageBox::information(NULL, "Wii-Remote",
                               "Unable to connect to Wii-Remote.\nMake sure it is in discoverable mode.");
    } 
    //else if (cwiid_set_mesg_callback(wiimote, &Avogadro::WiiMoteTool::cwiid_callback)) {
    else if (cwiid_set_mesg_callback(wiimote, &Avogadro::WiiMoteTool::cwiid_callback_wrapper)) {
      QMessageBox::information( NULL, "Wii-Remote",
    "Error setting callback.");

    if (cwiid_close(wiimote)) {
    QMessageBox::information( NULL, "Wii-Remote",
    "Error on disconnect.");
    }
    wiimote = NULL;
    }
    else {
      //connected
      ui.m_buttonConnect->setText("Connected");
      ui.m_buttonConnect->setEnabled(false);
      ui.m_buttonDisconnect->setEnabled(true);
      if (cwiid_get_acc_cal(wiimote, CWIID_EXT_NONE, &wm_cal)) {
        QMessageBox::information( NULL, "Wii-Remote",
                                  "Unable to retrieve accelerometer.\n");
      }
      ///// Thread for uinput 
      //m_uinputThread = new WiiMoteUInputThread(wiimote);
      //m_uinputThread->start();

      load_conf();

      ///////////////////////////////////////////////////////////////
      /* UInput */
      char *uinput_filename[] = {"/dev/uinput", "/dev/input/uinput",
        "/dev/misc/uinput"};
        int uinputfilenamecount = 3;

        int i;
        int j;

        /* Open uinput device */
        for (i=0; i < uinputfilenamecount; i++) {
          m_conf.fd = open(uinput_filename[i], O_RDWR);
          if (m_conf.fd >= 0)
            break;
        }

        if (m_conf.fd < 0) {
          cout << "Unable to open uinput" << endl;
        }

        if (write(m_conf.fd, &m_conf.dev, sizeof m_conf.dev) != sizeof m_conf.dev) {
          cout << "error on uinput device setup" << endl;
          close(m_conf.fd);
        }

        if (m_conf.ff) {
          if (ioctl(m_conf.fd, UI_SET_EVBIT, EV_FF) < 0) {
            cout << "error on uinput ioctl" << endl;
            close(m_conf.fd);
          }
          if (ioctl(m_conf.fd, UI_SET_FFBIT, FF_RUMBLE) < 0) {
            cout << "error on uinput ioctl" << endl;
            close(m_conf.fd);
          }
        }
  
        if (ioctl(m_conf.fd, UI_SET_EVBIT, EV_KEY) < 0) {
          cout << "error on uinput ioctl";
          close(m_conf.fd);
        } 

      for (i=0; i < CONF_WM_BTN_COUNT; i++) {
        if (m_conf.wiimote_bmap[i].active) {
          if (ioctl(m_conf.fd, UI_SET_KEYBIT, m_conf.wiimote_bmap[i].action)
          < 0) {
            cout << "error on uinput ioctl";
            close(m_conf.fd);
          }
        } 
      }

      if (ioctl(m_conf.fd, UI_DEV_CREATE) < 0) {
        cout << "Error on uinput dev create";
        close(m_conf.fd);
      }

      set_report_mode();
      cwiid_request_status(wiimote);

    }

   if (reset_bdaddr) { bdaddr = *BDADDR_ANY; }

  } 
}

void WiiMoteTool::load_conf()
{
  memset(&m_conf.dev, 0, sizeof m_conf.dev);
  strncpy(m_conf.dev.name, UINPUT_NAME, UINPUT_MAX_NAME_SIZE);
  m_conf.dev.id.bustype = UINPUT_BUSTYPE;
  m_conf.dev.id.vendor = UINPUT_VENDOR;
  m_conf.dev.id.product = UINPUT_PRODUCT;
  m_conf.dev.id.version = UINPUT_VERSION;
  for (int i=0; i < ABS_MAX; i++) {
    m_conf.dev.absmax[i] = -1;
    m_conf.dev.absmin[i] = -1;
    m_conf.dev.absfuzz[i] = -1;
    m_conf.dev.absflat[i] = -1;
  }
  for (int i=0; i < CONF_WM_BTN_COUNT; i++) {
    m_conf.wiimote_bmap[i].active = 0;
  }
  //Assign key wii mote masks
  m_conf.wiimote_bmap[CONF_WM_BTN_UP].mask = CWIID_BTN_UP;
  m_conf.wiimote_bmap[CONF_WM_BTN_DOWN].mask = CWIID_BTN_DOWN;
  m_conf.wiimote_bmap[CONF_WM_BTN_LEFT].mask = CWIID_BTN_LEFT;
  m_conf.wiimote_bmap[CONF_WM_BTN_RIGHT].mask = CWIID_BTN_RIGHT;
  m_conf.wiimote_bmap[CONF_WM_BTN_A].mask = CWIID_BTN_A;
  m_conf.wiimote_bmap[CONF_WM_BTN_B].mask = CWIID_BTN_B;
  m_conf.wiimote_bmap[CONF_WM_BTN_MINUS].mask = CWIID_BTN_MINUS;
  m_conf.wiimote_bmap[CONF_WM_BTN_PLUS].mask = CWIID_BTN_PLUS;
  m_conf.wiimote_bmap[CONF_WM_BTN_HOME].mask = CWIID_BTN_HOME;
  m_conf.wiimote_bmap[CONF_WM_BTN_1].mask = CWIID_BTN_1;
  m_conf.wiimote_bmap[CONF_WM_BTN_2].mask = CWIID_BTN_2;
  //Assign key codes
  m_conf.wiimote_bmap[CONF_WM_BTN_A].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_A].action = BTN_LEFT;
  m_conf.wiimote_bmap[CONF_WM_BTN_B].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_B].action = BTN_RIGHT;

  m_conf.wiimote_bmap[CONF_WM_BTN_1].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_1].action = KEY_A;
  m_conf.wiimote_bmap[CONF_WM_BTN_2].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_2].action = KEY_B;
}

void WiiMoteTool::set_report_mode()
{
  uint8_t rpt_mode;
  
  rpt_mode = CWIID_RPT_STATUS | CWIID_RPT_BTN;

  //if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(chkIR))) {
    rpt_mode |= CWIID_RPT_IR;
  //}
  //if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(chkAcc))) {
    rpt_mode |= CWIID_RPT_ACC;
  /**}
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(chkExt))) {
    rpt_mode |= CWIID_RPT_EXT;
  }**/
  if (cwiid_set_rpt_mode(wiimote, rpt_mode)) {
    QMessageBox::information( NULL, "Wii-Remote",
                              "Error setting report mode.\n");
  }
}

void WiiMoteTool::disconnectClicked()
{
  if (cwiid_close(wiimote)) {
        QMessageBox::information(NULL, "Wii-Remote",
                                 "Error on disconnect.\n");
  }
  
  cleanup();
  
  wiimote = NULL;
  ui.m_buttonConnect->setText("Connect");
  ui.m_buttonConnect->setEnabled(true);
  ui.m_buttonDisconnect->setEnabled(false);
  //clear_widgets();
  //set_gui_state();
}

void WiiMoteTool::detach() const
{
  m_detached = true;
}

WiiMoteUInputThread *WiiMoteTool::uinputThread() const
{
  return m_uinputThread;
}

void WiiMoteTool::cleanup()
{
  if ( !m_detached ) {
    if ( m_uinputThread->isRunning() ) {
      m_uinputThread->stop();
      m_uinputThread->wait();
    }
    delete m_uinputThread;
  }
}

// ##########  settingsWidget  ##########

QWidget *WiiMoteTool::settingsWidget()
{

  if(!m_settingsWidget) {
    m_settingsWidget = new QWidget;
    
    ui.setupUi(m_settingsWidget);
    
    // Connect the connect button
    connect(ui.m_buttonConnect, SIGNAL(clicked()),
        this, SLOT(connectClicked()));
    
    // Connect the disconnect button
    connect(ui.m_buttonDisconnect, SIGNAL(clicked()),
        this, SLOT(disconnectClicked()));

    connect(m_settingsWidget, SIGNAL(destroyed()),
        this, SLOT(settingsWidgetDestroyed()));
  }

  return m_settingsWidget;
}


// ##########  settingsWidgetDestroyed  ##########

void WiiMoteTool::settingsWidgetDestroyed() {
  m_settingsWidget = 0;
}

// #########################  WiiMoteUInputThread  #########################
WiiMoteUInputThread::WiiMoteUInputThread(cwiid_wiimote_t *wiimote, QObject *parent ) : 
    QThread( parent )
{
  m_stop = false;
  /* UInput */
  char *uinput_filename[] = {"/dev/uinput", "/dev/input/uinput",
    "/dev/misc/uinput"};
  int uinputfilenamecount = 3;

  int i;
  int j;
  int request;
  
  /* Open uinput device */
  for (i=0; i < uinputfilenamecount; i++) {
    conf.fd = open(uinput_filename[i], O_RDWR);
    if (conf.fd >= 0)
      break;
  }
  
  if (conf.fd < 0) {
    cout << "Unable to open uinput" << endl;
  }

  if (write(conf.fd, &conf.dev, sizeof conf.dev) != sizeof conf.dev) {
    cout << "error on uinput device setup" << endl;
    close(conf.fd);
  }
  
  if (conf.ff) {
    if (ioctl(conf.fd, UI_SET_EVBIT, EV_FF) < 0) {
      cout << "error on uinput ioctl" << endl;
      close(conf.fd);
    }
    if (ioctl(conf.fd, UI_SET_FFBIT, FF_RUMBLE) < 0) {
      cout << "error on uinput ioctl" << endl;
      close(conf.fd);
    }
  }
  
  if (ioctl(conf.fd, UI_SET_EVBIT, EV_KEY) < 0) {
    cout << "error on uinput ioctl";
    close(conf.fd);
  } 
  /**
#define CONF_WM_BTN_COUNT 11
  for (i=0; i < CONF_WM_BTN_COUNT; i++) {
    //if (conf->wiimote_bmap[i].active) {
      if (ioctl(conf.fd, UI_SET_KEYBIT, conf.wiimote_bmap[i].action)
          < 0) {
        wminput_err("error on uinput ioctl");
        close(conf->fd);
        return -1;
          }
    //}     
  }**/
  
  
  if (ioctl(conf.fd, UI_DEV_CREATE) < 0) {
    cout << "Error on uinput dev create";
    close(conf.fd);
  }


}

int WiiMoteUInputThread::send_event(__u16 type, __u16 code, __s32 value)
{
  struct input_event event;

  memset(&event, 0, sizeof(event));
  event.type = type;
  event.code = code;
  event.value = value;
  
  if (write(conf.fd, &event, sizeof(event)) != sizeof(event)) {
    cout << "Error on send_event" << endl;
    return -1;
  }

  return 0;
}

void WiiMoteUInputThread::run()
{
  /**
  //thread body
  size_t len;
  struct input_event event;
  struct uinput_ff_upload upload;
  struct uinput_ff_erase erase;

  do {
    if ((len = read(conf.fd, &m_event, sizeof m_event)) !=
         sizeof m_event) {
      cout << "Error on WiiMoteUInputThread read" << endl;
      continue;
         }

         switch (m_event.type) {
           case EV_UINPUT:
             switch (m_event.code) {
               case UI_FF_UPLOAD:
                 erase.request_id = m_event.value;
                 if (ioctl(conf.fd, UI_BEGIN_FF_UPLOAD, &upload) < 0) {
                   cout << "Error on ff upload begin" << endl;
                 }
                 if (cwiid_set_rumble(m_wiimote, 1)) {
                   cout << "Error setting rumble" << endl;
                 }
                 if (ioctl(conf.fd, UI_END_FF_UPLOAD, &upload) < 0) {
                   cout << "Error on ff upload end" << endl;
                 }
                 break;
               case UI_FF_ERASE:
                 erase.request_id = event.value;
                 if (ioctl(conf.fd, UI_BEGIN_FF_ERASE, &erase) < 0) {
                   cout << "Error on ff erase begin" << endl;
                 }
                 if (cwiid_set_rumble(m_wiimote, 0)) {
                   cout << "Error clearing rumble" << endl;
                 }
                 if (ioctl(conf.fd, UI_END_FF_ERASE, &erase) < 0) {
                   cout << "Error on ff erase end" << endl;
                 }
                 break;
               default:
                 break;
             }
             break;
           default:
             break;
         }
  } while (-1);
  //emit stepsTaken(..) to emit a signal
  **/
}

void WiiMoteUInputThread::setEvent(input_event event)
{
  m_event = event;
}

void WiiMoteUInputThread::stop()
{
  m_stop = true;
}
// #########################  WiiMoteMoveCommand  ##########################

// ##########  Constructor  ##########

WiiMoteMoveCommand::WiiMoteMoveCommand(Molecule *molecule,
                                               QUndoCommand *parent)
                      : QUndoCommand(parent), m_molecule(0)
{
  // Store the molecule - this call won't actually move an atom
  setText(QObject::tr("Bond Centric Manipulation"));
  m_moleculeCopy = *molecule;
  m_molecule = molecule;
  m_atomIndex = 0;
  undone = false;
}

// ##########  Constructor  ##########

WiiMoteMoveCommand::WiiMoteMoveCommand(Molecule *molecule,
                                               Atom *atom, Eigen::Vector3d pos,
                                               QUndoCommand *parent)
                      : QUndoCommand(parent), m_molecule(0)
{
  // Store the original molecule before any modifications are made
  setText(QObject::tr("Bond Centric Manipulation"));
  m_moleculeCopy = *molecule;
  m_molecule = molecule;
  m_atomIndex = atom->GetIdx();
  m_pos = pos;
  undone = false;
}

// ##########  redo  ##########

void WiiMoteMoveCommand::redo()
{
  // Move the specified atom to the location given
  if (undone)
  {
    Molecule newMolecule = *m_molecule;
    *m_molecule = m_moleculeCopy;
    m_moleculeCopy = newMolecule;
  }
  else if (m_atomIndex)
  {
    m_molecule->BeginModify();
    Atom *atom = static_cast<Atom *>(m_molecule->GetAtom(m_atomIndex));
    atom->setPos(m_pos);
    m_molecule->EndModify();
    atom->update();
  }

  QUndoCommand::redo();
}

// ##########  undo  ##########

void WiiMoteMoveCommand::undo()
{
  // Restore our original molecule
  Molecule newMolecule = *m_molecule;
  *m_molecule = m_moleculeCopy;
  m_moleculeCopy = newMolecule;
  undone = true;
}

// ##########  mergeWith  ##########

bool WiiMoteMoveCommand::mergeWith (const QUndoCommand *)
{
  return false;
}

// ##########  id  ##########

int WiiMoteMoveCommand::id() const
{
  //changed from 26011980[manipulatetool]
  return 26011981;
}

#include "wiimotetool.moc"

Q_EXPORT_PLUGIN2(wiimotetool, WiiMoteToolFactory)
