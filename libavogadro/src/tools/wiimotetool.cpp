/**********************************************************************
  WiiMoteTool - Manipulation Tool using WiiMote for Avogadro

  Copyright (C) 2007 by Shahzad Ali

  Parts in this file have been used from cwiid project.

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

#include "drawcommand.h"
#include <avogadro/undosequence.h>

#include <avogadro/primitive.h>
#include <avogadro/color.h>
#include <avogadro/glwidget.h>
#include <avogadro/camera.h>
#include <avogadro/toolgroup.h>

#include <openbabel/obiter.h>
#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

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
//struct m_accCalibration wm_cal, nc_cal;

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
                                                    m_editMode(true),
                                                    m_beginAtomAdded(false),
                                                    m_beginAtom(0),
                                                    m_endAtom(0),
                                                    m_element(6),
                                                    m_bond(0),
                                                    m_prevBond(0),
                                                    m_prevAtomElement(0),
                                                    m_block(false),
                                                    m_glwidget(NULL)
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
  m_forceField = OBForceField::FindForceField( "Ghemical" );
  //action->setShortcut(Qt::Key_F9);
}

// ##########  Desctructor  ##########

WiiMoteTool::~WiiMoteTool()
{
  delete m_referencePoint;
  m_referencePoint = NULL;
  delete m_currentReference;
  m_currentReference = NULL;

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
  m_glwidget = widget;
  m_undo = 0;

  m_lastDraggingPosition = event->pos();
  m_movedSinceButtonPressed = false;
  _buttons = event->buttons();
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
////////////
  if (m_leftButtonPressed & m_rightButtonPressed)
    m_lastRollMovePosition = m_roll;
/////////////

  if (m_editMode)
  {
    m_initialDragginggPosition = event->pos();
  //! List of hits from a selection/pick
    m_hits = widget->hits(event->pos().x()-SEL_BOX_HALF_SIZE,
                          event->pos().y()-SEL_BOX_HALF_SIZE,
                                     SEL_BOX_SIZE,
                                     SEL_BOX_SIZE);

    if(m_leftButtonPressed)
    {
      if (m_hits.size()>0)
      {
        if (m_hits[0].type() == Primitive::AtomType)
        {
      // "alchemy" -- change this atom to a new element
      // Make sure we call BeginModify / EndModify (e.g., PR#1720879)
        widget->molecule()->BeginModify();
        m_beginAtom = (Atom *)widget->molecule()->GetAtom(m_hits[0].name());
        m_prevAtomElement = m_beginAtom->GetAtomicNum();
        m_beginAtom->SetAtomicNum(m_element);
        widget->molecule()->EndModify();
        m_beginAtom->update(); // Make sure to call for a repaint(#1741653).
      // FIXME: This should really be something we can undo
        }
      }
      else
      {
        m_beginAtom = newAtom(widget, event->pos());
        m_beginAtomAdded = true;
        widget->updateGeometry();
        m_beginAtom->update();
      }
    }
    return 0;
  }

/////////////
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
      rollCalibrationVal = 0;

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

QUndoCommand* WiiMoteTool::mouseRelease(GLWidget *widget, const QMouseEvent* event)
{

  if (m_editMode)
  {
    m_undo = editModeMouseRelease(widget, event);
  }

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

QUndoCommand* WiiMoteTool::editModeMouseRelease(GLWidget *widget, const QMouseEvent* event)
{
  QUndoCommand *undo = 0;

  if(_buttons & Qt::LeftButton)
  {
    // we can have a beginAtom w/out bond or endAtom
    // we can hava bond w/out endAtom
    // we cannot have endAtom w/out bond
    // i go through a lot of testing to make the text look prettier and save memory.
    if(m_beginAtomAdded || m_bond)
    {
      AddAtomDrawCommand *beginAtomDrawCommand = 0;
      if(m_beginAtomAdded) {
        beginAtomDrawCommand = new AddAtomDrawCommand(widget->molecule(), m_beginAtom);
        beginAtomDrawCommand->setText(tr("Draw Atom"));
      }

      AddAtomDrawCommand *endAtomDrawCommand = 0;
      if(m_endAtom) {
        endAtomDrawCommand = new AddAtomDrawCommand(widget->molecule(), m_endAtom);
        endAtomDrawCommand->setText(tr("Draw Atom"));
      }

      AddBondDrawCommand *bondCommand = 0;
      if(m_bond) {
        bondCommand = new AddBondDrawCommand(widget->molecule(), m_bond);
        bondCommand->setText(tr("Draw Bond"));
      }

      if(endAtomDrawCommand || (bondCommand && beginAtomDrawCommand))
      {
        UndoSequence *seq = new UndoSequence();
        seq->setText(tr("Draw"));

        if(beginAtomDrawCommand) {
          seq->append(beginAtomDrawCommand);
        }
        if(endAtomDrawCommand) {
          seq->append(endAtomDrawCommand);
        }
        seq->append(bondCommand);
        undo = seq;
      }
      else if(bondCommand)
      {
        undo = bondCommand;
      }
      else
      {
        undo = beginAtomDrawCommand;
      }
    }

    m_beginAtom=0;
    m_bond=0;
    m_endAtom=0;
    m_prevBond=0;
    m_prevBondOrder=0;
    m_prevAtomElement=0;
    m_beginAtomAdded=false;

    // create the undo action for creating endAtom and bond
    //  pass along atom idx, element, vector, bond idx, order, start/end
  }
#ifdef Q_WS_MAC
  // On the Mac, either use a three-button mouse
  // or hold down the Command key (ControlModifier in Qt notation)
  else if( (_buttons & Qt::RightButton) ||
            ((_buttons & Qt::LeftButton) && (event->modifiers() == Qt::ControlModifier)) )
#else
  // Every other platform, use a three-button mouse
    else if(_buttons & Qt::RightButton)
#endif
  {
    m_hits = widget->hits(event->pos().x()-SEL_BOX_HALF_SIZE,
                          event->pos().y()-SEL_BOX_HALF_SIZE,
                                     SEL_BOX_SIZE,
                                     SEL_BOX_SIZE);
    if(m_hits.size())
    {
      qDebug() << m_hits[0].name() << " -- " << m_hits[0].type();

      qDebug() << "bondtype -- " << Primitive::BondType;
      // get our top hit
      if(m_hits[0].type() == Primitive::AtomType)
      {
        undo = new DeleteAtomDrawCommand(widget->molecule(), m_hits[0].name());
//         molecule->DeleteAtom(atom);
//         widget->updateGeometry();
//         molecule->update();
      }
      if(m_hits[0].type() == Primitive::BondType)
      {
        qDebug() << m_hits[0].name();
        undo = new DeleteBondDrawCommand(widget->molecule(), m_hits[0].name());
//         molecule->DeleteAtom(atom);
//         widget->updateGeometry();
//         molecule->update();
      }
    }
  }

  return undo;
}

void WiiMoteTool::translate(GLWidget *widget, const Eigen::Vector3d &what, const QPoint &from, const QPoint &to) const
{
  // Translate the selected atoms in the x and y sense of the view
  //WiiMoteMoveAtomCommand *cmd  = 0;
  //to.setY(10);
  Vector3d fromPos = widget->camera()->unProject(from, what);
  Vector3d toPos = widget->camera()->unProject(to, what);

  MatrixP3d atomTranslation;
  atomTranslation.loadTranslation(toPos - fromPos);

  if (m_clickedAtom)
  {
    widget->molecule()->BeginModify();
    m_clickedAtom->setPos(atomTranslation * m_clickedAtom->pos());
    widget->molecule()->EndModify();
    m_clickedAtom->update();
  }
  //delete newTo;
}

// ##########  mouseMove  ##########

QUndoCommand* WiiMoteTool::mouseMove(GLWidget *widget, const QMouseEvent *event)
{
  m_glwidget = widget;
  if (!m_glwidget->molecule()) {
    return 0;
  }

  if (m_editMode)
    return editModeMouseMove(m_glwidget, event);
  //m_undo = new WiiMoteMoveAtomCommand(widget->molecule());

  // Get the currently selected atoms from the view
  QList<Primitive *> currentSelection = m_glwidget->selectedPrimitives();

  QPoint deltaDragging = event->pos() - m_lastDraggingPosition;

  // Manipulation can be performed in two ways - centred on an individual atom
  if (m_clickedAtom)
  {
      if (m_leftButtonPressed && !m_rightButtonPressed)
      {
      // translate the molecule following mouse movement
      QPoint newTo;
      QPoint to = event->pos();
      newTo.setX(to.x());
      newTo.setY(m_lastDraggingPosition.y());
      translate(m_glwidget, m_clickedAtom->pos(), m_lastDraggingPosition, newTo);
      }
      else if (event->buttons() & Qt::MidButton)
      {
      /**
      if (deltaDragging.y() == 0)
        // Perform the rotation
        tilt(widget, m_clickedAtom->pos(), deltaDragging.x());
      else
        // Perform the zoom toward clicked atom
        zoom(widget, m_clickedAtom->pos(), deltaDragging.y());
      **/
      }
      else if (!m_leftButtonPressed && m_rightButtonPressed )
      {
      // translate the molecule following mouse movement
        QPoint newTo;
        QPoint to = event->pos();
        newTo.setX(m_lastDraggingPosition.x());
        newTo.setY(to.y());
        translate(m_glwidget, m_clickedAtom->pos(), m_lastDraggingPosition, newTo);
      }
  }
  else if (!m_clickedAtom)
  {
    // rotation around the center of the molecule
    Navigate::rotate(m_glwidget, m_glwidget->center(), deltaDragging.x(), deltaDragging.y());
  }

  m_lastDraggingPosition = event->pos();
  widget->update();

  return 0;
}

QUndoCommand* WiiMoteTool::editModeMouseMove(GLWidget *widget, const QMouseEvent *event)
{
  Molecule *molecule = widget->molecule();
  if(!molecule) {
    return 0;
  }

  if((_buttons & Qt::LeftButton) && m_beginAtom)
  {
    m_hits = widget->hits(event->pos().x()-SEL_BOX_HALF_SIZE,
                          event->pos().y()-SEL_BOX_HALF_SIZE,
                                     SEL_BOX_SIZE,
                                     SEL_BOX_SIZE);

    bool hitBeginAtom = false;
    Atom *existingAtom = 0;
    if(m_hits.size())
    {
      // parse our hits.  we want to know
      // if we hit another existingAtom that is not
      // the m_endAtom which we created
      for(int i=0; i < m_hits.size() && !hitBeginAtom; i++)
      {
        if(m_hits[i].type() == Primitive::AtomType)
        {
          // hit the same atom either moved here from somewhere else
          // or were already here.
          if(m_hits[i].name() == m_beginAtom->GetIdx())
          {
            hitBeginAtom = true;
          }
          else if(!m_endAtom)
          {
            existingAtom = (Atom *)molecule->GetAtom(m_hits[i].name());
          }
          else
          {
            if(m_hits[i].name() != m_endAtom->GetIdx())
            {
              existingAtom = (Atom *)molecule->GetAtom(m_hits[i].name());
            }
          }
        }
      }
    }
    if(hitBeginAtom)
    {
      if(m_endAtom)
      {
        molecule->DeleteAtom(m_endAtom);
        // widget->updateGeometry();
        m_bond = 0;
        m_endAtom = 0;
        m_prevAtomElement = m_beginAtom->GetAtomicNum();
        m_beginAtom->SetAtomicNum(m_element);
        // m_beginAtom->update();
      }
      else if(m_bond)
      {
        Atom *oldAtom = (Atom *)m_bond->GetEndAtom();
        oldAtom->DeleteBond(m_bond);
        molecule->DeleteBond(m_bond);
        m_bond=0;
        m_prevAtomElement = m_beginAtom->GetAtomicNum();
        m_beginAtom->SetAtomicNum(m_element);
        // m_beginAtom->update();
      }
    }
    else
    {
      if(m_prevAtomElement)
      {
        m_beginAtom->SetAtomicNum(m_prevAtomElement);
        m_prevAtomElement = 0;
      }

      // we hit an existing atom != m_endAtom

      if(existingAtom)
      {
        Bond *existingBond = (Bond *)molecule->GetBond(m_beginAtom, existingAtom);
        cout << "13" << endl;
        if(!existingBond) {
          if(m_prevBond)
          {
            cout << "14" << endl;
            m_prevBond->SetBondOrder(m_prevBondOrder);
            // m_prevBond->update();
            m_prevBond = 0;
            m_prevBondOrder = 0;
          }
          cout << "15" << endl;
          if(m_bond) {
            if(m_endAtom) {
              m_endAtom->DeleteBond(m_bond);
              molecule->DeleteAtom(m_endAtom);
              m_endAtom = 0;
            } else {
              Atom *oldAtom = (Atom *)m_bond->GetEndAtom();
              oldAtom->DeleteBond(m_bond);
            }
            m_bond->SetEnd(existingAtom);
            existingAtom->AddBond(m_bond);
            // m_bond->update();
          } else {
            m_bond = newBond(molecule, m_beginAtom, existingAtom);
            // m_bond->update();
          }
        }
        // (existingBond)
        else {
          if(m_prevBond && m_prevBond != existingBond) {
            m_prevBond->SetBondOrder(m_prevBondOrder);
            // m_prevBond->update();
            m_prevBond = 0;
            m_prevBondOrder = 0;
          }
          if(!m_prevBond) {
            m_prevBond = existingBond;
            m_prevBondOrder = existingBond->GetBO();
            existingBond->SetBondOrder(m_bondOrder);
            // existingBond->update();
          }

          if(m_bond && m_bond != existingBond) {
            if(m_endAtom) {
              // will delete bonds too (namely m_bond)
              molecule->DeleteAtom(m_endAtom);
              m_endAtom = 0;
            } else {
              molecule->DeleteBond(m_bond);
            }
            m_bond = 0;
          }
        }
      }
      // (!existingAtom && !hitBeginAtom)
      else if(!m_endAtom)
      {
        if(m_prevBond) {
          m_prevBond->SetBondOrder(m_prevBondOrder);
          // m_prevBond->update();
          m_prevBond = 0;
          m_prevBondOrder = 0;
        }
        m_endAtom = newAtom(widget, event->pos());
        if(!m_bond)
        {
          m_bond = newBond(molecule, m_beginAtom, m_endAtom);
        }
        else
        {
          Atom *oldAtom = (Atom *)m_bond->GetEndAtom();
          oldAtom->DeleteBond(m_bond);
          m_bond->SetEnd(m_endAtom);
          m_endAtom->AddBond(m_bond);
        }
        // m_bond->update();
        // m_endAtom->update();
        widget->updateGeometry();
      }
      else
      {
        moveAtom(widget, m_endAtom, event->pos());
        // widget->updateGeometry();
        // m_endAtom->update();
      }
    }
    molecule->update();
  }

  return 0;
}


Atom *WiiMoteTool::newAtom(GLWidget *widget, const QPoint& p)
{
  // GRH (for reasons I don't understand, calling Begin/EndModify here
  // causes crashes with multiple bond orders
  // (need to investigate, probable OB bug.

  widget->molecule()->BeginModify();
  Atom *atom = static_cast<Atom*>(widget->molecule()->NewAtom());
  moveAtom(widget, atom, p);
  atom->SetAtomicNum(element());
  widget->molecule()->EndModify();

  return atom;
}

void WiiMoteTool::setBondOrder( int index )
{
  m_bondOrder = index;
}

int WiiMoteTool::bondOrder() const
{
  return m_bondOrder;
}

void WiiMoteTool::setElement( int index )
{
  m_element = index;
}

int WiiMoteTool::element() const
{
  return m_element;
}

void WiiMoteTool::moveAtom(GLWidget *widget, Atom *atom, const QPoint& p)
{
  Eigen::Vector3d refPoint;
  if(m_beginAtom) {
    refPoint = m_beginAtom->pos();
  } else {
    refPoint = widget->center();
  }
  Eigen::Vector3d newAtomPos = widget->camera()->unProject(p, refPoint);

  atom->setPos(newAtomPos);
}


Bond *WiiMoteTool::newBond(Molecule *molecule, Atom *beginAtom, Atom *endAtom)
{
  molecule->BeginModify();
  Bond *bond = (Bond *)molecule->NewBond();
  bond->SetBondOrder(bondOrder());
  bond->SetBegin(beginAtom);
  bond->SetEnd(endAtom);
  beginAtom->AddBond(bond);
  endAtom->AddBond(bond);
  molecule->EndModify();

  return bond;
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
int m_count = 0;
QUndoCommand* WiiMoteTool::wiimoteRoll(__s32 roll, __s32 prevRoll)
{
  //cout << "hellow there " << roll << " " << prevRoll << endl;
  // Qt::LeftButton  == 1
  if (!m_glwidget)
    return 0;

  if (!m_glwidget->molecule()) {
    return 0;
  }

  m_roll = roll;
  double relRoll = 0;

  if (roll > 500)
    relRoll = 500;
  else if (roll < -500)
    relRoll = -500;
  else relRoll = roll;

  relRoll = relRoll / 500;

  // Mouse navigation has two modes - atom centred when an atom is clicked
  // and scene if no atom has been clicked.
  if (!m_editMode)
  if (m_leftButtonPressed)
  {
    if (m_clickedBond && m_selectedBond && m_referencePoint)
    {
      /*
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
      }*/
    }
    else if (m_clickedAtom && m_rightButtonPressed)
    {
      if (m_count < 5)
      {
        m_count++;
        return 0;
      }
      m_count = 0;
      //Do atom rotation.
      Vector3d pos = m_clickedAtom->pos();
      __s32 smallChange;
      //smallChange = deltaRoll / 50;
      smallChange = relRoll;
      Vector3d newPos = Vector3d(pos.x(), pos.y(), pos.z() + smallChange);
      m_clickedAtom->setPos(newPos);
    }
    else {
      // rotation around the center of the molecule
      //Navigate::rotate(m_glwidget, m_glwidget->center(), deltaDragging.x(),
        //               deltaDragging.y());
    }
  }
#ifdef Q_WS_MAC
  // On the Mac, either use a three-button mouse
  // or hold down the Option key (AltModifier in Qt notation)
//  else if ((event->buttons() & Qt::MidButton) || (event->buttons() &
  //          Qt::LeftButton && event->modifiers() & Qt::AltModifier))
#else
  else if (m_midButtonPressed)
#endif
  {/**
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
    }**/
  }
  else if (m_rightButtonPressed)
  {/*
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
    }**/
  }

  
  m_glwidget->update();

  return 0;
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

QUndoCommand* WiiMoteTool::wiimoteRumble(int rumble)
{
  if (rumble < 30)
    return 0;

  if (!m_glwidget)
    return 0;
  if (!m_glwidget->molecule())
    return 0;

  if (m_block)
    return 0;
  else
    m_block = true;

  if ( !m_forceField->Setup( *m_glwidget->molecule() ) ) {
    qWarning() << "GhemicalCommand: Could not set up force field on " << m_glwidget->molecule();
    m_block = false;
    return 0;
  }

  m_forceField->SteepestDescent(1,pow(10.0, -7 ), OBFF_NUMERICAL_GRADIENT);
  m_forceField->UpdateCoordinates( *m_glwidget->molecule() );
  m_glwidget->molecule()->update();
  m_glwidget->update();
  m_block = false;
  return 0;
}


void WiiMoteTool::cwiid_callback(cwiid_wiimote_t *wiimote, int mesg_count,
                    union cwiid_mesg mesg[], struct timespec *timestamp)
{
  int i;

  for (i=0; i < mesg_count; i++) {
    switch (mesg[i].type) {
      case CWIID_MESG_BTN:
        sendButtonEvent((struct cwiid_btn_mesg *) &mesg[i]);
        //cout << "Msg Btn" << endl;
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
  processIR(mesg_count, mesg);
  processAcc(mesg_count, mesg);
  sendEvent(m_conf, EV_SYN, SYN_REPORT, 0);
}

void WiiMoteTool::sendButtonEvent(struct cwiid_btn_mesg *mesg)
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
        //HACK ////////////////////
        if (pressed & m_conf.wiimote_bmap[CONF_WM_BTN_MINUS].mask)
        {
          m_editMode = true;
          cout << "Edit Mode" << endl;
        } else
          if (pressed & m_conf.wiimote_bmap[CONF_WM_BTN_PLUS].mask)
        {
          m_editMode = false;
          cout << "Move Mode" << endl;
        }
        ///////////////////////////
        sendEvent(m_conf, EV_KEY, m_conf.wiimote_bmap[i].action, 1);
      }
      else if (released & m_conf.wiimote_bmap[i].mask) {
        sendEvent(m_conf, EV_KEY, m_conf.wiimote_bmap[i].action, 0);
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
    sendEvent(&conf, conf.amap[CONF_WM_AXIS_DPAD_X].axis_type,
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
    sendEvent(&conf, conf.amap[CONF_WM_AXIS_DPAD_Y].axis_type,
                conf.amap[CONF_WM_AXIS_DPAD_Y].action, axis_value);
  }**/
}


int WiiMoteTool::sendEvent(struct WiiConfig conf, __u16 type, __u16 code, __s32 value)
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
    cout << "Error on sendEvent" << endl;
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
    //Cancelled
    return;
  }
  else
  {
    ptO2Object = this;
    char reset_bdaddr = 0;

    if (bacmp(&bdaddr, BDADDR_ANY) == 0) {
      reset_bdaddr = 1;
    }

    if ((wiimote = cwiid_open(&bdaddr, CWIID_FLAG_MESG_IFC)) == NULL) {
      displayError("Unable to connect to Wii-Remote.\nMake sure it is in discoverable mode.");
    } else if (cwiid_set_mesg_callback(wiimote, &Avogadro::WiiMoteTool::cwiid_callback_wrapper)) {
      displayError("Error setting callback.");

      if (cwiid_close(wiimote)) {
        displayError("Error on disconnect.");
    }
    wiimote = NULL;
    }
    else
    {
      int request;
      __u16 m_axisType = 3;
      __u16 actionX = 0;
      __u16 actionY = 1;
      __u16 actionZ = 2;

      //connected
      ui.m_buttonConnect->setText("Connected");
      ui.m_buttonConnect->setEnabled(false);
      ui.m_buttonDisconnect->setEnabled(true);

      if (cwiid_get_acc_cal(wiimote, CWIID_EXT_NONE, &m_accCalibration)) {
        displayError("Unable to retrieve accelerometer.");
      }

      loadWiiConfig();

      ///////////////////////////////////////////////////////////////
      /* Setup UInput */
      ///////////////////////////////////////////////////////////////
      char *uinput_filename[] = {"/dev/uinput", "/dev/input/uinput",
        "/dev/misc/uinput"};
        int uinputfilenamecount = 3;

        int i;
        int j;

        /* Open uinput device */
        for (i=0; i < uinputfilenamecount; i++) {
          m_conf.fd = open(uinput_filename[i], O_RDWR);
          if (m_conf.fd >= 0)
          {
            break;
          }
        }

        if (m_conf.fd < 0) {
          displayError("Unable to open uinput");
          return;
        }

        if (write(m_conf.fd, &m_conf.dev, sizeof m_conf.dev) != sizeof m_conf.dev) {
          displayError("error on uinput device setup");
          close(m_conf.fd);
          return;
        }

        if (m_conf.ff) {
          if (ioctl(m_conf.fd, UI_SET_EVBIT, EV_FF) < 0) {
            displayError("error on uinput ioctl");
            close(m_conf.fd);
            return;
          }
          if (ioctl(m_conf.fd, UI_SET_FFBIT, FF_RUMBLE) < 0) {
            displayError("error on uinput ioctl" );
            close(m_conf.fd);
            return;
          }
        }

        if (ioctl(m_conf.fd, UI_SET_EVBIT, EV_KEY) < 0) {
          displayError("error on uinput ioctl");
          close(m_conf.fd);
          return;
        }

      for (i=0; i < CONF_WM_BTN_COUNT; i++) {
        if (m_conf.wiimote_bmap[i].active) {
          if (ioctl(m_conf.fd, UI_SET_KEYBIT, m_conf.wiimote_bmap[i].action) < 0) {
            displayError("error on uinput ioctl");
            close(m_conf.fd);
            return;
          }
        }
      }

      if (ioctl(m_conf.fd, UI_SET_EVBIT, m_axisType) < 0) {
        displayError("Error uinput ioctl");
      }

      request = (m_axisType == EV_ABS) ? UI_SET_ABSBIT : UI_SET_RELBIT;
      if (ioctl(m_conf.fd, request, actionX) < 0) {
        displayError("Error uinput ioctl");
      }

      if (ioctl(m_conf.fd, request, actionY) < 0) {
        displayError("Error uinput ioctl");
      }

      if ((m_axisType == EV_ABS) && ((actionX == ABS_X) || (actionY == ABS_Y) || (actionZ == ABS_Z))) {
        if (ioctl(m_conf.fd, UI_SET_EVBIT, EV_REL) < 0) {
          displayError("error uinput ioctl");
          close(m_conf.fd);
          return;
        }

        if (ioctl(m_conf.fd, UI_SET_RELBIT, actionX) < 0) {
          displayError("error uinput ioctl");
          close(m_conf.fd);
          return;
        }

        if (ioctl(m_conf.fd, UI_SET_RELBIT, actionY) < 0) {
          displayError("error uinput ioctl");
          close(m_conf.fd);
          return;
        }
      }

      if (ioctl(m_conf.fd, UI_DEV_CREATE) < 0) {
        displayError("Error on uinput dev create");
        close(m_conf.fd);
      }

      setReportMode();
      cwiid_request_status(wiimote);

   if (reset_bdaddr) { bdaddr = *BDADDR_ANY; }

  }
  }
}

void WiiMoteTool::processAcc(int mesg_count, union cwiid_mesg mesg[])
{
  int i;
  struct ProcessedWiiMoteData *data = NULL;

  __s32 oldAcc = m_accData.axes[0].value;

  for (i=0; i < mesg_count; i++) {
    switch (mesg[i].type) {
      case CWIID_MESG_ACC:
        data = processAccData(&mesg[i].acc_mesg);
        if (data->axes[0].value != oldAcc)
        {
          //data->axes[0].value = data->axes[0].value | 0x40000000;
         // cout <<"_" << data->axes[0].value << endl;
         // sendEvent(m_conf, EV_ABS, ABS_X, data->axes[0].value);
          sendCustomEvent(data->axes[0].value, oldAcc);
        }
        if (data->axes[2].value > 25)
          sendRumbleEvent(data->axes[2].value);
        break;
      default:
        break;
    }
  }
}

void WiiMoteTool::sendCustomEvent(__s32 rollVal, __s32 prevRollVal)
{
  wiimoteRoll(rollVal, prevRollVal);
}

void WiiMoteTool::sendRumbleEvent(int rumble)
{
  wiimoteRumble(rumble);
}

struct ProcessedWiiMoteData *WiiMoteTool::processAccData(struct cwiid_acc_mesg *mesg)
{
  double a, a_x, a_y, a_z;
  double roll, pitch;
  double newAmount = 0.1;
  double oldAmount = 1 - newAmount;
  float Roll_Scale = 1.0;
  float Pitch_Scale = 1.0;
  float X_Scale = 1.0;
  float Y_Scale = 1.0;

  a_x = (((double)mesg->acc[CWIID_X] - m_accCalibration.zero[CWIID_X]) /
      (m_accCalibration.one[CWIID_X] - m_accCalibration.zero[CWIID_X]))*newAmount +
      a_x*oldAmount;
  a_y = (((double)mesg->acc[CWIID_Y] - m_accCalibration.zero[CWIID_Y]) /
      (m_accCalibration.one[CWIID_Y] - m_accCalibration.zero[CWIID_Y]))*newAmount +
      a_y*oldAmount;
  a_z = (((double)mesg->acc[CWIID_Z] - m_accCalibration.zero[CWIID_Z]) /
      (m_accCalibration.one[CWIID_Z] - m_accCalibration.zero[CWIID_Z]))*newAmount +
      a_z*oldAmount;

  a = sqrt(pow(a_x,2)+pow(a_y,2)+pow(a_z,2));
  roll = atan(a_x/a_z);
  if (a_z <= 0.0) {
    roll += PI * ((a_x > 0.0) ? 1 : -1);
  }

  pitch = atan(a_y/a_z*cos(roll));

  m_accData.axes[0].value = roll  * 1000 * Roll_Scale;
  m_accData.axes[1].value = pitch * 1000 * Pitch_Scale;
  m_accData.axes[2].value = a * 100;

/**  if ((a > 0.85) && (a < 1.15)) {
    if ((fabs(roll)*(180/PI) > 10) && (fabs(pitch)*(180/PI) < 80)) {
      m_accData.axes[2].valid = 1;
      m_accData.axes[2].value = roll * 5 * X_Scale;
    }
    else {
      m_accData.axes[2].valid = 0;
    }
    if (fabs(pitch)*(180/PI) > 10) {
      m_accData.axes[3].valid = 1;
      m_accData.axes[3].value = pitch * 10 * Y_Scale;
    }
    else {
      m_accData.axes[3].valid = 0;
    }
  }
  else {
    m_accData.axes[2].valid = 0;
    m_accData.axes[3].valid = 0;
  }**/
  return &m_accData;
}

void WiiMoteTool::processIR(int mesg_count, union cwiid_mesg mesg[])
{
  static union cwiid_mesg irMessage[CWIID_MAX_MESG_COUNT];
  int irMessageCount = 0;
  int i;
  uint8_t flag;
  ProcessedWiiMoteData *data = NULL;

  for (i=0; i < mesg_count; i++) {
    switch (mesg[i].type) {
      case CWIID_MESG_STATUS:
        flag = CWIID_RPT_STATUS;
        break;
      case CWIID_MESG_BTN:
        flag = CWIID_RPT_BTN;
        break;
      case CWIID_MESG_ACC:
        flag = CWIID_RPT_ACC;
        break;
      case CWIID_MESG_IR:
        flag = CWIID_RPT_IR;
        break;
      case CWIID_MESG_NUNCHUK:
        flag = CWIID_RPT_NUNCHUK;
        break;
      case CWIID_MESG_CLASSIC:
        flag = CWIID_RPT_CLASSIC;
        break;
      default:
        break;
    }
    if (m_reportMode & flag) {
      /* TODO: copy correct (smaller) message size */
      memcpy(&irMessage[irMessageCount++], &mesg[i], sizeof mesg[i]);
    }
  }

  if (irMessageCount > 0) {
    data = processIRData(mesg_count, mesg); 
    if (!data) {
      return;
    }

    sendEvent(m_conf, EV_ABS, ABS_X, data->axes[0].value);
    sendEvent(m_conf, EV_ABS, ABS_Y, data->axes[1].value);
    //sendEvent(m_conf, EV_ABS, ABS_Z, m_data.axes[2].value);
  }
}

struct ProcessedWiiMoteData *WiiMoteTool::processIRData(int mesg_count, union cwiid_mesg mesg[])
{
  double newAmount = 0.1;
  double oldAmount = 1 - newAmount;

  int src_index = -1;
  int debounce = 0;
  uint8_t old_flag;

  int i;
  uint8_t flag;
  struct cwiid_ir_mesg *ir_mesg;

  ir_mesg = NULL;

  for (i=0; i < mesg_count; i++) {
    if (mesg[i].type == CWIID_MESG_IR) {
      ir_mesg = &mesg[i].ir_mesg;
    }
  }

  if (!ir_mesg) {
    return NULL;
  }

  /* invalidate src index if source is no longer present */
  if ((src_index != -1) && !ir_mesg->src[src_index].valid) {
    if (debounce > DEBOUNCE_THRESHOLD) {
      src_index = -1;
    }
    else {
      debounce++;
    }
  }
  else {
    debounce = 0;
  }

  /* of not set, pick largest available source */
  if (src_index == -1) {
    for (i=0; i < CWIID_IR_SRC_COUNT; i++) {
      if (ir_mesg->src[i].valid) {
        if ((src_index == -1) ||
             (ir_mesg->src[i].size > ir_mesg->src[src_index].size)) {
          src_index = i;
             }
      }
    }
  }

  /* LEDs */
  switch (src_index) {
    case 0:
      flag = CWIID_LED1_ON;
      break;
    case 1:
      flag = CWIID_LED2_ON;
      break;
    case 2:
      flag = CWIID_LED3_ON;
      break;
    case 3:
      flag = CWIID_LED4_ON;
      break;
    default:
      flag = 0;
      break;
  }
  if (flag != old_flag) {
    cwiid_set_led(wiimote, flag);
    old_flag = flag;
  }

  if ((src_index == -1) || !ir_mesg->src[src_index].valid) {
    m_irData.axes[0].valid = m_irData.axes[1].valid = m_irData.axes[2].valid = 0;
  }
  else {
    m_irData.axes[0].valid = m_irData.axes[1].valid = m_irData.axes[2].valid = 1;
    m_irData.axes[0].value = newAmount * (CWIID_IR_X_MAX -
        ir_mesg->src[src_index].pos[CWIID_X])
        + oldAmount * m_irData.axes[0].value;
    m_irData.axes[1].value = newAmount * ir_mesg->src[src_index].pos[CWIID_Y]
        + oldAmount * m_irData.axes[1].value;


    /*
    int l1 = 0;
    int l2 = 0;
    //get the two largest source size
    for (i=0; i < CWIID_IR_SRC_COUNT; i++) {
      if (ir_mesg->src[i].valid) {
        if (ir_mesg->src[i].size > ir_mesg->src[l1].size) {
        l2 = l1;
        l1 = i;
        }
      }
    }
    //the increase and decrease of distance between the two
    //sources determines the Z axes position.
    //Hence, there must be atleast two ir source to work with
    ////////// This technique does not work /////////////////
    if (l1 != l2) {
      m_data.axes[2].value = abs((( newAmount * ir_mesg->src[l1].pos[CWIID_X])/CWIID_IR_X_MAX *
          ( newAmount * ir_mesg->src[l1].pos[CWIID_Y]) / CWIID_IR_Y_MAX) - 
          ((newAmount * ir_mesg->src[l2].pos[CWIID_X]) / CWIID_IR_X_MAX *
          ( newAmount * ir_mesg->src[l2].pos[CWIID_Y]) / CWIID_IR_Y_MAX));
    }
    cout << ir_mesg->src[src_index].pos[CWIID_X]
    << " " << ir_mesg->src[src_index].pos[CWIID_Y]
        << " " << m_data.axes[2].value << endl;
    **/

    if (m_irData.axes[0].value > CWIID_IR_X_MAX - X_EDGE) {
      m_irData.axes[0].value = CWIID_IR_X_MAX - X_EDGE;
    }
    else if (m_irData.axes[0].value < X_EDGE) {
      m_irData.axes[0].value = X_EDGE;
    }
    if (m_irData.axes[1].value > CWIID_IR_Y_MAX - Y_EDGE) {
      m_irData.axes[1].value = CWIID_IR_Y_MAX - Y_EDGE;
    }
    else if (m_irData.axes[1].value < Y_EDGE) {
      m_irData.axes[1].value = Y_EDGE;
    }
  }
  return &m_irData;
}

void WiiMoteTool::displayError(QString message)
{
  QMessageBox::information( NULL, "Wii-Remote", message);
}

void WiiMoteTool::loadWiiConfig()
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
  
  m_conf.dev.absmax[0] = CWIID_IR_X_MAX - X_EDGE;
  m_conf.dev.absmax[1] = CWIID_IR_Y_MAX - Y_EDGE;
  m_conf.dev.absmin[0] = X_EDGE;
  m_conf.dev.absmin[1] = Y_EDGE;
  m_conf.dev.absfuzz[0] = 0;
  m_conf.dev.absfuzz[1] = 0;
  m_conf.dev.absflat[0] = 0;
  m_conf.dev.absflat[1] = 0;

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

  m_conf.wiimote_bmap[CONF_WM_BTN_MINUS].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_MINUS].action = KEY_1;

  m_conf.wiimote_bmap[CONF_WM_BTN_PLUS].active = 1;
  m_conf.wiimote_bmap[CONF_WM_BTN_PLUS].action = KEY_2;

}

void WiiMoteTool::setReportMode()
{ 
  m_reportMode = CWIID_RPT_STATUS | CWIID_RPT_BTN;
  m_reportMode |= CWIID_RPT_IR;
  m_reportMode |= CWIID_RPT_ACC;
  /*
  if (gtk_check_menu_item_get_active(GTK_CHECK_MENU_ITEM(chkExt))) {
    rpt_mode |= CWIID_RPT_EXT;
  }**/
  if (cwiid_set_rpt_mode(wiimote, m_reportMode)) {
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
  ui.m_buttonConnect->setText("Connect");
  ui.m_buttonConnect->setEnabled(true);
  ui.m_buttonDisconnect->setEnabled(false);
}

void WiiMoteTool::cleanup()
{
  close(m_conf.fd);
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
