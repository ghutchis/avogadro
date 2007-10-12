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

#ifndef __WIIMOTETOOL_H
#define __WIIMOTETOOL_H

#include "skeletontree.h"
#include "ui_wiimotetool.h"
#include "cwiid.h"
#include <signal.h>

#include <avogadro/glwidget.h>
#include <avogadro/tool.h>

#include <openbabel/mol.h>
#include <openbabel/forcefield.h>

#include <QGLWidget>
#include <QObject>
#include <QStringList>
#include <QImage>
#include <QAction>
#include <QUndoCommand>

#include <QLabel>
#include <QSpinBox>
#include <QCheckBox>
#include <QGridLayout>

#include <stdint.h>
#include <stdio.h>
#include <linux/input.h>
#include <linux/uinput.h>

#define CONF_WM_BTN_COUNT 11
#define WMPLUGIN_MAX_AXIS_COUNT 6

#define DEBOUNCE_THRESHOLD  50
#define X_EDGE  50
#define Y_EDGE  50
#define PI  3.14159265358979323

namespace Avogadro {

  struct ButtonsMap {
    unsigned char active;
    uint16_t mask;
    __u16 action;
  };

  struct WiiConfig {
    int fd;
    struct uinput_user_dev dev;
    unsigned char ff;
    struct ButtonsMap wiimote_bmap[CONF_WM_BTN_COUNT];
  };

  struct AxesData {
    char valid;
    __s32 value;
  };

  struct ProcessedWiiMoteData {
    uint16_t buttons;
    struct AxesData axes[WMPLUGIN_MAX_AXIS_COUNT];
  };

  /**
   * @class WiiMoteTool
   * @brief Manipulation Tool using Wii-Mote
   * @author Shahzad Ali, Ross Braithwaite, James Bunt
   *
   * This class is a molecule manipulation system based on bond-centric
   * design as apposed to points in free space design.  It is based off
   * the NavigationTool class by Marcus D. Hanwell.
   */
  class AddAtomCommand;
  class WiiMoteTool : public Tool
  {
    Q_OBJECT

    public:
      //! Constructor
      WiiMoteTool(QObject *parent = 0);
      //! Deconstructor
      virtual ~WiiMoteTool();

      //! \name Description methods
      //@{
      //! Tool Name (ie Draw)
      virtual QString name() const { return(tr("WiiMote")); }
      //! Tool Description (ie. Draws atoms and bonds)
      virtual QString description() const { 
      return(tr("Wii Mote Manipulation Tool")); }
      //@}

      //! \name Tool Methods
      //@{
      //! \brief Callback methods for ui.actions on the canvas.
      virtual QUndoCommand* mousePress(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* mouseRelease(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* mouseMove(GLWidget *widget, const QMouseEvent *event);
      virtual QUndoCommand* wheel(GLWidget *widget, const QWheelEvent *event);

      //FIXME HACK 
      virtual QUndoCommand* wiimoteRoll(__s32 roll, __s32 prevRoll);
      virtual QUndoCommand* wiimoteRumble(int rumble);
      QUndoCommand* editModeMouseRelease(GLWidget *widget, const QMouseEvent* event);
      QUndoCommand* editModeMouseMove(GLWidget *widget, const QMouseEvent *event);
      //@}

      virtual int usefulness() const;

      virtual bool paint(GLWidget *widget);

      virtual QWidget *settingsWidget();

      void setBondOrder(int i);
      int bondOrder() const;
      void setElement(int i);
      int element() const;
      
    public Q_SLOTS:
      /**
       * Sets the snap-to angle to a given angle in degrees.
       *
       * @param newAngle The new value for the snap-to angle.
       */
      void snapToAngleChanged(int newAngle);

      /**
       * Sets whether or not snap-to is enabled.
       *
       * @param state The state of the check box relating to whether or not
       *              snap-to is enabled.
       *
       *              Qt:Checked - enable snap-to.
       *              Qt:Unchecked - disable snap-to.
       */
      void snapToCheckBoxChanged(int state);

      /**
       * Sets whether or not to show angles.
       *
       * @param state The state of the check box relating to whether or not
       *              to show angles.
       *
       *              Qt:Checked - show angles.
       *              Qt:Unchecked - don't show angles.
       */
      void showAnglesChanged(int state);
      
      /**
       * Function to be called when connect button is clicked.
       *
       *
       */
      void connectClicked();
      /**
       * Function to be called when connect button is clicked.
       *
       *
       */
      void disconnectClicked();


    protected:
      GLWidget *          m_glwidget;
      QWidget *           m_settingsWidget;

      Ui::WiiMoteSettings ui;

      Atom *              m_clickedAtom;
      Bond *              m_clickedBond;
      Bond *              m_selectedBond;

      SkeletonTree *      m_skeleton;

      Eigen::Vector3d *   m_referencePoint;
      Eigen::Vector3d *   m_currentReference;
      bool                m_snapped;
      ToolGroup *         m_toolGroup;

      QUndoCommand *      m_undo; // The current undo command

      bool                m_leftButtonPressed;  // rotation
      bool                m_midButtonPressed;   // scale / zoom
      bool                m_rightButtonPressed; // translation
      bool                m_movedSinceButtonPressed;

      bool                m_showAngles;
      bool                m_snapToEnabled;

      int                 m_snapToAngle; // In degrees

      QPoint              m_lastDraggingPosition;

      QLabel *            m_snapToAngleLabel;
      QLabel *            m_spacer;
      QCheckBox *         m_showAnglesBox;
      QCheckBox *         m_snapToCheckBox;
      QSpinBox *          m_snapToAngleBox;
      QGridLayout *       m_layout;

      Qt::MouseButtons _buttons;
      QPoint              m_initialDragginggPosition;
      bool m_beginAtomAdded;
      Atom *m_beginAtom;
      Atom *m_endAtom;
      int m_element;
      Bond *m_bond;
      int m_bondOrder;
      int m_prevAtomElement;
      Bond *m_prevBond;
      int m_prevBondOrder;
      QList<GLHit> m_hits;

      OpenBabel::OBForceField*  m_forceField;
      bool m_block;

      double m_avgChange;
      int m_valuesAveraged;
      double rollCalibrationVal;
      double m_lastRollMovePosition;
      double m_roll;
      bool m_editMode;

      cwiid_wiimote_t *wiimote;
      WiiConfig m_conf;
      struct input_event m_event;
      struct ProcessedWiiMoteData m_irData;
      struct ProcessedWiiMoteData m_accData;
      struct acc_cal m_accCalibration;

      uint8_t m_reportMode;

      int sendEvent(struct WiiConfig conf, __u16 type, __u16 code, __s32 value);
      void sendButtonEvent(struct cwiid_btn_mesg *mesg);
      void cwiid_callback(cwiid_wiimote_t *wiimote, int mesg_count,
                          union cwiid_mesg mesg[], struct timespec *timestamp);
      static void cwiid_callback_wrapper(cwiid_wiimote_t *wiimote, int mesg_count,
                                         union cwiid_mesg mesg[], struct timespec *timestamp);
      ProcessedWiiMoteData *processIRData(int mesg_count, union cwiid_mesg mesg[]);
      ProcessedWiiMoteData *processAccData(struct cwiid_acc_mesg *mesg);

      void processIR(int mesg_count, union cwiid_mesg mesg[]);
      void processAcc(int mesg_count, union cwiid_mesg mesg[]);
      void setReportMode();
      void loadWiiConfig();
      void cleanup();
      void displayError(QString message);

      //#######Hack######################################################
      void sendCustomEvent(__s32 rollVal, __s32 prevRollVal);
      void sendRumbleEvent(int rumble);
      ///////////////////////////////////////////////////////////////////
      void translate(GLWidget *widget, const Eigen::Vector3d &what, const QPoint &from, const QPoint &to) const;
      
      //ATOM ADD METHODS////////////////////////////////////////////////
      Atom *newAtom(GLWidget *widget, const QPoint& p);
      void moveAtom(GLWidget *widget, Atom *atom, const QPoint& p);
      Bond *newBond(Molecule *molecule, Atom *beginAtom, Atom *endAtom);
      //////////////////////////////////////////////////////////////////

      //! \name Construction Plane/Angles Methods
      //@{
      //! \brief Methods used to construct and draw the angle-sectors, the construction plane, and the rotation-sphere

      /**
       * Checks whether a given atom is at either end of a given bond.
       *
       * @param atom The atom that is being examined for membership of the given bond.
       * @param bond The bond that is being examined to see if the given atom is
       *             attached to it.
       *
       * @return True if the given atom is the begin or end atom of the given
       *         bond, false otherwise, or if either of the pointers point to NULL.
       */
      bool isAtomInBond(Atom *atom, Bond *bond);
      /**
       * Draws a sector that shows the angle between two lines from a given origin.
       *
       * @param widget The widget this angle-sector will be drawn on.
       * @param origin The origin around which this angle is being calculated.
       * @param direction1 A vector that defines the line from the given origin
       *                   through this point.
       * @param direction2 A vector that defines the line from the given origin
       *                   through this second point.
       */
      void drawAngleSector(GLWidget *widget, Eigen::Vector3d origin,
                           Eigen::Vector3d direction1, Eigen::Vector3d direction2);

      /**
       * Draws sectors around a given atom representing the angles between neighbouring
       * atoms bonded with this atom.
       *
       * @param widget The widget the angle-sectors will be drawn on.
       * @param atom The atom whose angles are being drawn.
       */
      void drawAtomAngles(GLWidget *widget, Atom *atom);

      /**
       * Draws sectors around a given atom representing the angles between neighbouring
       * atoms bonded with this atom and an atom bonded to this atom by a given bond.
       *
       * @param widget The widget the angle-sectors will be drawn on.
       * @param atom The atom whose angles are being drawn.
       * @param bond The bond attached to the given atom that will be used as a reference
       *             point for all the angles.
       *
       * @pre The given atom must be either the begin or end atom of the given bond.
       */
      void drawAngles(GLWidget *widget, Atom *atom, Bond *bond);

      /**
       * Draws sectors around the root Atom of a given SkeletonTree based on the root
       * Bond of the tree and whether or not adjacent Atoms form a part of the skeleton
       * or not.
       *
       * @param widget The widget the angle-sectors will be drawn on.
       * @param skeleton The SkeletonTree whose root Atom's angles are to be drawn.
       */
      void drawSkeletonAngles(GLWidget *widget, SkeletonTree *skeleton);

      /**
       * Calculates whether the manipulation plane is close enough to any atoms (that
       * are 1 bond away from either of the atoms attached to the given bond) to
       * 'snap-to' them.
       *
       * NOTE: Any atoms that lie along the same line as the bond are disregarded in
       * the calculations otherwise the plane would always try snap-to them as their
       * angle is 0.
       *
       * @param widget The widget the molecule and construction plane are on.
       * @param bond The bond through which the manipulation plane lies.
       * @param referencePoint The current reference point that defines the manipulation
       *                       plane.
       * @param maximumAngle The maximum angle between the current reference point
       *                     and any atom that determines whether or not the plane is
       *                     close enough to snap-to the atom.
       *
       * @return A vector representing the closest Atom to the manipulation plane, to
       *         be used as the reference point for drawing the plane, if any atom is
       *         close enough.  If no atom is close enough to 'snap-to', NULL is
       *         returned.
       */
      Eigen::Vector3d* calculateSnapTo(GLWidget *widget, Bond *bond, 
                                       Eigen::Vector3d *referencePoint, 
                                       double maximumAngle);

      /**
       * Draws a rectangle through a bond that can be used as a construction plane to
       * manipulate the bond itself, or the atoms at either end of the bond.
       *
       * @param widget The widget the rectangle will be drawn on.
       * @param bond The bond through which the rectangle will be drawn.
       * @param referencePoint A point orthagonal to the bond that defines the plane
       *                       the rectangle will be drawn on.
       * @param rgb An array of doubles representing the red/green/blue values of the
       *            color for the rectangle.
       */
      void drawManipulationRectangle(GLWidget *widget, Bond *bond, 
                                     Eigen::Vector3d *referencePoint, double rgb[3]);

      /**
       * Draws a sphere of a given radius around a given vector.
       *
       * @param widget The widget the sphere will be drawn on.
       * @param center The center of the sphere.
       * @param radius The radius of the sphere.
       * @param alpha The alpha value that determines the opacity of the sphere.
       */
      void drawSphere(GLWidget *widget, const Eigen::Vector3d &center, double radius,
                      float alpha);
      //@}

      /**
       * Connects this tool to the widget's ToolGroup so as to detect the signal
       * emitted when the tool changes.
       *
       * @param widget The GLWidget containing the ToolGroup being connected to.
       * @param toolGroup A pointer that will be (or is already) set to the
       *                  ToolGroup being connected to.
       */
      void connectToolGroup(GLWidget *widget, ToolGroup *toolGroup);

      /**
       * Clears any data and frees up any memory that is used by the tool.  This
       * procedure should be used when the tool is changed, the molecule cleared,
       * or the program exits etc.
       */
      void clearData();

      /**
       * Performs a rotation on a vector.
       *
       * @param angle The angle to rotate by in radians.
       * @param rotationVector The Vector3d to rotate around, must be a unit vector.
       * @param centerVector The Vector3d postion around which to rotate.
       * @param postionVector The Vector3d postion of the vector to rotate.
       *
       * @return A Vector3d with the final postion after the rotation is performed.
       *
       * @pre rotationVector must be a unit vector (of length 1).
       */
      Eigen::Vector3d performRotation(double angle, Eigen::Vector3d rotationVector,
                                      Eigen::Vector3d centerVector,
                                      Eigen::Vector3d positionVector);

    private Q_SLOTS:
     
       /**
       * Function to be called when the tool is changed.
       *
       * @param tool The currently selected tool.
       */
      void toolChanged(Tool* tool);

      /**
       * Functnion to be called when the molecule is changed.
       *
       * @param previous The previous Molecule.
       * @param next The new Molecule.
       */
      void moleculeChanged(Molecule* previous, Molecule* next);

      /**
       * Function to be called when a primitive is removed.
       *
       * @param primitive The primitive that was removed.
       */
      void primitiveRemoved(Primitive* primitive);

      /**
       * Function to be called when the settings widget is destroyed.
       */
      void settingsWidgetDestroyed();

  };

  /**
   * @class WiiMoteMoveCommand
   * @brief An implementation of QUndoCommand used to undo bond centric manipulations.
   * @author Shahzad Ali, Ross Braithwaite, James Bunt
   *
   * This class is an implementation of QUndoCommand that can be used to allow
   * the two types of bond-centric manipulation to be undone.  These two types
   * of manipulation are:
   *  - Adjusting bond length.
   *  - Adjusting bond angles.
   */
  class WiiMoteMoveCommand : public QUndoCommand
  {
    public:
      //!Constructor
      /**
       * Creates an undo command.
       *
       * @param molecule The molecule to store for undoing.
       * @param parent The parent undo command, or nothing.
       */
          WiiMoteMoveCommand(Molecule *molecule, QUndoCommand *parent = 0);

      //!Constructor
      /**
       * Creates an undo command.
       *
       * @param molecule The molecule to store for undoing.
       * @param atom The atom that has been moved.
       * @param pos The new position of the atom.
       * @param parent The parent undo command or null.
       */
          WiiMoteMoveCommand(Molecule *molecule, Atom *atom,
                             Eigen::Vector3d pos, QUndoCommand *parent = 0);

      /**
       * Redo move.
       */
      void redo();

      /**
       * Undo move
       */
      void undo();

      /**
       * returns if undo commands are merged together to one command.
       *
       * @param command undo command to merge
       *
       * @return false
       */
      bool mergeWith(const QUndoCommand * command);

      /**
       * returns id of this undo command
       *
       * @return id of this undo command
       */
      int id() const;

    private:
      Molecule m_moleculeCopy;
      Molecule *m_molecule;
      int m_atomIndex;
      Eigen::Vector3d m_pos;
      bool undone;
  };


  class WiiMoteToolFactory : public QObject, public ToolFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::ToolFactory)

    public:
      Tool *createInstance(QObject *parent = 0) { 
            return new WiiMoteTool(parent); }
  };
/*
  * Copyright (c) 1999-2002 Vojtech Pavlik
  *
  * This program is free software; you can redistribute it and/or modify it
  * under the terms of the GNU General Public License version 2 as published by
  * the Free Software Foundation.
 */

#define KEY_ESC   1
#define KEY_1     2
#define KEY_2     3
#define KEY_3     4
#define KEY_4     5
#define KEY_5     6
#define KEY_6     7
#define KEY_7     8
#define KEY_8     9
#define KEY_9     10
#define KEY_0     11
#define KEY_MINUS   12
#define KEY_EQUAL   13
#define KEY_BACKSPACE   14
#define KEY_TAB     15
#define KEY_Q     16
#define KEY_W     17
#define KEY_E     18
#define KEY_R     19
#define KEY_T     20
#define KEY_Y     21
#define KEY_U     22
#define KEY_I     23
#define KEY_O     24
#define KEY_P     25
#define KEY_LEFTBRACE   26
#define KEY_RIGHTBRACE    27
#define KEY_ENTER   28
#define KEY_LEFTCTRL    29
#define KEY_A     30
#define KEY_S     31
#define KEY_D     32
#define KEY_F     33
#define KEY_G     34
#define KEY_H     35
#define KEY_J     36
#define KEY_K     37
#define KEY_L     38
#define KEY_SEMICOLON   39
#define KEY_APOSTROPHE    40
#define KEY_GRAVE   41
#define KEY_LEFTSHIFT   42
#define KEY_BACKSLASH   43
#define KEY_Z     44
#define KEY_X     45
#define KEY_C     46
#define KEY_V     47
#define KEY_B     48
#define KEY_N     49
#define KEY_M     50
#define KEY_COMMA   51
#define KEY_DOT     52
#define KEY_SLASH   53
#define KEY_RIGHTSHIFT    54
#define KEY_KPASTERISK    55
#define KEY_LEFTALT   56
#define KEY_SPACE   57
#define KEY_CAPSLOCK    58
#define KEY_F1      59
#define KEY_F2      60
#define KEY_F3      61
#define KEY_F4      62
#define KEY_F5      63
#define KEY_F6      64
#define KEY_F7      65
#define KEY_F8      66
#define KEY_F9      67
#define KEY_F10     68
#define KEY_NUMLOCK   69
#define KEY_SCROLLLOCK    70
#define KEY_KP7     71
#define KEY_KP8     72
#define KEY_KP9     73
#define KEY_KPMINUS   74
#define KEY_KP4     75
#define KEY_KP5     76
#define KEY_KP6     77
#define KEY_KPPLUS    78
#define KEY_KP1     79
#define KEY_KP2     80
#define KEY_KP3     81
#define KEY_KP0     82
#define KEY_KPDOT   83

#define KEY_ZENKAKUHANKAKU  85
#define KEY_102ND   86
#define KEY_F11     87
#define KEY_F12     88
#define KEY_RO      89
#define KEY_KATAKANA    90
#define KEY_HIRAGANA    91
#define KEY_HENKAN    92
#define KEY_KATAKANAHIRAGANA  93
#define KEY_MUHENKAN    94
#define KEY_KPJPCOMMA   95
#define KEY_KPENTER   96
#define KEY_RIGHTCTRL   97
#define KEY_KPSLASH   98
#define KEY_SYSRQ   99
#define KEY_RIGHTALT    100
#define KEY_LINEFEED    101
#define KEY_HOME    102
#define KEY_UP      103
#define KEY_PAGEUP    104
#define KEY_LEFT    105
#define KEY_RIGHT   106
#define KEY_END     107
#define KEY_DOWN    108
#define KEY_PAGEDOWN    109
#define KEY_INSERT    110
#define KEY_DELETE    111
#define KEY_MACRO   112
#define KEY_MUTE    113
#define KEY_VOLUMEDOWN    114
#define KEY_VOLUMEUP    115
#define KEY_POWER   116
#define KEY_KPEQUAL   117
#define KEY_KPPLUSMINUS   118
#define KEY_PAUSE   119

#define KEY_KPCOMMA   121
#define KEY_HANGUEL   122
#define KEY_HANJA   123
#define KEY_YEN     124
#define KEY_LEFTMETA    125
#define KEY_RIGHTMETA   126
#define KEY_COMPOSE   127

#define KEY_STOP    128
#define KEY_AGAIN   129
#define KEY_PROPS   130
#define KEY_UNDO    131
#define KEY_FRONT   132
#define KEY_COPY    133
#define KEY_OPEN    134
#define KEY_PASTE   135
#define KEY_FIND    136
#define KEY_CUT     137
#define KEY_HELP    138
#define KEY_MENU    139
#define KEY_CALC    140
#define KEY_SETUP   141
#define KEY_SLEEP   142
#define KEY_WAKEUP    143
#define KEY_FILE    144
#define KEY_SENDFILE    145
#define KEY_DELETEFILE    146
#define KEY_XFER    147
#define KEY_PROG1   148
#define KEY_PROG2   149
#define KEY_WWW     150
#define KEY_MSDOS   151
#define KEY_COFFEE    152
#define KEY_DIRECTION   153
#define KEY_CYCLEWINDOWS  154
#define KEY_MAIL    155
#define KEY_BOOKMARKS   156
#define KEY_COMPUTER    157
#define KEY_BACK    158
#define KEY_FORWARD   159
#define KEY_CLOSECD   160
#define KEY_EJECTCD   161
#define KEY_EJECTCLOSECD  162
#define KEY_NEXTSONG    163
#define KEY_PLAYPAUSE   164
#define KEY_PREVIOUSSONG  165
#define KEY_STOPCD    166
#define KEY_RECORD    167
#define KEY_REWIND    168
#define KEY_PHONE   169
#define KEY_ISO     170
#define KEY_CONFIG    171
#define KEY_HOMEPAGE    172
#define KEY_REFRESH   173
#define KEY_EXIT    174
#define KEY_MOVE    175
#define KEY_EDIT    176
#define KEY_SCROLLUP    177
#define KEY_SCROLLDOWN    178
#define KEY_KPLEFTPAREN   179
#define KEY_KPRIGHTPAREN  180
#define KEY_NEW     181
#define KEY_REDO    182

#define KEY_F13     183
#define KEY_F14     184
#define KEY_F15     185
#define KEY_F16     186
#define KEY_F17     187
#define KEY_F18     188
#define KEY_F19     189
#define KEY_F20     190
#define KEY_F21     191
#define KEY_F22     192
#define KEY_F23     193
#define KEY_F24     194

#define KEY_PLAYCD    200
#define KEY_PAUSECD   201
#define KEY_PROG3   202
#define KEY_PROG4   203
#define KEY_SUSPEND   205
#define KEY_CLOSE   206
#define KEY_PLAY    207
#define KEY_FASTFORWARD   208
#define KEY_BASSBOOST   209
#define KEY_PRINT   210
#define KEY_HP      211
#define KEY_CAMERA    212
#define KEY_SOUND   213
#define KEY_QUESTION    214
#define KEY_EMAIL   215
#define KEY_CHAT    216
#define KEY_SEARCH    217
#define KEY_CONNECT   218
#define KEY_FINANCE   219
#define KEY_SPORT   220
#define KEY_SHOP    221
#define KEY_ALTERASE    222
#define KEY_CANCEL    223
#define KEY_BRIGHTNESSDOWN  224
#define KEY_BRIGHTNESSUP  225
#define KEY_MEDIA   226

#define KEY_SWITCHVIDEOMODE 227
#define KEY_KBDILLUMTOGGLE  228
#define KEY_KBDILLUMDOWN  229
#define KEY_KBDILLUMUP    230

#define KEY_SEND    231
#define KEY_REPLY   232
#define KEY_FORWARDMAIL   233
#define KEY_SAVE    234
#define KEY_DOCUMENTS   235

#define KEY_BATTERY   236

#define KEY_UNKNOWN   240

#define BTN_MISC    0x100
#define BTN_0     0x100
#define BTN_1     0x101
#define BTN_2     0x102
#define BTN_3     0x103
#define BTN_4     0x104
#define BTN_5     0x105
#define BTN_6     0x106
#define BTN_7     0x107
#define BTN_8     0x108
#define BTN_9     0x109

#define BTN_MOUSE   0x110
#define BTN_LEFT    0x110
#define BTN_RIGHT   0x111
#define BTN_MIDDLE    0x112
#define BTN_SIDE    0x113
#define BTN_EXTRA   0x114
#define BTN_FORWARD   0x115
#define BTN_BACK    0x116
#define BTN_TASK    0x117

#define BTN_JOYSTICK    0x120
#define BTN_TRIGGER   0x120
#define BTN_THUMB   0x121
#define BTN_THUMB2    0x122
#define BTN_TOP     0x123
#define BTN_TOP2    0x124
#define BTN_PINKIE    0x125
#define BTN_BASE    0x126
#define BTN_BASE2   0x127
#define BTN_BASE3   0x128
#define BTN_BASE4   0x129
#define BTN_BASE5   0x12a
#define BTN_BASE6   0x12b
#define BTN_DEAD    0x12f

#define BTN_GAMEPAD   0x130
#define BTN_A     0x130
#define BTN_B     0x131
#define BTN_C     0x132
#define BTN_X     0x133
#define BTN_Y     0x134
#define BTN_Z     0x135
#define BTN_TL      0x136
#define BTN_TR      0x137
#define BTN_TL2     0x138
#define BTN_TR2     0x139
#define BTN_SELECT    0x13a
#define BTN_START   0x13b
#define BTN_MODE    0x13c
#define BTN_THUMBL    0x13d
#define BTN_THUMBR    0x13e

#define BTN_DIGI    0x140
#define BTN_TOOL_PEN    0x140
#define BTN_TOOL_RUBBER   0x141
#define BTN_TOOL_BRUSH    0x142
#define BTN_TOOL_PENCIL   0x143
#define BTN_TOOL_AIRBRUSH 0x144
#define BTN_TOOL_FINGER   0x145
#define BTN_TOOL_MOUSE    0x146
#define BTN_TOOL_LENS   0x147
#define BTN_TOUCH   0x14a
#define BTN_STYLUS    0x14b
#define BTN_STYLUS2   0x14c
#define BTN_TOOL_DOUBLETAP  0x14d
#define BTN_TOOL_TRIPLETAP  0x14e

#define BTN_WHEEL   0x150
#define BTN_GEAR_DOWN   0x150
#define BTN_GEAR_UP   0x151

#define KEY_OK      0x160
#define KEY_SELECT    0x161
#define KEY_GOTO    0x162
#define KEY_CLEAR   0x163
#define KEY_POWER2    0x164
#define KEY_OPTION    0x165
#define KEY_INFO    0x166
#define KEY_TIME    0x167
#define KEY_VENDOR    0x168
#define KEY_ARCHIVE   0x169
#define KEY_PROGRAM   0x16a
#define KEY_CHANNEL   0x16b
#define KEY_FAVORITES   0x16c
#define KEY_EPG     0x16d
#define KEY_PVR     0x16e
#define KEY_MHP     0x16f
#define KEY_LANGUAGE    0x170
#define KEY_TITLE   0x171
#define KEY_SUBTITLE    0x172
#define KEY_ANGLE   0x173
#define KEY_ZOOM    0x174
#define KEY_MODE    0x175
#define KEY_KEYBOARD    0x176
#define KEY_SCREEN    0x177
#define KEY_PC      0x178
#define KEY_TV      0x179
#define KEY_TV2     0x17a
#define KEY_VCR     0x17b
#define KEY_VCR2    0x17c
#define KEY_SAT     0x17d
#define KEY_SAT2    0x17e
#define KEY_CD      0x17f
#define KEY_TAPE    0x180
#define KEY_RADIO   0x181
#define KEY_TUNER   0x182
#define KEY_PLAYER    0x183
#define KEY_TEXT    0x184
#define KEY_DVD     0x185
#define KEY_AUX     0x186
#define KEY_MP3     0x187
#define KEY_AUDIO   0x188
#define KEY_VIDEO   0x189
#define KEY_DIRECTORY   0x18a
#define KEY_LIST    0x18b
#define KEY_MEMO    0x18c
#define KEY_CALENDAR    0x18d
#define KEY_RED     0x18e
#define KEY_GREEN   0x18f
#define KEY_YELLOW    0x190
#define KEY_BLUE    0x191
#define KEY_CHANNELUP   0x192
#define KEY_CHANNELDOWN   0x193
#define KEY_FIRST   0x194
#define KEY_LAST    0x195
#define KEY_AB      0x196
#define KEY_NEXT    0x197
#define KEY_RESTART   0x198
#define KEY_SLOW    0x199
#define KEY_SHUFFLE   0x19a
#define KEY_BREAK   0x19b
#define KEY_PREVIOUS    0x19c
#define KEY_DIGITS    0x19d
#define KEY_TEEN    0x19e
#define KEY_TWEN    0x19f

#define KEY_DEL_EOL   0x1c0
#define KEY_DEL_EOS   0x1c1
#define KEY_INS_LINE    0x1c2
#define KEY_DEL_LINE    0x1c3

#define KEY_FN      0x1d0
#define KEY_FN_ESC    0x1d1
#define KEY_FN_F1   0x1d2
#define KEY_FN_F2   0x1d3
#define KEY_FN_F3   0x1d4
#define KEY_FN_F4   0x1d5
#define KEY_FN_F5   0x1d6
#define KEY_FN_F6   0x1d7
#define KEY_FN_F7   0x1d8
#define KEY_FN_F8   0x1d9
#define KEY_FN_F9   0x1da
#define KEY_FN_F10    0x1db
#define KEY_FN_F11    0x1dc
#define KEY_FN_F12    0x1dd
#define KEY_FN_1    0x1de
#define KEY_FN_2    0x1df
#define KEY_FN_D    0x1e0
#define KEY_FN_E    0x1e1
#define KEY_FN_F    0x1e2
#define KEY_FN_S    0x1e3
#define KEY_FN_B    0x1e4

#define KEY_BRL_DOT1    0x1f1
#define KEY_BRL_DOT2    0x1f2
#define KEY_BRL_DOT3    0x1f3
#define KEY_BRL_DOT4    0x1f4
#define KEY_BRL_DOT5    0x1f5
#define KEY_BRL_DOT6    0x1f6
#define KEY_BRL_DOT7    0x1f7
#define KEY_BRL_DOT8    0x1f8

/*
  * Relative axes
 */

#define REL_X     0x00
#define REL_Y     0x01
#define REL_Z     0x02
#define REL_RX      0x03
#define REL_RY      0x04
#define REL_RZ      0x05
#define REL_HWHEEL    0x06
#define REL_DIAL    0x07
#define REL_WHEEL   0x08
#define REL_MISC    0x09

/*
  * Absolute axes
 */

#define ABS_X     0x00
#define ABS_Y     0x01
#define ABS_Z     0x02
#define ABS_RX      0x03
#define ABS_RY      0x04
#define ABS_RZ      0x05
#define ABS_THROTTLE    0x06
#define ABS_RUDDER    0x07
#define ABS_WHEEL   0x08
#define ABS_GAS     0x09
#define ABS_BRAKE   0x0a
#define ABS_HAT0X   0x10
#define ABS_HAT0Y   0x11
#define ABS_HAT1X   0x12
#define ABS_HAT1Y   0x13
#define ABS_HAT2X   0x14
#define ABS_HAT2Y   0x15
#define ABS_HAT3X   0x16
#define ABS_HAT3Y   0x17
#define ABS_PRESSURE    0x18
#define ABS_DISTANCE    0x19
#define ABS_TILT_X    0x1a
#define ABS_TILT_Y    0x1b
#define ABS_TOOL_WIDTH    0x1c
#define ABS_VOLUME    0x20
#define ABS_MISC    0x28

////////Wii Mote Buttons//////////
#define CONF_WM_BTN_UP    0
#define CONF_WM_BTN_DOWN  1
#define CONF_WM_BTN_LEFT  2
#define CONF_WM_BTN_RIGHT 3
#define CONF_WM_BTN_A   4
#define CONF_WM_BTN_B   5
#define CONF_WM_BTN_MINUS 6
#define CONF_WM_BTN_PLUS  7
#define CONF_WM_BTN_HOME  8
#define CONF_WM_BTN_1   9
#define CONF_WM_BTN_2   10

///////UINPUT///////////////////
#define UINPUT_NAME   "Nintendo Wiimote"
#ifdef BUS_BLUETOOTH
#define UINPUT_BUSTYPE  BUS_BLUETOOTH
#else
#define UINPUT_BUSTYPE BUS_USB
#endif
#define UINPUT_VENDOR 0x0001
#define UINPUT_PRODUCT  0x0001
#define UINPUT_VERSION  0x0001

} // end namespace Avogadro


#endif /*__WIIMOTETOOL_H*/
