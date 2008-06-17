
#include "animatemolextension.h"

#include <QAction>

using namespace std;
using namespace OpenBabel;

namespace Avogadro
{

  // this is a trick to identify what action we are taking
  enum AnimateMolExtensionIndex {
    FirstAction = 0,
    SecondAction
  };

  AnimateMolExtension::AnimateMolExtension( QObject *parent ) : Extension( parent ), m_molecule(0), m_animateMolDialog(0), m_timeLine(0), m_frameCount(0)
  {
    
    QAction *action = new QAction(this);
    action->setText(tr("Animate Mol..."));
    m_actions.append(action);
    
    action = new QAction( this );
    action->setSeparator(true);
    m_actions.append(action);

  }

  AnimateMolExtension::~AnimateMolExtension()
  {
  }

  QList<QAction *> AnimateMolExtension::actions() const
  {
    return m_actions;
  }

  // allows us to set the intended menu path for each action
  QString AnimateMolExtension::menuPath(QAction *action) const
  {
    int i = action->data().toInt();

    switch ( i ) {
      case FirstAction:
        return tr("&Extensions") + '>' + tr("&Animate Mol Extension");
        break;
      case SecondAction:
        return tr("&Edit") + '>' + tr("&Animate Mol Extension");
        break;
    }
    return "";
  }

  QDockWidget * AnimateMolExtension::dockWidget()
  {
    // if we need a dock widget we can set one here
    return 0;
  }

  void AnimateMolExtension::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QUndoCommand* AnimateMolExtension::performAction(QAction *, GLWidget*)
  {
    qDebug( "AnimateMolExtension::performAction()" );

    if (!m_animateMolDialog)
      {
	m_timeLine = new QTimeLine;
	m_animateMolDialog = new AnimateMolDialog;
	
	connect(m_animateMolDialog, SIGNAL(fileName(QString)), this, SLOT(loadFile(QString)));
	connect(m_animateMolDialog, SIGNAL(sliderChanged(int)), this, SLOT(setFrame(int)));
	connect(m_animateMolDialog, SIGNAL(fpsChanged(int)), this, SLOT(setDuration(int)));
	connect(m_animateMolDialog, SIGNAL(loopChanged(int)), this, SLOT(setLoop(int)));
	
	connect(m_timeLine, SIGNAL(frameChanged(int)), this, SLOT(setFrame(int)));
	connect(m_animateMolDialog, SIGNAL(play()), m_timeLine, SLOT(start()));
	connect(m_animateMolDialog, SIGNAL(pause()), m_timeLine, SLOT(stop()));
	connect(m_animateMolDialog, SIGNAL(stop()), this, SLOT(stop()));
	
	} 
    m_animateMolDialog->show();
    return 0;
  }


  void AnimateMolExtension::setDuration(int i)
  {
    int interval = 1000 / i;
    m_timeLine->setUpdateInterval(interval);
    int duration = interval * m_frameCount;
    m_timeLine->setDuration(duration);
  }

  void AnimateMolExtension::setLoop(int state)
  {
    if (state == Qt::Checked) {
      m_timeLine->setLoopCount(0);
    } else {
      m_timeLine->setLoopCount(1);
    }
  }

  void AnimateMolExtension::setFrame(int i)
  {
    // if (m_timeLine->state() != QTimeLine::Running)
    //  m_timeLine->setCurrentTime(m_timeLine->updateInterval() * i);

    m_animateMolDialog->setFrame(i);
    m_molecule->SetConformer(i - 1);
    m_molecule->update();
  }

  void AnimateMolExtension::stop()
  {
    m_timeLine->stop();
    m_timeLine->setCurrentTime(0);
    setFrame(1);
  }
}

#include "animatemolextension.moc"


