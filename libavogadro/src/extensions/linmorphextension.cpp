
#include "linmorphextension.h"

#include <QAction>

#include <QMessageBox>

#include <openbabel/obconversion.h>
#include <avogadro/povpainter.h>
#include <QInputDialog>

using namespace std;
using namespace OpenBabel;

namespace Avogadro
{

  //  LinMorphExtension::LinMorphExtension( QObject *parent ) : AnimateMolExtension( parent ), m_numFrames(100)  
  LinMorphExtension::LinMorphExtension( QObject *parent ) :Extension( parent ), m_molecule(0), m_animateMolDialog(0), m_timeLine(0), m_frameCount(100) 
  {
    
    QAction *action = new QAction(this);
    action->setText(tr("Lin Morph..."));
    m_actions.append(action);
    
    action = new QAction( this );
    action->setSeparator(true);
    m_actions.append(action);
  }

  LinMorphExtension::~LinMorphExtension()
  {
    if (m_secondMolecule)
      delete m_secondMolecule;
  }

  QList<QAction *> LinMorphExtension::actions() const
  {
    return m_actions;
  }


  void LinMorphExtension::loadFile(QString file)
  {
    qDebug("LinMorphExtension::loadFile()");

    if (file.isEmpty())
      return;

    //if (!m_secondMolecule) 
      m_secondMolecule = new Molecule;
    //Molecule* mTemp = new Molecule;
    OBConversion conv;
    OBFormat *inFormat = OBConversion::FormatFromExt(( file.toAscii() ).data() );        
    if ( !inFormat || !conv.SetInFormat( inFormat ) ) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Cannot read file format of file %1." )
          .arg( file ) );
      return;
    }

    if (!conv.ReadFile(m_secondMolecule, file.toStdString())) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
	  tr( "Read mol file %1 failed." )
          .arg( file ) );
      return;
    }

    qDebug("LinMorphExtension::loadFile complete");
    computeConformers(m_secondMolecule);

    
    //    m_frameCount = m_molecule->NumConformers();
    //m_animateMolDialog->setFrameCount(m_frameCount);
    m_animateMolDialog->setFrame(1);
    m_timeLine->setFrameRange(1, m_frameCount);
    setDuration(m_animateMolDialog->fps());

   

  }

  // allows us to set the intended menu path for each action
  QString LinMorphExtension::menuPath(QAction *) const
  {
    return tr("&Extensions");

  }


     //! compute the conformers and set in molecule
  void LinMorphExtension::computeConformers(Molecule* conformer2Mol){
    
    //copied from OBForceField::SetConfor
    int k,l;
    vector<double*> conf;
    double* xyz = NULL;

    if (conformer2Mol->NumAtoms() != m_molecule->NumAtoms()) {
      
      QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Two molecules have different number atoms %1 %2" )
			    .arg(m_molecule->NumAtoms()).arg( conformer2Mol->NumAtoms()));
      return;
  }


    double* initCoords = m_molecule->GetCoordinates();
    double* finalCoords = conformer2Mol->GetCoordinates();


    for (k=0 ; k < m_frameCount; ++k) {
      xyz = new double [3*m_molecule->NumAtoms()];
      for (l=0 ; l<(int) (3*m_molecule->NumAtoms()) ; ++l) {
	double dCoord = (finalCoords[l] - initCoords[l])/m_frameCount;
	xyz[l]=initCoords[l] + dCoord*k;
      }
      conf.push_back(xyz);
    }
    
    m_molecule->SetConformers(conf);
        
    qDebug("Num conformers set: %d", m_molecule->NumConformers());
  }


  QDockWidget * LinMorphExtension::dockWidget()
  {
    // if we need a dock widget we can set one here
    return 0;
  }

  void LinMorphExtension::setMolecule(Molecule *molecule)
  {
    m_molecule = molecule;
  }

  QUndoCommand* LinMorphExtension::performAction(QAction *, GLWidget* widget)
  {
    qDebug( "LinMorphExtension::performAction()" );
    
    m_widget = widget;

    if (!m_animateMolDialog)
      {
	m_timeLine = new QTimeLine;
	m_animateMolDialog = new LinMorphDialog;
	
	connect(m_animateMolDialog, SIGNAL(fileName(QString)), this, SLOT(loadFile(QString)));
	connect(m_animateMolDialog, SIGNAL(snapshotsPrefix(QString)), this, SLOT(savePovSnapshots(QString)));
	connect(m_animateMolDialog, SIGNAL(sliderChanged(int)), this, SLOT(setFrame(int)));
	connect(m_animateMolDialog, SIGNAL(fpsChanged(int)), this, SLOT(setDuration(int)));
	connect(m_animateMolDialog, SIGNAL(loopChanged(int)), this, SLOT(setLoop(int)));
	
	connect(m_timeLine, SIGNAL(frameChanged(int)), this, SLOT(setFrame(int)));
	


	connect(m_animateMolDialog, SIGNAL(play()), m_timeLine, SLOT(start()));
	connect(m_animateMolDialog, SIGNAL(pause()), m_timeLine, SLOT(stop()));
	connect(m_animateMolDialog, SIGNAL(stop()), this, SLOT(stop()));
	connect(m_animateMolDialog, SIGNAL(frameCountChanged(int)), this, SLOT(setFrameCount(int)));
	
	} 
    m_animateMolDialog->show();


        
    return 0;
  } 


  void LinMorphExtension::setDuration(int i)
  {
    int interval = 1000 / i;
    m_timeLine->setUpdateInterval(interval);
    int duration = interval * m_frameCount;
    m_timeLine->setDuration(duration);
  }


  void LinMorphExtension::setFrameCount(int i)
  {
    m_frameCount = i;
    m_animateMolDialog->setFrameCount(i);
    if (m_secondMolecule)
      computeConformers(m_secondMolecule);
  }

  void LinMorphExtension::setLoop(int state)
  {
    if (state == Qt::Checked) {
      m_timeLine->setLoopCount(0);
    } else {
      m_timeLine->setLoopCount(1);
    }
  }

  void LinMorphExtension::setFrame(int i)
  {
    // if (m_timeLine->state() != QTimeLine::Running)
    //  m_timeLine->setCurrentTime(m_timeLine->updateInterval() * i);

    m_animateMolDialog->setFrame(i);
    m_molecule->SetConformer(i - 1);
    m_molecule->update();
  }

  void LinMorphExtension::stop()
  {
    m_timeLine->stop();
    m_timeLine->setCurrentTime(0);
    setFrame(1);
  }


  
  void LinMorphExtension::saveGlSnapshots(QString prefix)
  {
    // This function does not work.  NKF - 6/30/2008 
    if (!m_widget) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "GL widget was not correctly initialized in order to save snapshots" ));
      return;
    }
      
    //computeConformers(m_secondMolecule);
    for (int i=1; i<=m_frameCount; i++) {
      setFrame(i);
      QImage exportImage = m_widget->grabFrameBuffer( true );
      QString ssfileName = prefix + QString::number(i) + ".png";
      if ( !exportImage.save( ssfileName ) ) {
	QMessageBox::warning( NULL, tr( "Avogadro" ),
			      tr( "Cannot save file %1." ).arg( ssfileName ) );
	return;
      }
    }
  }


  void LinMorphExtension::savePovSnapshots(QString prefix)
  {
    // use the current glWidge for things like camera
    if (!m_widget) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "GL widget was not correctly initialized in order to save snapshots" ));
      return;
    }

    computeConformers(m_secondMolecule);
    if (m_frameCount != m_molecule->NumConformers()){
      QMessageBox::warning( NULL, tr( "Avogadro" ),
			    tr( "m_frameCount != numConformers" ) );
      return;
    }
      
  //
    bool ok;
    int w = m_widget->width();
    int h = m_widget->height();
    double defaultAspectRatio = static_cast<double>(w)/h;
    double aspectRatio =
      QInputDialog::getDouble(0,
			      QObject::tr("Set Aspect Ratio"),
			      QObject::tr("The current Avogadro scene is %1x%2 pixels large, "
					  "and therefore has aspect ratio %3.\n"
					  "You may keep this value, for example if you "
					  "intend to use POV-Ray\n"
					  "to produce an image of %4x1000 pixels, "
					  "or you may enter any other positive value,\n"
					  "for example 1 if you intend to use POV-Ray to "
					  "produce a square image, "
					  "like 1000x1000 pixels.")
			      .arg(w).arg(h).arg(defaultAspectRatio)
			      .arg(static_cast<int>(1000*defaultAspectRatio)),
			      defaultAspectRatio,
			      0.1,
			      10,
			      6,
			      &ok);
    
    for (int i=1; i<=m_frameCount; i++) {
      
      setFrame(i);
      QString ssfileName = prefix + QString::number(i) + ".pov";
      
            
      if(ok)
	POVPainterDevice pd( ssfileName, aspectRatio, m_widget );
      
    }
  }

    
  void LinMorphExtension::saveTrajFile(QString file) 
  {
    OBConversion conv;
    conv.FormatFromExt((file.toAscii() ).data());
    if (!conv.WriteFile(m_molecule, file.toStdString())) {
       QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Write trajectory file %1 failed." )
          .arg( file ) );
    }
  }


}
#include "linmorphextension.moc"

Q_EXPORT_PLUGIN2(linmorphextension, Avogadro::LinMorphExtensionFactory)
