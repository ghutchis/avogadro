
#include "linmorphextension.h"

#include <QAction>

#include <QMessageBox>

#include <openbabel/obconversion.h>


using namespace std;
using namespace OpenBabel;

namespace Avogadro
{

  // this is a trick to identify what action we are taking
  enum LinMorphExtensionIndex {
    FirstAction = 0,
    SecondAction
  };

  LinMorphExtension::LinMorphExtension( QObject *parent ) : AnimateMolExtension( parent ), m_numFrames(100)  
  {
    //do nothing
   

    // open a molecule file and compute the lin morph between that molecule
    // and the other one
        
    
  }

  LinMorphExtension::~LinMorphExtension()
  {
  }


  void LinMorphExtension::loadFile(QString file)
  {
    qDebug("LinMorphExtension::loadFile()");

    if (file.isEmpty())
      return;

    Molecule* mTemp = new Molecule;
    OBConversion conv;
    OBFormat *inFormat = OBConversion::FormatFromExt(( file.toAscii() ).data() );        
    if ( !inFormat || !conv.SetInFormat( inFormat ) ) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
          tr( "Cannot read file format of file %1." )
          .arg( file ) );
      return;
    }

    if (!conv.ReadFile(mTemp, file.toStdString())) {
      QMessageBox::warning( NULL, tr( "Avogadro" ),
	  tr( "Read mol file %1 failed." )
          .arg( file ) );
      return;
    }

    qDebug("LinMorphExtension::loadFile complete");
    computeConformers(mTemp);

    m_frameCount = m_molecule->NumConformers();
    m_animateMolDialog->setFrameCount(m_frameCount);
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


    for (k=0 ; k < m_numFrames; ++k) {
      xyz = new double [3*m_molecule->NumAtoms()];
      for (l=0 ; l<(int) (3*m_molecule->NumAtoms()) ; ++l) {
	double dCoord = (finalCoords[l] - initCoords[l])/m_numFrames;
	xyz[l]=initCoords[l] + dCoord*k;
      }
      conf.push_back(xyz);
    }
    
    m_molecule->SetConformers(conf);
        
    qDebug("Num conformers set: %d", m_molecule->NumConformers());
  }


}
#include "linmorphextension.moc"

Q_EXPORT_PLUGIN2(linmorphextension, Avogadro::LinMorphExtensionFactory)
