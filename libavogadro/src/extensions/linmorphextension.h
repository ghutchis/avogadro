
#ifndef LINMORPHEXTENSION_H
#define LINMORPHEXTENSION_H

#include "animatemolextension.h"
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

namespace Avogadro {
  
  class LinMorphExtension : public AnimateMolExtension
  {
    Q_OBJECT
      public:
      //! Constructor
      LinMorphExtension(QObject *parent=0);
      //! Deconstructor
      virtual ~LinMorphExtension();

      virtual QString name() const { return QObject::tr("Animate Mol Lin"); }
      virtual QString description() const { return QObject::tr("Lin Morph Extension."); };

      virtual QString menuPath(QAction *action) const;

      //! the number of frames to be animated... i.e. the timestep dt
      int m_numFrames;
      
      //!the current frame
      int m_frameCount;

      
      
  private:
      virtual void loadFile(QString file);
      virtual void computeConformers(Molecule* conformer2Mol);
      
            
  };

  class LinMorphExtensionFactory : public QObject, public ExtensionFactory
  {
    Q_OBJECT
    Q_INTERFACES(Avogadro::ExtensionFactory)

    public:
      Extension *createInstance(QObject *parent = 0) { return new LinMorphExtension(parent); }
  };

}
#endif
