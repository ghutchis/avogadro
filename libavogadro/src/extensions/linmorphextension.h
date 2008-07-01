
#ifndef LINMORPHEXTENSION_H
#define LINMORPHEXTENSION_H


#include <avogadro/extension.h>
#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>
#include <QTimeLine>
//#include "animatemoldialog.h"
#include "linmorphdialog.h"

#include <avogadro/primitive.h>
#include <avogadro/glwidget.h>

namespace Avogadro {
  
  class LinMorphExtension : public Extension
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

      
      //!the current frame?
      int m_frameCount;




      virtual QList<QAction *> actions() const;


      virtual QDockWidget * dockWidget();
      virtual QUndoCommand* performAction(QAction *action, GLWidget *widget);

      virtual void setMolecule(Molecule *molecule);
      


  protected:
      Molecule *m_molecule;
      Molecule *m_secondMolecule;
      
      // hate hate hate to put this here, but not sure what else to do
      // this will be the glwidget passed from performAction
      GLWidget *m_widget;
      

      QList<QAction *> m_actions;
      LinMorphDialog *m_animateMolDialog;
      QTimeLine *m_timeLine;

    
   protected Q_SLOTS:
      void saveGlSnapshots(QString prefix);
      void savePovSnapshots(QString prefix);
      void setDuration(int i);
      void setLoop(int state);
      void setFrame(int i);
      void setFrameCount(int i);
      void stop();
      virtual void loadFile(QString file);

      
  private:
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
