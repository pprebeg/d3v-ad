import os.path

from PySide6.QtWidgets import QApplication, QMenu, QFormLayout,QWidget,QHeaderView,QSlider,QLineEdit
from PySide6.QtWidgets import QDialog, QPushButton,QGridLayout,QVBoxLayout,QHBoxLayout,QTableView,QTextEdit,QLabel
from PySide6.QtCore import Slot,Qt,QMetaObject,QCoreApplication
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QRect
from PySide6.QtGui import QColor, QPainter
#from PySide6.QtCharts import QtCharts
from PySide6.QtWidgets import QFileDialog,QDialogButtonBox,QProgressDialog
from PySide6.QtCore import SIGNAL,SLOT
from typing import Dict,List
from uuid import uuid4


try:
    # d3v imports
    from signals import Signals
    from commands import Command
    from iohandlers import IOHandler
    from core import geometry_manager as manager
    from core.geometry import Geometry
    # d3v-ad
    from aircraftcomponent import *
    from bramoflaw import BraMoFlyingWing

except BaseException as error:
    print('An exception occurred: {}'.format(error))
except:
    print('Unknown exception occurred during signals connection')

class AircraftDesignCommand(Command):
    def __init__(self):
        super().__init__()
        self._app = QApplication.instance()
        self._app.registerIOHandler(AircraftDesignImporter())
        self.aircraft_components:Dict[uuid4, AircraftComponent]={}
        self.selected_component=None
        self.active_component = None
        self._last_active_guid = None

        mb = self.mainwin.menuBar()
        self.menuMain = QMenu("&Aircraft Geometry")
        mb.addMenu(self.menuMain)

        self.menuMain.addSeparator()

        menuImport = self.menuMain.addAction("&Import")
        menuImport.triggered.connect(self.onImportComponent)
        menuExport = self.menuMain.addAction("&Export")
        menuExport.triggered.connect(self.onExportComponent)


        try:
            manager.selected_geometry_changed.connect(self.onSelectedGeometryChanged)
            manager.geometry_created.connect(self.onGeometryCreated)
            manager.geometry_removed.connect(self.onGeometryRemoved)
            manager.visible_geometry_changed.connect(self.onVisibleGeometryChanged)
        except BaseException as error:
            print('An exception occurred: {}'.format(error))
        except:
            print('Unknown exception occurred during signals connection')


    @Slot()
    def onImportComponent(self):
        fname = QFileDialog.getOpenFileName(self.mainwin,
                                            'Select CPACS file for import','../../../../examples/ad',
                                            "hull form files (*.xml *.obj *.stl *txt)")
        fname = fname[0]
        if fname != "":
            hfi = AircraftDesignImporter(True)
            hf = hfi.importGeometry(fname)
            if hf is not None:
                hf.emit_geometry_built()

    def onExportComponent(self):
        available_export = ""
        if isinstance(self.active_component,TiglComponent):
            available_export = "available hull form types (*.hgf *.obj *.stl)"
        elif isinstance(self.active_component, AircraftComponentFromMesh):
            available_export = "available hull form types (*.obj *.stl)"

        fname = QFileDialog.getSaveFileName(self.mainwin,
                                            'Export {0} form as'.format(self.active_component.name),
                                            '../../../../examples/kyrenia', available_export)
        fileName = fname[0]
        self.active_component.exportGeometry(fileName)


    def on_active_component_changed(self):
        pass


    @property
    def num_components(self):
        return len(self.aircraft_components)

    @Slot()
    def onVisibleGeometryChanged(self, visible:List[Geometry], loaded:List[Geometry], selected:List[Geometry]):
        for g in visible:
            pass

    @Slot()
    def onSelectedGeometryChanged(self, visible: List[Geometry], loaded: List[Geometry], selected: List[Geometry]):
        if len(selected) == 1:
            if isinstance(selected[0],AircraftComponent):
                self.selected_component = selected[0]
                self.active_component = self.selected_component
                self.on_active_component_changed()
        elif len(selected) == 0:
            self.selected_component = None

    @Slot()
    def onGeometryCreated(self, geometries:List[Geometry]):
        nlast=self.num_components
        for g in geometries:
            if isinstance(g,AircraftComponent):
                self.aircraft_components[g.guid]=g
        if nlast == 0 and self.num_components==1:
            self.active_component=geometries[0]
            self.on_active_component_changed()
        elif self._last_active_guid == geometries[0].guid:
            self.active_component=geometries[0]
            self.on_active_component_changed()


    @Slot()
    def onGeometryRemoved(self, geometries:List[Geometry]):
        for g in geometries:
            if isinstance(g, AircraftComponent):
                self.aircraft_components.pop(g.guid)
                if self.active_component is g:
                    self.active_component = None
                    self.on_active_component_changed()
                if self.referent_component is g:
                    self.referent_component = None
                    self.on_referent_component_changed()



    @property
    def app(self):
        return self._app

    @property
    def mainwin(self):
        return self.app.mainFrame

    @property
    def glwin(self):
        return self.mainwin.glWin

class AircraftDesignImporter(IOHandler):
    def __init__(self,force_import=False):
        super().__init__()
        self.force_import=force_import

    def importGeometry(self, fileName):
        if len(fileName) < 1:
            return
        filename_no_ext, file_extension = os.path.splitext(os.path.basename(fileName))
        ta=None

        if self.force_import:
            if file_extension == ".xml" or file_extension == ".xml":
                ta = TiglAircraft(fileName)
            elif file_extension == ".stl" or file_extension == ".obj":
                ta = AircraftComponentFromMesh(fileName,filename_no_ext)
            elif file_extension == ".txt":
                ta = BraMoFlyingWing(fileName,filename_no_ext)
        if ta is not None:
            return ta

    def exportGeometry(self, fileName, geometry2export):
        if isinstance(geometry2export,AircraftComponent):
            geometry2export.exportGeometry(fileName)
        om.write_mesh(geometry2export.mesh,fileName)
        pass

    def getExportFormats(self):
       return (".xml")

    def getImportFormats(self):
        return (".xml")



def createCommand():
    return AircraftDesignCommand()