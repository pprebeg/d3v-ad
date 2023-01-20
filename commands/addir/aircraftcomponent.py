from geometry_extend import GeometryExtension
import openmesh as om
import numpy as np
import os
from tigl3.tigl3wrapper import Tigl3, TiglBoolean
from tixi3.tixi3wrapper import Tixi3
from tigl3.configuration import CCPACSConfiguration,CTiglUIDManager,CCPACSConfigurationManager_get_instance
from tigl3.geometry import ITiglGeometricComponent, CNamedShape
from OCC.Core.TopoDS import TopoDS_Shape
import sys
from typing import Dict,List
from addir.topodsmesher import get_open_mesh_from_TopoDS_using_shape_tesselator

class AircraftComponent(GeometryExtension):
    def __init__(self, name=''):
        super().__init__(name)


    def regenerate_component(self):
        # to be overriden in extended class
        pass

    def exportGeometry(self, fileName):
        # to be overriden in extended class
        pass


class AircraftComponentFromMesh(AircraftComponent):
    def __init__(self, fileName, name='', translation=np.zeros(3)):
        self._filename = fileName
        self._translation = translation
        self._original_mesh = self.read_file()
        super().__init__(name)

    @staticmethod
    def export_component(fileName: str, component: AircraftComponent):
        filename_no_ext, file_extension = os.path.splitext(fileName)
        if file_extension == ".stl" or file_extension == ".obj":
            om.write_mesh(fileName, component.mesh)

    @property
    def filename(self):
        return self._filename

    def read_file(self):
        return om.read_trimesh(self.filename)

    def translate(self, translate_vector):
        translate_vector = np.array(translate_vector)
        self._translation += translate_vector
        self.regenerate_component()

    def transform_mesh(self, initial_mesh: om.TriMesh):
        fvi = initial_mesh.fv_indices().copy()
        points = initial_mesh.points().copy()
        points += self._translation
        mesh = om.TriMesh(points, fvi)
        return mesh

    def regenerate_component(self):
        if self._original_mesh is not None:
            self.mesh = self.transform_mesh(self._original_mesh)

    def exportGeometry(self, fileName):
        AircraftComponentFromMesh.export_component(fileName, self)

class TiglComponent(AircraftComponent):
    def __init__(self, config,UID):
        name = UID
        super().__init__()
        self._config = config
        self._uid = UID
        self.name = self.UID

    @property
    def config(self) -> CCPACSConfiguration:
        return self._config

    @property
    def uid_mgr(self)->CTiglUIDManager:
        return self.config.get_uidmanager()

    @property
    def UID(self):
        return self._uid

    @property
    def component(self)->ITiglGeometricComponent:
        return self.uid_mgr.get_geometric_component(self.UID)

    @property
    def loft(self)-> CNamedShape:
        return self.component.get_loft()

    @property
    def shape(self) -> TopoDS_Shape:
        return self.loft.shape()

    def regenerate_component(self):
        shape=self.shape
        self.mesh = om.TriMesh()
        if shape is not None:
            self.mesh = get_open_mesh_from_TopoDS_using_shape_tesselator(shape)

class TiglComponentWing(TiglComponent):
    def __init__(self, config,UID):
        super().__init__(config,UID)
        self.regenerate_component()


class TiglComponentFuselage(TiglComponent):
    def __init__(self, config,UID):
        super().__init__(config,UID)
        self.regenerate_component()


class TiglAircraft(TiglComponent):
    def __init__(self, fileName):
        self._filename = fileName
        tigl,tixi,config = self.initialize_tigl()
        self._tixi: Tixi3 = tixi
        self._tigl: Tigl3 = tigl
        super().__init__(config, config.get_uid())
        self._components:List[TiglComponent]=[]
        self.get_tigl_aircraft_components()
        self.regenerate_component()
        self.emit_aircraft_componets_build()


    def __del__(self):
        self._tigl.logToFileDisabled()
        self._tigl.close()
        self._tixi.close()
        self._tixi.cleanup()
    # body of destructor

    def initialize_tigl(self):
        filename = self._filename
        tixi = Tixi3()
        tigl = Tigl3()
        # open cpacs xml with tixi
        try:
            tixi.open(filename)
        except:
            print('Error opening cpacs file with TiXI')
            return None

        # enable logging of errors and warnings into file
        tigl.logSetFileEnding('txt')
        tigl.logSetTimeInFilenameEnabled(TiglBoolean.TIGL_FALSE)
        tigl.logToFileEnabled('adlog')
        # open cpacs with tigl
        try:
            tigl.open(tixi, '')
        except:
            print('Error opening cpacs file with TiGL')
            return None
        # get the configuration manager
        mgr = CCPACSConfigurationManager_get_instance()

        # get the CPACS configuration, defined by the tigl handle
        # we need to access the underlying tigl handle (that is used in the C/C++ API)
        config = mgr.get_configuration(tigl._handle.value)

        return tigl,tixi,config

    def get_tigl_aircraft_components(self):
        # query number of wings and fuselages and their names
        tigl=self._tigl
        nWings = tigl.getWingCount()
        nFuselages = tigl.getFuselageCount()
        for i in range(1, nWings + 1):
            uid = tigl.wingGetUID(i)
            self._components.append(TiglComponentWing(self.config,uid))

        for i in range(1, nFuselages + 1):
            uid = tigl.fuselageGetUID(i)
            self._components.append(TiglComponentFuselage(self.config, uid))



    def prepare_for_regeneration(self):
        pass

    def regenerate_component(self):
        self.sub_geometry.clear()
        for component in self._components:
            self.sub_geometry.append(component)
        pass

    def emit_aircraft_componets_build(self):
        self.emit_geometry_built()
        for ac in self._components:
            ac.emit_geometry_built()
        pass



