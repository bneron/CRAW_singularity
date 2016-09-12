import shutil
import tempfile
import os

from tests import Gen2DockTest
from gensoft2docker.package import Package, Reference
from gensoft2docker.config import Config
from gensoft2docker.error import Gensoft2DockerError


class Test(Gen2DockTest):

    def setUp(self):
        if 'GEN2DOCK_HOME' in os.environ:
            self.gen2dock_home = os.environ['GEN2DOCK_HOME']
            self.local_install = True
        else:
            self.local_install = False
            self.gen2dock_home = os.path.normpath(os.path.abspath(os.path.join(os.path.dirname(__file__), '..' '..')))
        self.tmp_dir = tempfile.gettempdir()
        self.out_dir = os.path.join(self.tmp_dir, 'Gensof2Docker_test')
        os.makedirs(self.out_dir)
        self.cfg = Config()
        self._set_config()

    def _set_config(self):
        self.cfg.set('server', 'http://my-domain.ext')
        self.cfg.set('pack_sch', 'package_schema.yml')
        self.cfg.set('pack_vers_sch', 'package_version_schema.yml')
        self.cfg.set('docker_tpl', 'Dockerfile.j2')
        self.cfg.set('build_tpl', 'buildDeps.j2')
        self.cfg.set('log_level', 'WARNING')
        self.cfg.set('inst_dir', os.path.join(self._data_dir, 'inst'))
        self.cfg.set('schemas_dir', os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', 'schemas')))
        self.cfg.set('templates_dir', os.path.normpath(os.path.join(os.path.dirname(__file__), '..', '..', 'templates')))
        self.cfg.set('dest_dir', self.out_dir)


    def tearDown(self):
        try:
            shutil.rmtree(self.out_dir)
        except:
            pass


    def test_title(self):
        ref = Reference()
        self.assertIsNone(ref.title)
        title = "this is a title"
        ref = Reference(title=title)
        self.assertEqual(ref.title, title)

    def test_idType(self):
        ref = Reference()
        self.assertIsNone(ref.idType)
        idType = "doi"
        ref = Reference(idType=idType)
        self.assertEqual(ref.idType, idType)

    def test_ID(self):
        ref = Reference()
        self.assertIsNone(ref.ID)
        ID = "12345"
        ref = Reference(ID=ID)
        self.assertEqual(ref.ID, ID)

    def test_get_url(self):
        ref = Reference(title='title',
                        idType='doi',
                        ID='12345')
        self.assertEqual(ref.get_url(), 'http://wwww.doi.org/12345')
        ref = Reference(title='title',
                        idType='pmid',
                        ID='12345')
        self.assertEqual(ref.get_url(), 'http://www.ncbi.nlm.nih.gov/pubmed/12345')
        ref = Reference(title='title',
                        idType='pmicid',
                        ID='12345')
        self.assertEqual(ref.get_url(), 'http://www.ncbi.nlm.nih.gov/pmc/articles/12345')
        ref = Reference(title='title',
                        idType='whatever',
                        ID='12345')
        with self.assertRaises(Gensoft2DockerError) as ctx:
            ref.get_url()
        self.assertEqual(ctx.exception.message, "idType: 'whatever' NOT managed")

    def test_str(self):
        ref = Reference(title='title',
                        idType='doi',
                        ID='12345')
        self.assertEqual(str(ref), 'title http://wwww.doi.org/12345')

        ref = Reference(title='title',
                        ID='12345')
        self.assertEqual(str(ref), 'title')
