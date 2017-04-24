try:
    from tests import CRAWTest
except ImportError as err:
    msg = "Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err)
    raise ImportError("Cannot import craw, check your installation or your CRAW_HOME variable : {0!s}".format(err))

from craw.util import progress


class TestUtil(CRAWTest):


    def test_progress(self):
        import io
        stream = io.StringIO()
        progress(10, 100, status="", file=stream)
        self.assertEqual('[======>                                                     ] 10.0% ... of 100 annotations\r',
                         stream.getvalue())

