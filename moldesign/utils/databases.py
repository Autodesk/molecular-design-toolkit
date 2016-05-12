import os
import dumbdbm
import json
import zlib

class CompressedJsonDbm(object):
    """ Quick-and-dirty interface to a DBM file
    """
    def __init__(self, filename, flag='r', dbm=dumbdbm):
        self.dbm = dbm
        if dbm is dumbdbm:  # dumbdbm corrupts databases even with 'r' flag, so make it read-only
            # TODO: this is bad and non-portable, needs to change before release
            os.system('chmod a-w {0}.dir {0}.dat'.format(filename))

        self.db = self.dbm.open(filename, flag)

    def __getattr__(self, item):
        return getattr(self.db, item)

    def __dir__(self):
        return self.__dict__.keys() + dir(self.db)

    def __len__(self):
        return len(self.db)

    def __getitem__(self, key):
        gzvalue = self.db[key]
        return json.loads(zlib.decompress(gzvalue))

    def __setitem__(self, key, value):
        gzvalue = zlib.compress(json.dumps(value))
        self.db[key] = gzvalue


class Hdf5Dbm(object):
    """ A simplified DBM interface backed by HDF5.

    Like DBMs, it only stores strings. Unlike DBMs, the strings must be completely ascii

    Note:
        if something other than a string is passed, it will be converted using json.dumps

    This was written only because dumbdbm has a tendency to get corrupted (because it doesn't
    support read-only mode), none of the other standard library DBMs are portable, and other
    alternatives (like semidbm) play nice with git.
    """
    def __init__(self, filename, mode='r', compression=10):
        """Initialization:

        Args:
            filename (str): path to file (including suffix)
            mode (str): ``mode`` argument for ``h5py.File``
            compression (str or int): compression for datasets (if int: gzip level)
        """
        import h5py

        self.f = h5py.File(filename, mode)
        self._strtype = h5py.special_dtype(vlen=bytes)
        self.compression = compression

    @classmethod
    def open(cls, *args, **kwargs):
        """Synonym for __init__ for compatibility with stdlib DBM modules"""
        return cls(*args, **kwargs)

    def __getitem__(self, item):
        ds = self.f[item]
        return json.loads(ds.value)

    def __setitem__(self, item, value):
        jsval = json.dumps(value)
        if item in self.f:
            del self.f[item]
        self.f.create_dataset(item, data=jsval,
                              dtype=self._strtype, compression=self.compression)

    def close(self):
        return self.f.close()






