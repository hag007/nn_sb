import multiprocessing
from multiprocessing.pool import Pool

def func_star(a_b):
    """Convert `f([1,2])` to `f(1,2)` call."""
    return a_b[0](*a_b[1])


class NoDaemonProcess(multiprocessing.Process):
    # make 'daemon' attribute always return False
    def _get_daemon(self):
        return False

    def _set_daemon(self, value):
        pass

    daemon = property(_get_daemon, _set_daemon)


# We sub-class multiprocessing.pool.Pool instead of multiprocessing.Pool
# because the latter is only a wrapper function, not a proper class.
class MyPool(Pool):
    Process = NoDaemonProcess


