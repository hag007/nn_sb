import datetime


class Stopwatch(object):
    """A simple timer class"""

    def __init__(self):
        pass

    def start(self):
        """Starts the timer"""
        self._start = datetime.datetime.now()
        return self._start

    def stop(self, message="Total: "):
        """Stops the timer.  Returns the time elapsed"""
        self._stop = datetime.datetime.now()
        return message + str(self._stop - self._start)

    def now(self, message="Now: "):
        """Returns the current time with a message"""
        return message + ": " + str(datetime.datetime.now())

    def elapsed(self, message="Elapsed: "):
        """Time elapsed since start was called"""
        return message + str(datetime.datetime.now() - self._start)

    def split(self, message="Split started at: "):
        """Start a split timer"""
        self._split_start = datetime.datetime.now()
        return message + str(self._split_start)

    def unsplit(self, message="Unsplit: "):
        """Stops a split. Returns the time elapsed since split was called"""
        return message + str(datetime.datetime.now() - self._split_start)