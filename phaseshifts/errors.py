"""Module defining common exception classes."""

class CoordinatesError(Exception):
    """Coordinate exception to raise and log duplicate coordinates."""

    def __init__(self, msg, *args, **kwargs):
        super(CoordinatesError, self).__init__(msg, *args, **kwargs)
        self.msg = "CoordinatesError: %s" % msg

    def __str__(self):
        return self.msg

    def __unicode__(self):
        return self.msg
