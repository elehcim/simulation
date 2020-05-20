from pynbody import units
from pynbody.filt import Filter


class Slit(Filter):
    """Emulate a Slit, so taking a slice at z=const of a snap"""
    def __init__(self, x1, y1=None, x2=None, y2=None):
        super().__init__()
        self._descriptor = "slit"
        x1, y1, x2, y2 = [units.Unit(x) if isinstance(x, str) else x for x in (x1, y1, x2, y2)]
        if y1 is None:
            y1 = x1
        if x2 is None:
            x2 = -x1
        if y2 is None:
            y2 = -y1

        self.x1, self.y1, self.x2, self.y2 = x1, y1, x2, y2

    def __call__(self, sim):
        x1, y1, x2, y2 = [x.in_units(sim["pos"].units, **sim["pos"].conversion_context())
                                  if units.is_unit_like(x) else x
                                  for x in (self.x1, self.y1, self.x2, self.y2)]

        return (sim["x"] > x1) * (sim["x"] < x2) * (sim["y"] > y1) * (sim["y"] < y2)

    def __repr__(self):
        x1, y1, x2, y2 = ["'%s'" % str(x)
                                  if units.is_unit_like(x) else x
                                  for x in (self.x1, self.y1, self.x2, self.y2)]
        return "Slit(%s, %s, %s, %s)" % (x1, y1, x2, y2)

