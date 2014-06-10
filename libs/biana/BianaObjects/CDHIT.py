# CD HIT related classes

class CDHITCluster(object):

    def __init__(self, representant = None, matches = None):

        self.representant = representant
        
        if matches is None:
            self.matches = []
        else:
            self.matches = matches

    def set_representant(self, representant):
        self.representant = representant

    def get_representant(self):
        return self.representant

    def add_match(self, match):
        self.matches.append(match)
        
    def get_matches(self):
        return self.matches


class CDHITMatch(object):

    def __init__(self, sequenceID, start_pos, end_pos, repr_start_pos, repr_end_pos, identity):

        self.sequenceID = sequenceID
        self.start_pos = start_pos
        self.end_pos = end_pos
        self.repr_start_pos = repr_start_pos
        self.repr_end_pos = repr_end_pos
        self.identity = identity


