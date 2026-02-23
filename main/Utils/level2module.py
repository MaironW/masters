# Define the basis for each Level-2 module class

class Level2Module():
    def __init__(self, name, par):
        self.name = name
        self.par  = par
        self.is_dynamic = False

    def initialize(self, states):
        pass

    def update_algebraic(self, t, parent_states, states, inputs=None):
        pass

    def derivatives(self, t, states):
        return []

    def get_state(self):
        return []

    def set_state(self, vec):
        pass
