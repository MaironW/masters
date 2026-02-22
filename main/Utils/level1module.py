# Define the basis for each Level-1 module class

class Level1Module():
    def __init__(self, par_override=None):
        self.modules = []

    # Organize dynamic modules
    @property
    def dynamic_modules(self):
        return [m for m in self.modules if m.is_dynamic]

    # Update time-dependent, non-integrated Level-2 modules
    def update_algebraic(self, t, parent_states, inputs):
        for m in self.modules:
            m.update_algebraic(t, parent_states.snapshot(), self.snapshot(), inputs)

    # Return a dict of modules that have dynamic (integrated) states
    def get_dynamic_modules(self):
        return self.dynamic_modules

    def get_state(self):
        return [m.get_state() for m in self.dynamic_modules]

    def set_state(self, vecs):
        for m, v in zip(self.dynamic_modules, vecs):
            m.set_state(v)

    def derivatives(self, t):
        return [m.derivatives(t, self.snapshot()) for m in self.dynamic_modules]

    def snapshot(self):
        return {m.name: m.state.copy() for m in self.modules}