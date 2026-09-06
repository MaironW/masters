import time
from collections import defaultdict

class Profiler:
    def __init__(self):
        self.times = defaultdict(float)
        self.calls = defaultdict(int)

    def start(self):
        return time.perf_counter()

    def stop(self, name, t0):
        self.times[name] += time.perf_counter() - t0
        self.calls[name] += 1

    def report(self):
        total = sum(self.times.values())

        print("\n" + "=" * 70)
        print("RUNTIME PROFILER")
        print("=" * 70)

        print(f"{'Function':35s} {'Time [s]':>12s} {'%':>8s} {'Calls':>12s}")
        print("-" * 70)

        for name, elapsed in sorted(
            self.times.items(),
            key=lambda x: x[1],
            reverse=True
        ):
            percentage = 100.0 * elapsed / total if total > 0 else 0.0

            print(
                f"{name:35s} "
                f"{elapsed:12.3f} "
                f"{percentage:7.2f}% "
                f"{self.calls[name]:12d}"
            )

        print("-" * 70)
        print(f"{'Measured time':35s} {total:12.3f} s")
        print("=" * 70)