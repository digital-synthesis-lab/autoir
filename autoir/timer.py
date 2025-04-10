from functools import wraps
from time import perf_counter, time


class Timer:
    def __enter__(self):
        self.start = perf_counter()
        return self

    def __exit__(self, type, value, traceback):
        self.time = perf_counter() - self.start

    def __str__(self):
        return f"Time: {self.time:.3f} seconds"
