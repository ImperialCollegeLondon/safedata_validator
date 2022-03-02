# A type for representing NA - this isn't meant to be a true singleton,
# just a consistent value used across modules and which works in isinstance

class NA:
    def __init__(self):
        pass

NA_type = type(NA())
