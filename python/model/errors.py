class BIGPNError(Exception):
    pass


class InvalidDatasetError(BIGPNError, ValueError):
    pass


class InvalidOptionError(BIGPNError, ValueError):
    pass


class SingularPropagationError(BIGPNError, ArithmeticError):
    pass


class TrainingFailedError(BIGPNError, RuntimeError):
    pass
