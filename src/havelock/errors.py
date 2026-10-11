"""Exception hierarchy (mirrors bolzakit's style; nothing is imported from bolzakit)."""


class HavelockError(Exception):
    """Base class for every error raised by havelock."""


class InvalidInput(HavelockError, ValueError):
    """Unknown model, N < 2, a parameter outside its domain, or a malformed argument.

    Raised before any computation starts."""


class OutOfValidityRange(HavelockError, ValueError):
    """The question is well posed but the answer would be wrong, not merely slow, with the methods available
    (for example the Montaldi–Tokieda family at n = 3, where c4 = 0 and the case rests on MT's own argument)."""


class CertificationFailed(HavelockError):
    """An exact certification could not be completed (precision ladder exhausted, size cap exceeded, or a
    sign pattern that the certificate format cannot express). Raised only when strict=True."""
