"""
Functions to work with the IUPAC glycan nomenclature.
"""
import re
from typing import Any


def reformat_link(string):
    """
    Reformat a linkage string from the IUPAC format to the biobuild (CHARMM force field) format.

    Parameters
    ----------
    string : str
        The linkage string in IUPAC format.

    Returns
    -------
    str
        The linkage string in biobuild format.
    """

    type, a, b = re.match("([ab])(\d+)-(\d+)", string).groups()
    return a + b + type * 2


class IUPACParser:
    """
    A parser for condensed IUPAC glycan nomenclature strings. This class will generate
    a list of connecting glycan segments from a string from which a Molecule can be built.
    """

    def __init__(self):
        self._string = ""
        self.reset()

    @property
    def _current(self):
        return self._string[-self._idx]

    @property
    def _next(self):
        return self._string[-self._idx - 1]

    @property
    def _can_shift(self):
        return self._idx < len(self._string)

    @property
    def _is_at_bracket(self):
        return self._current in ("(", ")")

    @property
    def _is_at_open_square_bracket(self):
        return self._current == "]"

    @property
    def _is_at_close_square_bracket(self):
        return self._current == "["

    @property
    def _can_store(self):
        return (
            not self._is_still_parsing
            and len(self._latest_residue) >= 1
            and len(self._latest_linkage) >= 3
            and len(self._second_latest_residue) >= 1
        )

    @property
    def _is_still_parsing(self):
        if not self._can_shift:
            return False
        return self._current not in ("(", ")", "[", "]")

    def parse(self, string):
        """
        Parse a string of IUPAC glycan nomenclature into a list of glycan segments.

        Parameters
        ----------
        string : str
            The IUPAC glycan nomenclature string.

        Returns
        -------
        list
            A list of tuples where each segment is a tuple of (residue1, residue2, linkage).
        """
        self._string = self._prep_greek_letters(string)
        self.reset()
        self._parse()
        return self._glycan

    def reset(self):
        """
        Reset the parser.
        """
        self._glycan = []
        self._idx = 1
        self._residue_counts = {}
        self._latest_residue = ""
        self._latest_conformation = ""
        self._latest_linkage = ""
        self._second_latest_residue = ""
        self._second_latest_conformation = ""
        self._latest_residue_before_square_bracket = ""
        self._latest_conformation_before_square_bracket = ""

    def _shift(self):
        self._idx += 1

    def _shift_residue(self):
        self._second_latest_residue = self._latest_residue
        self._latest_residue = ""
        self._second_latest_conformation = self._latest_conformation
        self._latest_conformation = ""

    def _parse(self):
        self._crop_end()
        while self._can_shift:
            if self._can_store:
                self._store()
                continue
            if self._is_at_bracket:
                self._shift_residue()
                self._latest_linkage = self._parse_linkage()
                self._shift()
                continue
            if self._is_at_open_square_bracket:
                self._latest_residue_before_square_bracket = self._fit_residue(
                    self._latest_residue, increase_count=False
                )
                self._latest_residue = self._latest_residue_before_square_bracket
                self._latest_conformation_before_square_bracket = (
                    self._latest_conformation
                )
                self._shift()
                continue
            if self._is_at_close_square_bracket:
                self._latest_residue = self._latest_residue_before_square_bracket
                self._latest_conformation = (
                    self._latest_conformation_before_square_bracket
                )
                self._second_latest_residue = ""
                self._second_latest_conformation = ""
                self._shift()
                continue
            self._latest_residue += self._current
            self._shift()
        self._latest_residue += self._current
        self._store()

    def _store(self):
        second = self._second_latest_residue
        if "@" not in second:
            second = self._fit_residue(second)
            self._second_latest_residue = second

        latest = self._latest_residue
        if not "@" in latest:
            latest = self._fit_residue(latest)
            self._latest_residue = latest

        branch = (second, latest, self._reformat_link(self._latest_linkage))
        self._glycan.append(branch)
        self._latest_linkage = ""

    def _fit_residue(self, r, increase_count=True):
        if "@" in r:
            return r

        r = r[::-1]
        if r in self._residue_counts:
            if increase_count:
                self._residue_counts[r] += 1
            r += "@" + str(self._residue_counts[r])
        else:
            self._residue_counts[r] = 1
            r += "@1"
        return r

    def _parse_linkage(self):
        self._shift()
        linkage = ""
        while self._can_shift and not self._is_at_bracket:
            linkage += self._current
            self._shift()
        self._latest_conformation = linkage[-1]
        linkage = linkage[:-1]
        return linkage

    def _crop_end(self):
        if self._string[-1] == "-":
            while self._next != "(":
                self._shift()
            self._latest_conformation = self._current
            self._shift()
            self._string = self._string[: -self._idx]
            self._idx = 1

    def _reformat_link(self, link):
        link = link[::-1].replace("-", "")
        link = link + self._second_latest_conformation + self._latest_conformation
        return link

    def _prep_greek_letters(self, string):
        string = string.replace("α", "a").replace("β", "b")
        return string

    def __call__(self, *args: Any, **kwds: Any) -> Any:
        return self.parse(*args, **kwds)


if __name__ == "__main__":
    string2 = "Man(b1-6)[Man(b1-3)]BMA(b1-4)GlcNAc(b1-4)GlcNAc(b1-"
    string = "F(b1-4)[E(a2-3)D(a1-4)]C(a1-6)B(b1-4)A(a1-"

    p = IUPACParser()
    g = p.parse(string2)
    print(g)
