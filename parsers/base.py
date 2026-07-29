"""Abstract base class for identification result parsers."""

from abc import ABC, abstractmethod
from typing import List
from ..models import Identification


class BaseIdentificationParser(ABC):
    """Abstract base for parsing search engine identification results."""

    @abstractmethod
    def parse(self, path: str) -> List[Identification]:
        """Parse an identification result file.

        Returns
        -------
        list of Identification
            Standardized identification entries.
        """
        ...

    @staticmethod
    def get_parser(parser_name: str) -> 'BaseIdentificationParser':
        """Factory: return appropriate parser by name."""
        from .psimxl_parser import PsimxlParser
        from .plink_parser import PlinkParser
        from .pfind_parser import PfindParser
        parsers = {
            'psimxl': PsimxlParser,
            'simxl': PsimxlParser,
            'plink': PlinkParser,
            'pfind': PfindParser,
        }
        if parser_name.lower() in parsers:
            return parsers[parser_name.lower()]()
        else:
            raise ValueError(
                f"Unknown parser: '{parser_name}'. "
                f"Available: {list(parsers.keys())}"
            )
