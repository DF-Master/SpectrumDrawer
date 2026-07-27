"""Configuration manager: load YAML config, merge with CLI overrides."""

import os
import yaml
from typing import Any, Dict, Optional


class ConfigManager:
    """Load and manage configuration from YAML."""

    DEFAULT_CONFIG_DIR = os.path.join(os.path.dirname(__file__))
    DEFAULT_CONFIG_PATH = os.path.join(DEFAULT_CONFIG_DIR, 'default_config.yaml')

    def __init__(self, config_path: Optional[str] = None):
        """Load configuration from file.

        Parameters
        ----------
        config_path : str or None
            Path to a custom YAML config file. If None, uses default.
        """
        self._data = self._load_default()
        if config_path and os.path.isfile(config_path):
            custom = self._load_yaml(config_path)
            self._deep_update(self._data, custom)

    def _load_default(self) -> Dict[str, Any]:
        return self._load_yaml(self.DEFAULT_CONFIG_PATH)

    @staticmethod
    def _load_yaml(path: str) -> Dict[str, Any]:
        with open(path, 'r', encoding='utf-8') as f:
            return yaml.safe_load(f) or {}

    @staticmethod
    def _deep_update(base: dict, update: dict):
        """Recursively update base dict with values from update."""
        for key, value in update.items():
            if isinstance(value, dict) and isinstance(base.get(key), dict):
                ConfigManager._deep_update(base[key], value)
            else:
                base[key] = value

    def get(self, *keys: str, default: Any = None) -> Any:
        """Traverse nested config: config.get('colors', 'b_ion')."""
        node = self._data
        for key in keys:
            if isinstance(node, dict):
                node = node.get(key)
                if node is None:
                    return default
            else:
                return default
        return node

    def apply_cli_overrides(self, **overrides):
        """Apply individual CLI overrides to the config.

        Keys use dotted notation: 'processing.tol_ppm'.
        """
        for dotted_key, value in overrides.items():
            if value is None:
                continue
            parts = dotted_key.split('.')
            node = self._data
            for part in parts[:-1]:
                if part not in node:
                    node[part] = {}
                node = node[part]
            node[parts[-1]] = value

    @property
    def tol_ppm(self) -> float:
        return self.get('processing', 'tol_ppm', default=20.0)

    @property
    def ion_types(self) -> str:
        return self.get('processing', 'ion_types', default='by')

    @property
    def max_charge(self) -> int:
        return self.get('processing', 'max_charge', default=2)

    @property
    def colors(self) -> Dict[str, str]:
        return self.get('colors', default={})

    @property
    def figure_config(self) -> Dict[str, Any]:
        return self.get('figure', default={})

    @property
    def ladder_config(self) -> Dict[str, Any]:
        return self.get('ladder', default={})

    @property
    def spectrum_config(self) -> Dict[str, Any]:
        return self.get('spectrum', default={})

    @property
    def mass_error_config(self) -> Dict[str, Any]:
        return self.get('mass_error', default={})

    @property
    def fix_mod_names(self) -> list:
        return self.get('modifications', 'fixed', default=['Carbamidomethyl[C]'])

    @property
    def var_mod_names(self) -> list:
        """Variable modification names from config (reserved interface)."""
        return self.get('modifications', 'variable', default=[])
