"""Logging functions of the package."""
# src/jazzy/logging.py
import logging


logger = logging.getLogger(__name__)
logger.setLevel(logging.WARNING)
formatter = logging.Formatter("Jazzy %(levelname)s: [%(asctime)s] %(message)s")
formatter.datefmt = "%H:%M:%S"
handler = logging.StreamHandler()
handler.setFormatter(formatter)
logger.addHandler(handler)
