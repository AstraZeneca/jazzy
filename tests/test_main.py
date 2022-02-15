"""Test cases for the __main__ module."""
import pytest
from click.testing import CliRunner

from jazzy import __main__


@pytest.fixture
def runner() -> CliRunner:
    """Fixture for invoking command-line interfaces."""
    return CliRunner()


def test_cli_succeeds(runner: CliRunner) -> None:
    """It exits with a status code of zero."""
    result = runner.invoke(__main__.cli)
    assert result.exit_code == 0


def test_vec_succeeds(runner: CliRunner) -> None:
    """It exits with a status code of zero."""
    smiles = "CC"
    result = runner.invoke(__main__.vec, [smiles])
    assert result.exit_code == 0


def test_vec_fails_without_smiles(runner: CliRunner) -> None:
    """It exits with no SMILES given."""
    smiles = None
    result = runner.invoke(__main__.vec, [smiles])
    assert result.exit_code == 1


def test_vec_fails_with_invalid_minimisation_method(runner: CliRunner) -> None:
    """It fails with an invalid minimisation method."""
    smiles = "CC"
    result = runner.invoke(__main__.vec, ["--opt", "invalid", smiles])
    assert result.exit_code == 1


def test_vec_fails_with_invalid_smiles(runner: CliRunner) -> None:
    """It fails with an invalid smiles."""
    smiles = "xxx"
    result = runner.invoke(__main__.vec, [smiles])
    assert "{'__status': 'error', 'smiles': 'xxx'}" in result.output


def test_vis_succeeds(runner: CliRunner) -> None:
    """It exits with a status code of zero."""
    smiles = "CC"
    result = runner.invoke(__main__.vis, ["--opt", "MMFF94", smiles])
    assert result.exit_code == 0


def test_vis_fails_with_invalid_smiles(runner: CliRunner) -> None:
    """It fails with an invalid smiles."""
    smiles = "xxx"
    result = runner.invoke(__main__.vis, ["--opt", "MMFF94", smiles])
    assert "{'__status': 'error', 'smiles': 'xxx'}" in result.output


def test_vis_fails_with_invalid_minimisation_method(runner: CliRunner) -> None:
    """It fails with an invalid minimisation method."""
    smiles = "CC"
    result = runner.invoke(__main__.vis, ["--opt", "invalid", smiles])
    assert result.exit_code == 1
