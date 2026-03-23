"""Unit tests for the coevo CLI (coevo/cli.py)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from typer.testing import CliRunner

from coevo.cli import app

runner = CliRunner()

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _write_taxids(path: Path, taxids: set[int]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(str(t) for t in sorted(taxids)) + "\n")


def _write_tsv(path: Path, content: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(content)


# ---------------------------------------------------------------------------
# Tests for draw-tree command
# ---------------------------------------------------------------------------


class TestDrawTree:
    """Tests for the ``coevo draw-tree`` CLI command."""

    def _minimal_config(self, tmp_path: Path) -> tuple[str, Path]:
        """Write a minimal config.yaml and populate taxa files; return (config_path, results_dir)."""
        results_dir = tmp_path / "results"
        config_path = tmp_path / "config.yaml"
        config_path.write_text(
            f"output:\n  results_dir: {results_dir}\n"
        )
        taxa_dir = results_dir / "taxa"
        _write_taxids(taxa_dir / "protein_taxa.txt", {1, 2, 3})
        _write_taxids(taxa_dir / "rna_taxa.txt", {3, 4, 5})
        return str(config_path), results_dir

    # ------------------------------------------------------------------
    # Help text
    # ------------------------------------------------------------------

    def test_help(self) -> None:
        result = runner.invoke(app, ["draw-tree", "--help"])
        assert result.exit_code == 0
        assert "draw-tree" in result.output.lower() or "circular" in result.output.lower()

    # ------------------------------------------------------------------
    # Successful run – output file
    # ------------------------------------------------------------------

    def test_draw_tree_with_output(self, tmp_path: Path) -> None:
        config_path, results_dir = self._minimal_config(tmp_path)
        output_file = str(tmp_path / "tree.png")

        mock_draw = MagicMock(return_value=MagicMock())
        with patch("coevo.analysis.phylogeny.draw_circular_tree", mock_draw):
            result = runner.invoke(
                app,
                [
                    "draw-tree",
                    "--config", config_path,
                    "--output", output_file,
                ],
            )

        assert result.exit_code == 0, result.output
        mock_draw.assert_called_once()
        call_kwargs = mock_draw.call_args[1]
        assert call_kwargs["output_file"] == output_file
        assert call_kwargs["protein_taxids"] == {1, 2, 3}
        assert call_kwargs["rna_taxids"] == {3, 4, 5}

    # ------------------------------------------------------------------
    # Successful run – interactive (no output file)
    # ------------------------------------------------------------------

    def test_draw_tree_interactive(self, tmp_path: Path) -> None:
        config_path, _ = self._minimal_config(tmp_path)

        mock_draw = MagicMock(return_value=MagicMock())
        with patch("coevo.analysis.phylogeny.draw_circular_tree", mock_draw):
            result = runner.invoke(
                app,
                ["draw-tree", "--config", config_path],
            )

        assert result.exit_code == 0, result.output
        mock_draw.assert_called_once()
        call_kwargs = mock_draw.call_args[1]
        assert call_kwargs["output_file"] is None

    # ------------------------------------------------------------------
    # Custom max_nodes and collapse_threshold
    # ------------------------------------------------------------------

    def test_custom_options(self, tmp_path: Path) -> None:
        config_path, _ = self._minimal_config(tmp_path)

        mock_draw = MagicMock(return_value=MagicMock())
        with patch("coevo.analysis.phylogeny.draw_circular_tree", mock_draw):
            result = runner.invoke(
                app,
                [
                    "draw-tree",
                    "--config", config_path,
                    "--max-nodes", "50",
                    "--collapse-threshold", "0.75",
                ],
            )

        assert result.exit_code == 0, result.output
        call_kwargs = mock_draw.call_args[1]
        assert call_kwargs["max_nodes"] == 50
        assert call_kwargs["collapse_threshold"] == pytest.approx(0.75)

    # ------------------------------------------------------------------
    # Motif highlights loaded when auxiliary files are present
    # ------------------------------------------------------------------

    def test_motif_taxids_loaded_when_files_present(self, tmp_path: Path) -> None:
        config_path, results_dir = self._minimal_config(tmp_path)

        # seq_index 1 and 2 have motif; map seq_index -> taxid via metadata
        _write_tsv(
            results_dir / "analysis" / "motif_results.tsv",
            "sequence_id\tmotif_present\nseq_1\tTrue\nseq_2\tFalse\n",
        )
        _write_tsv(
            results_dir / "alignment" / "16s_metadata.tsv",
            "seq_index\ttaxid\n1\t10\n2\t20\n",
        )

        mock_draw = MagicMock(return_value=MagicMock())
        with patch("coevo.analysis.phylogeny.draw_circular_tree", mock_draw):
            result = runner.invoke(
                app,
                ["draw-tree", "--config", config_path],
            )

        assert result.exit_code == 0, result.output
        call_kwargs = mock_draw.call_args[1]
        # seq_index 1 has motif_present=True, taxid=10
        assert call_kwargs["motif_taxids"] == {10}

    # ------------------------------------------------------------------
    # Motif files absent → motif_taxids is None
    # ------------------------------------------------------------------

    def test_no_motif_highlights_when_files_absent(self, tmp_path: Path) -> None:
        config_path, _ = self._minimal_config(tmp_path)

        mock_draw = MagicMock(return_value=MagicMock())
        with patch("coevo.analysis.phylogeny.draw_circular_tree", mock_draw):
            result = runner.invoke(
                app,
                ["draw-tree", "--config", config_path],
            )

        assert result.exit_code == 0, result.output
        call_kwargs = mock_draw.call_args[1]
        assert call_kwargs["motif_taxids"] is None

    # ------------------------------------------------------------------
    # Missing taxa file → non-zero exit
    # ------------------------------------------------------------------

    def test_missing_taxa_file_raises(self, tmp_path: Path) -> None:
        results_dir = tmp_path / "results"
        config_path = tmp_path / "config.yaml"
        config_path.write_text(f"output:\n  results_dir: {results_dir}\n")
        # Intentionally do NOT create the taxa files

        result = runner.invoke(
            app,
            ["draw-tree", "--config", str(config_path)],
        )

        assert result.exit_code != 0

    # ------------------------------------------------------------------
    # draw-tree is listed in the top-level help
    # ------------------------------------------------------------------

    def test_draw_tree_in_app_help(self) -> None:
        result = runner.invoke(app, ["--help"])
        assert result.exit_code == 0
        assert "draw-tree" in result.output


# ---------------------------------------------------------------------------
# Tests for analyse command
# ---------------------------------------------------------------------------


class TestAnalyse:
    """Tests for the ``coevo analyse`` CLI command."""

    def _minimal_config(self, tmp_path: Path) -> tuple[str, Path]:
        results_dir = tmp_path / "results"
        config_path = tmp_path / "config.yaml"
        config_path.write_text(f"output:\n  results_dir: {results_dir}\n")
        taxa_dir = results_dir / "taxa"
        _write_taxids(taxa_dir / "protein_taxa.txt", {1, 2})
        _write_taxids(taxa_dir / "rna_taxa.txt", {2, 3})
        return str(config_path), results_dir

    def test_help(self) -> None:
        result = runner.invoke(app, ["analyse", "--help"])
        assert result.exit_code == 0

    def test_analyse_writes_output_files(self, tmp_path: Path) -> None:
        """analyse command must write enterobacteriaceae_summary.tsv and motif_histogram.tsv."""
        from unittest.mock import patch, MagicMock
        config_path, results_dir = self._minimal_config(tmp_path)

        with patch("coevo.analysis.phylogeny.HAS_ETE4", False):
            result = runner.invoke(app, ["analyse", "--config", config_path])

        assert result.exit_code == 0, result.output
        assert (results_dir / "analysis" / "enterobacteriaceae_summary.tsv").exists()
        assert (results_dir / "analysis" / "motif_histogram.tsv").exists()

    def test_analyse_enterobacteriaceae_summary_has_columns(self, tmp_path: Path) -> None:
        """enterobacteriaceae_summary.tsv must contain the expected column headers."""
        config_path, results_dir = self._minimal_config(tmp_path)

        with patch("coevo.analysis.phylogeny.HAS_ETE4", False):
            result = runner.invoke(app, ["analyse", "--config", config_path])

        assert result.exit_code == 0, result.output
        content = (results_dir / "analysis" / "enterobacteriaceae_summary.tsv").read_text()
        assert "rank" in content
        assert "name" in content

    def test_analyse_histogram_when_motif_file_present(self, tmp_path: Path) -> None:
        """When motif_results.tsv exists with offset data, motif_histogram.tsv is populated."""
        config_path, results_dir = self._minimal_config(tmp_path)
        # Write a motif_results.tsv with offset data
        _write_tsv(
            results_dir / "analysis" / "motif_results.tsv",
            "sequence_id\tmotif_present\tmotif_offset\nseq_1\tTrue\t0\nseq_2\tTrue\t1\nseq_3\tFalse\t\n",
        )

        with patch("coevo.analysis.phylogeny.HAS_ETE4", False):
            result = runner.invoke(app, ["analyse", "--config", config_path])

        assert result.exit_code == 0, result.output
        content = (results_dir / "analysis" / "motif_histogram.tsv").read_text()
        # Histogram should have rows for offsets 0 and 1
        assert "0" in content
        assert "1" in content
