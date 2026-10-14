from __future__ import annotations

import json
from types import SimpleNamespace

import pandas as pd
import pytest

from cellpilot.interpretation import (
    ClusterInterpretation,
    InterpretationBundle,
    build_evidence_packet,
    deterministic_interpret,
    openai_interpret,
    validate_bundle,
)


def _packet() -> dict:
    markers = pd.DataFrame(
        {
            "cluster": ["0", "0", "0", "1", "1", "1"],
            "rank": [1, 2, 3, 1, 2, 3],
            "names": ["MS4A1", "CD79A", "CD74", "NKG7", "GNLY", "PRF1"],
            "scores": [10, 9, 8, 10, 9, 8],
            "logfoldchanges": [3, 2, 2, 3, 2, 2],
            "pvals_adj": [1e-8] * 6,
        }
    )
    counts = pd.DataFrame({"cluster": ["0", "1"], "n_cells": [40, 60]})
    registry = {
        "source": {"title": "test registry"},
        "panels": {
            "B cell": {"markers": ["MS4A1", "CD79A", "CD74"]},
            "NK cell": {"markers": ["NKG7", "GNLY", "PRF1"]},
        },
    }
    return build_evidence_packet(markers, counts, marker_registry=registry)


def test_deterministic_interpretation_is_grounded() -> None:
    packet = _packet()
    bundle = deterministic_interpret(packet)
    validate_bundle(bundle, packet)
    labels = {item.cluster: item.proposed_label for item in bundle.interpretations}
    assert labels == {"0": "B cell", "1": "NK cell"}
    assert all(item.evidence_ids for item in bundle.interpretations)


def test_cluster_counts_survive_pandas_numeric_row_coercion() -> None:
    markers = pd.DataFrame(
        {
            "cluster": [0, 0, 1, 1],
            "rank": [1, 2, 1, 2],
            "names": ["MS4A1", "CD79A", "NKG7", "GNLY"],
        }
    )
    # The float fraction column makes iterrows coerce integral cluster IDs to
    # floats in affected pandas versions (0 -> 0.0).
    counts = pd.DataFrame(
        {"cluster": [0, 1], "n_cells": [40, 60], "fraction": [0.4, 0.6]}
    )
    registry = {
        "source": {"title": "test", "url": "https://example.org/evidence"},
        "panels": {
            "B cell": {"markers": ["MS4A1", "CD79A"]},
            "NK cell": {"markers": ["NKG7", "GNLY"]},
        },
    }
    packet = build_evidence_packet(markers, counts, marker_registry=registry)

    assert [(item["cluster"], item["n_cells"]) for item in packet["clusters"]] == [
        ("0", 40),
        ("1", 60),
    ]


def test_unknown_evidence_id_is_rejected() -> None:
    packet = _packet()
    bundle = InterpretationBundle(
        interpretations=[
            ClusterInterpretation(
                cluster="0",
                proposed_label="B cell",
                confidence="low",
                evidence_ids=["invented:evidence"],
                uncertainties=["test"],
            ),
            ClusterInterpretation(
                cluster="1",
                proposed_label="NK cell",
                confidence="low",
                evidence_ids=["cluster:1:marker:NKG7"],
                uncertainties=["test"],
            ),
        ]
    )
    with pytest.raises(ValueError, match="Unknown evidence IDs"):
        validate_bundle(bundle, packet)


def test_openai_backend_without_key_makes_no_client_or_request(monkeypatch) -> None:
    """An explicitly requested OpenAI run must fail closed before client creation."""

    import openai

    constructed = False

    class ForbiddenClient:
        def __init__(self, **_kwargs) -> None:
            nonlocal constructed
            constructed = True
            raise AssertionError("OpenAI client must not be constructed without a key")

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setattr(openai, "OpenAI", ForbiddenClient)

    with pytest.raises(
        RuntimeError,
        match=r"OPENAI_API_KEY is not set.*No external API request was made",
    ):
        openai_interpret(_packet(), model="test-model")

    assert constructed is False


def test_openai_backend_honors_custom_api_key_environment_name(monkeypatch) -> None:
    import openai

    captured: dict = {}
    expected = deterministic_interpret(_packet())

    class FakeResponses:
        def parse(self, **_kwargs):
            return SimpleNamespace(output_parsed=expected, id="resp_custom_env")

    class FakeClient:
        def __init__(self, *, api_key: str) -> None:
            captured["api_key"] = api_key
            self.responses = FakeResponses()

    monkeypatch.delenv("OPENAI_API_KEY", raising=False)
    monkeypatch.setenv("CELLPILOT_TEST_OPENAI_KEY", "custom-environment-key")
    monkeypatch.setattr(openai, "OpenAI", FakeClient)

    _, metadata = openai_interpret(
        _packet(),
        model="test-model",
        api_key_env="CELLPILOT_TEST_OPENAI_KEY",
    )

    assert captured["api_key"] == "custom-environment-key"
    assert metadata["external_api_called"] is True


def test_mocked_openai_structured_output_is_locally_validated(monkeypatch) -> None:
    """Exercise the Responses parse contract without making a network request."""

    import openai

    packet = _packet()
    expected = deterministic_interpret(packet)
    captured: dict = {}

    class FakeResponses:
        def parse(self, **kwargs):
            captured.update(kwargs)
            return SimpleNamespace(output_parsed=expected, id="resp_mock_123")

    class FakeClient:
        def __init__(self, *, api_key: str) -> None:
            captured["api_key"] = api_key
            self.responses = FakeResponses()

    monkeypatch.setattr(openai, "OpenAI", FakeClient)

    actual, metadata = openai_interpret(
        packet,
        model="test-structured-output-model",
        api_key="not-a-real-key",
    )

    assert actual == expected
    validate_bundle(actual, packet)
    assert captured["model"] == "test-structured-output-model"
    assert captured["text_format"] is InterpretationBundle
    assert json.loads(captured["input"][1]["content"]) == packet
    assert "only the supplied JSON evidence" in captured["input"][0]["content"]
    assert metadata == {
        "backend": "openai_responses_structured_output",
        "model": "test-structured-output-model",
        "response_id": "resp_mock_123",
        "prompt_sha256": metadata["prompt_sha256"],
        "validated": True,
        "external_api_called": True,
    }
    assert len(metadata["prompt_sha256"]) == 64


@pytest.mark.parametrize(
    ("bad_label", "bad_evidence", "message"),
    [
        ("invented cell", "cluster:0:marker:MS4A1", "Unsupported label"),
        ("B cell", "invented:evidence", "Unknown evidence IDs"),
    ],
)
def test_mocked_openai_invalid_grounding_is_rejected(
    monkeypatch,
    bad_label: str,
    bad_evidence: str,
    message: str,
) -> None:
    """Structured parsing is necessary but not sufficient: grounding is local."""

    import openai

    packet = _packet()
    invalid = InterpretationBundle(
        interpretations=[
            ClusterInterpretation(
                cluster="0",
                proposed_label=bad_label,
                confidence="high",
                evidence_ids=[bad_evidence],
                uncertainties=[],
            ),
            ClusterInterpretation(
                cluster="1",
                proposed_label="NK cell",
                confidence="low",
                evidence_ids=["cluster:1:marker:NKG7"],
                uncertainties=["mock"],
            ),
        ]
    )

    class FakeResponses:
        def parse(self, **_kwargs):
            return SimpleNamespace(output_parsed=invalid, id="resp_invalid")

    class FakeClient:
        def __init__(self, *, api_key: str) -> None:
            assert api_key == "not-a-real-key"
            self.responses = FakeResponses()

    monkeypatch.setattr(openai, "OpenAI", FakeClient)

    with pytest.raises(ValueError, match=message):
        openai_interpret(
            packet,
            model="test-structured-output-model",
            api_key="not-a-real-key",
        )
