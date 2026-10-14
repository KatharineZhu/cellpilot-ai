"""Evidence-bounded interpretation with deterministic and optional OpenAI backends."""

from __future__ import annotations

import hashlib
import json
import os
from pathlib import Path
from typing import Any, Literal

import pandas as pd
import yaml
from pydantic import BaseModel, ConfigDict, Field


CELLTYPIST_MODEL_REFERENCE = {
    "title": "Healthy_COVID19_PBMC CellTypist model source (Stephenson et al., 2021)",
    "url": "https://doi.org/10.1038/s41591-021-01329-2",
}


class ClusterInterpretation(BaseModel):
    model_config = ConfigDict(extra="forbid")

    cluster: str
    proposed_label: str
    confidence: Literal["low", "moderate", "high"]
    evidence_ids: list[str]
    uncertainties: list[str]


class InterpretationBundle(BaseModel):
    model_config = ConfigDict(extra="forbid")

    interpretations: list[ClusterInterpretation]


def _clean_number(value: Any) -> float | None:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return None
    if not pd.notna(number):
        return None
    return number


def _cluster_identifier(value: Any) -> str:
    """Normalize identifiers that pandas may coerce from ``0`` to ``0.0``."""

    if pd.isna(value):
        raise ValueError("cluster identifiers must not be missing")
    text = str(value).strip()
    try:
        numeric = float(text)
    except ValueError:
        return text
    if numeric.is_integer():
        return str(int(numeric))
    return text


def load_marker_registry(path: Path) -> dict[str, Any]:
    payload = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(payload, dict) or "panels" not in payload:
        raise ValueError(f"Invalid marker registry: {path}")
    return payload


def build_evidence_packet(
    markers: pd.DataFrame,
    cluster_counts: pd.DataFrame,
    *,
    marker_registry: dict[str, Any],
    annotation_consensus: pd.DataFrame | None = None,
    top_n: int = 10,
) -> dict[str, Any]:
    """Build a JSON-safe packet whose IDs are the only allowed claim evidence."""

    counts = {
        _cluster_identifier(row["cluster"]): int(row["n_cells"])
        for _, row in cluster_counts.iterrows()
    }
    consensus_by_cluster: dict[str, dict[str, Any]] = {}
    if annotation_consensus is not None and not annotation_consensus.empty:
        for _, row in annotation_consensus.iterrows():
            cluster = _cluster_identifier(
                row.get("cluster", row.get("leiden", ""))
            )
            label = row.get("predicted_label", row.get("majority_label"))
            if cluster and pd.notna(label):
                consensus_by_cluster[cluster] = {
                    "label": str(label),
                    "fraction": _clean_number(
                        row.get("fraction", row.get("label_fraction", 0.0))
                    ),
                    "mean_confidence": _clean_number(row.get("mean_confidence")),
                }

    panels = marker_registry["panels"]
    packet_clusters: list[dict[str, Any]] = []
    marker_cluster_ids = markers["cluster"].map(_cluster_identifier)
    for cluster, frame in markers.groupby(marker_cluster_ids, sort=True):
        frame = frame.sort_values("rank").head(top_n)
        evidence: list[dict[str, Any]] = []
        top_genes: list[str] = []
        for _, row in frame.iterrows():
            gene = str(row["names"])
            top_genes.append(gene)
            evidence.append(
                {
                    "evidence_id": f"cluster:{cluster}:marker:{gene}",
                    "kind": "computed_marker",
                    "gene": gene,
                    "rank": int(row["rank"]),
                    "score": _clean_number(row.get("scores")),
                    "logfoldchange": _clean_number(row.get("logfoldchanges")),
                    "adjusted_p_value": _clean_number(row.get("pvals_adj")),
                    "pct_group": _clean_number(row.get("pct_nz_group")),
                    "pct_reference": _clean_number(row.get("pct_nz_reference")),
                }
            )

        upper_genes = {gene.upper() for gene in top_genes}
        panel_matches: dict[str, list[str]] = {}
        for label, panel in panels.items():
            matches = [
                gene
                for gene in panel.get("markers", [])
                if str(gene).upper() in upper_genes
            ]
            if matches:
                panel_matches[str(label)] = matches

        allowed_labels = ["unresolved", *sorted(panel_matches)]
        annotation = consensus_by_cluster.get(str(cluster))
        if annotation:
            allowed_labels.append(annotation["label"])
            evidence.append(
                {
                    "evidence_id": f"cluster:{cluster}:celltypist",
                    "kind": "celltypist_cluster_consensus",
                    **annotation,
                }
            )
        allowed_labels = list(dict.fromkeys(allowed_labels))
        packet_clusters.append(
            {
                "cluster": str(cluster),
                "n_cells": counts.get(str(cluster)),
                "top_genes": top_genes,
                "panel_matches": panel_matches,
                "allowed_labels": allowed_labels,
                "evidence": evidence,
            }
        )

    external_references = [marker_registry.get("source", {})]
    if consensus_by_cluster:
        external_references.append(CELLTYPIST_MODEL_REFERENCE)
    return {
        "schema_version": 1,
        "scientific_status": "exploratory interpretation; not biological validation",
        "marker_registry": marker_registry.get("source", {}),
        "external_references": external_references,
        "clusters": packet_clusters,
    }


def deterministic_interpret(packet: dict[str, Any]) -> InterpretationBundle:
    """Render only observations encoded in the packet, with explicit uncertainty."""

    results: list[ClusterInterpretation] = []
    for cluster in packet["clusters"]:
        annotation = next(
            (
                item
                for item in cluster["evidence"]
                if item["kind"] == "celltypist_cluster_consensus"
            ),
            None,
        )
        if annotation and (annotation.get("fraction") or 0) >= 0.5:
            label = str(annotation["label"])
            evidence_ids = [annotation["evidence_id"]]
            confidence = "moderate"
            uncertainties = [
                "CellTypist agreement is a model prediction, not experimental ground truth."
            ]
        else:
            matches = cluster.get("panel_matches", {})
            ordered = sorted(matches.items(), key=lambda item: (-len(item[1]), item[0]))
            if ordered and len(ordered[0][1]) >= 2:
                label, genes = ordered[0]
                evidence_ids = [
                    f"cluster:{cluster['cluster']}:marker:{gene}" for gene in genes
                ]
                confidence = "low"
                uncertainties = [
                    "The label is a heuristic marker-panel match and requires orthogonal validation."
                ]
            else:
                label = "unresolved"
                evidence_ids = [
                    item["evidence_id"]
                    for item in cluster["evidence"][: min(3, len(cluster["evidence"]))]
                ]
                confidence = "low"
                uncertainties = [
                    "Available evidence did not support a unique broad cell-type label."
                ]
        results.append(
            ClusterInterpretation(
                cluster=str(cluster["cluster"]),
                proposed_label=label,
                confidence=confidence,
                evidence_ids=evidence_ids,
                uncertainties=uncertainties,
            )
        )
    return InterpretationBundle(interpretations=results)


def validate_bundle(bundle: InterpretationBundle, packet: dict[str, Any]) -> None:
    """Reject labels or citations not explicitly available in the evidence packet."""

    packet_by_cluster = {item["cluster"]: item for item in packet["clusters"]}
    seen: set[str] = set()
    for result in bundle.interpretations:
        if result.cluster in seen:
            raise ValueError(f"Duplicate interpretation for cluster {result.cluster}")
        seen.add(result.cluster)
        if result.cluster not in packet_by_cluster:
            raise ValueError(f"Unknown cluster in interpretation: {result.cluster}")
        cluster = packet_by_cluster[result.cluster]
        if result.proposed_label not in cluster["allowed_labels"]:
            raise ValueError(
                f"Unsupported label for cluster {result.cluster}: {result.proposed_label}"
            )
        valid_ids = {item["evidence_id"] for item in cluster["evidence"]}
        invalid = set(result.evidence_ids) - valid_ids
        if invalid:
            raise ValueError(
                f"Unknown evidence IDs for cluster {result.cluster}: {sorted(invalid)}"
            )
        if not result.evidence_ids:
            raise ValueError(f"Cluster {result.cluster} has no evidence IDs")
    expected = set(packet_by_cluster)
    if seen != expected:
        raise ValueError(f"Missing cluster interpretations: {sorted(expected - seen)}")


def openai_interpret(
    packet: dict[str, Any],
    *,
    model: str,
    api_key: str | None = None,
    api_key_env: str = "OPENAI_API_KEY",
) -> tuple[InterpretationBundle, dict[str, Any]]:
    """Use Responses Structured Outputs, then enforce the local evidence contract."""
    if not api_key_env:
        raise ValueError("api_key_env must be non-empty")
    key = api_key or os.getenv(api_key_env)
    if not key:
        raise RuntimeError(
            f"OpenAI interpretation was requested, but {api_key_env} is not set. "
            "No external API request was made."
        )
    if not model:
        raise ValueError("An explicit OpenAI model is required")
    try:
        from openai import OpenAI
    except ImportError as exc:  # pragma: no cover - dependency is present in CI
        raise RuntimeError(
            "OpenAI interpretation requires the 'openai' Python package. "
            "No external API request was made."
        ) from exc
    instructions = (
        "Interpret each cluster using only the supplied JSON evidence. Choose a proposed_label "
        "only from that cluster's allowed_labels and cite only its evidence_id values. "
        "Do not invent genes, measurements, citations, or biological facts. Use low confidence "
        "and explicit uncertainty when the evidence is ambiguous. Return every cluster exactly once."
    )
    prompt = json.dumps(packet, sort_keys=True, separators=(",", ":"))
    client = OpenAI(api_key=key)
    response = client.responses.parse(
        model=model,
        input=[
            {"role": "system", "content": instructions},
            {"role": "user", "content": prompt},
        ],
        text_format=InterpretationBundle,
    )
    bundle = response.output_parsed
    if bundle is None:
        raise RuntimeError("OpenAI response did not contain parsed structured output")
    validate_bundle(bundle, packet)
    metadata = {
        "backend": "openai_responses_structured_output",
        "model": model,
        "response_id": response.id,
        "prompt_sha256": hashlib.sha256(prompt.encode("utf-8")).hexdigest(),
        "validated": True,
        "external_api_called": True,
    }
    return bundle, metadata


def write_interpretation_artifacts(
    run_dir: Path,
    packet: dict[str, Any],
    bundle: InterpretationBundle,
    *,
    metadata: dict[str, Any],
) -> None:
    output_dir = run_dir / "interpretation"
    output_dir.mkdir(parents=True, exist_ok=True)
    validate_bundle(bundle, packet)
    (output_dir / "evidence_packet.json").write_text(
        json.dumps(packet, indent=2) + "\n", encoding="utf-8"
    )
    payload = {"metadata": metadata, **bundle.model_dump()}
    (output_dir / "interpretation.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8"
    )
    lines = [
        "# Evidence-grounded cluster interpretation",
        "",
        "> Exploratory output only. Labels are hypotheses derived from computed markers and, when available, a CellTypist prediction; they are not experimental ground truth.",
        "",
    ]
    for result in bundle.interpretations:
        lines.extend(
            [
                f"## Cluster {result.cluster}: {result.proposed_label}",
                "",
                f"Confidence: **{result.confidence}**",
                "",
                "Evidence IDs: " + ", ".join(f"`{item}`" for item in result.evidence_ids),
                "",
                "Uncertainty: " + " ".join(result.uncertainties),
                "",
            ]
        )
    references = [
        item
        for item in packet.get("external_references", [])
        if isinstance(item, dict) and item.get("title") and item.get("url")
    ]
    if references:
        lines.extend(["## Traceable external sources", ""])
        lines.extend(
            f"- [{item['title']}]({item['url']})" for item in references
        )
        lines.append("")
    (output_dir / "interpretation.md").write_text("\n".join(lines), encoding="utf-8")
