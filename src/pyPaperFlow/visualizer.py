"""Simple data visualization exporters for analysis outputs."""

from __future__ import annotations

import json
from typing import Dict, List, Tuple


class Visualizer:
    """Export analysis results to interoperable text formats."""

    def to_markdown_table(self, title: str, headers: List[str], rows: List[List[str]]) -> str:
        lines = [f"## {title}", ""]
        lines.append("| " + " | ".join(headers) + " |")
        lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
        for row in rows:
            lines.append("| " + " | ".join(str(c) for c in row) + " |")
        lines.append("")
        return "\n".join(lines)

    def export_cytoscape_json(self, edges: List[Tuple[str, str]], output_path: str) -> None:
        payload = {
            "nodes": [],
            "edges": [],
        }
        seen = set()
        for src, dst in edges:
            if src not in seen:
                payload["nodes"].append({"data": {"id": src, "label": src}})
                seen.add(src)
            if dst not in seen:
                payload["nodes"].append({"data": {"id": dst, "label": dst}})
                seen.add(dst)
            payload["edges"].append({"data": {"source": src, "target": dst}})

        with open(output_path, "w", encoding="utf-8") as f:
            json.dump(payload, f, ensure_ascii=False, indent=2)
