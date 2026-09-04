"""Self-contained, read-only HTML report for contrastive conservation."""

from __future__ import annotations

import html
import json
from dataclasses import asdict

from qpcr_pipeline.contrastive_conservation import (
    ContrastCandidateRegion,
    ContrastWindowEvidence,
    DatasetWindowEvidence,
)


def _safe_json(value: object) -> str:
    return (
        json.dumps(value, separators=(",", ":"), sort_keys=True)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
    )


def _cell(value: object) -> str:
    return html.escape("" if value is None else str(value), quote=True)


def render_contrastive_html(
    *,
    target_name: str,
    reference_id: str | None,
    windows: tuple[ContrastWindowEvidence, ...],
    dataset_evidence: tuple[DatasetWindowEvidence, ...],
    candidates: tuple[ContrastCandidateRegion, ...],
) -> str:
    """Render deterministic HTML with no network dependencies or edit controls."""
    target_points = [
        {
            "start": row.reference_start,
            "end": row.reference_end,
            "x": row.worst_similarity,
            "y": row.target_mean_conservation,
            "critical": row.worst_critical_similarity,
            "important": row.worst_important_similarity,
            "dataset": row.worst_dataset_name,
        }
        for row in windows
    ]
    target_lookup = {
        (row.reference_start, row.reference_end): row.target_mean_conservation
        for row in windows
    }
    challenge_points = [
        {
            "start": row.reference_start,
            "end": row.reference_end,
            "x": row.similarity,
            "y": target_lookup.get((row.reference_start, row.reference_end), 0.0),
            "dataset": row.dataset_name,
            "criticality": row.criticality,
            "sequence": row.best_sequence_id,
            "orientation": row.best_orientation,
        }
        for row in dataset_evidence
    ]

    candidate_rows = "".join(
        "<tr>"
        f"<td>{_cell(item.region.region_id)}</td>"
        f"<td>{item.region.rank}</td>"
        f"<td>{item.region.reference_start}-{item.region.reference_end}</td>"
        f"<td>{_cell(item.contributing_windows)}</td>"
        f"<td>{_cell(item.worst_dataset_name)}</td>"
        f"<td>{_cell(item.worst_dataset_criticality)}</td>"
        f"<td>{_cell(item.worst_similarity)}</td>"
        "</tr>"
        for item in candidates
    ) or '<tr><td colspan="7">No candidate regions.</td></tr>'

    dataset_rows = "".join(
        "<tr>"
        f"<td>{row.reference_start}-{row.reference_end}</td>"
        f"<td>{_cell(row.dataset_name)}</td>"
        f"<td>{_cell(row.criticality)}</td>"
        f"<td>{row.sequence_count}</td>"
        f"<td>{_cell(row.best_sequence_id)}</td>"
        f"<td>{_cell(row.best_orientation)}</td>"
        f"<td>{row.similarity:.4f}</td>"
        "</tr>"
        for row in dataset_evidence
    ) or '<tr><td colspan="7">No dataset evidence.</td></tr>'

    return f"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Target vs non-target contrast</title>
<style>
:root{{font-family:system-ui,sans-serif;color:#171717;background:#fff}}
body{{max-width:1180px;margin:0 auto;padding:28px;line-height:1.45}}
h1,h2{{letter-spacing:-.02em}} .muted{{color:#666}} .grid{{display:grid;grid-template-columns:repeat(auto-fit,minmax(220px,1fr));gap:12px}}
.card{{border:1px solid #ddd;border-radius:12px;padding:14px}} canvas{{width:100%;height:320px;border:1px solid #ddd;border-radius:10px}}
table{{width:100%;border-collapse:collapse;font-size:14px}} th,td{{text-align:left;border-bottom:1px solid #e7e7e7;padding:8px;vertical-align:top}}
#hover{{min-height:2em;padding:8px 0;color:#333}}
</style>
</head>
<body>
<h1>Target vs non-target contrast</h1>
<p class="muted">Target: {_cell(target_name)} · Reference: {_cell(reference_id or "n/a")}</p>
<div class="grid">
<div class="card"><strong>Target windows</strong><div>{len(windows)}</div></div>
<div class="card"><strong>Challenge evidence rows</strong><div>{len(dataset_evidence)}</div></div>
<div class="card"><strong>Candidate regions</strong><div>{len(candidates)}</div></div>
</div>
<h2>Contrast overview</h2>
<canvas id="quadrant" width="1000" height="360"></canvas>
<h2>Reference track</h2>
<canvas id="track" width="1000" height="300"></canvas>
<div id="hover" aria-live="polite">Move over a point for read-only evidence details.</div>
<h2>Candidate regions</h2>
<table><thead><tr><th>ID</th><th>Rank</th><th>Reference</th><th>Contributing windows</th><th>Worst challenge</th><th>Criticality</th><th>Similarity</th></tr></thead><tbody>{candidate_rows}</tbody></table>
<h2>Per-dataset evidence</h2>
<table><thead><tr><th>Window</th><th>Dataset</th><th>Criticality</th><th>Sequences</th><th>Best sequence</th><th>Orientation</th><th>Similarity</th></tr></thead><tbody>{dataset_rows}</tbody></table>
<script>
const targetPoints={_safe_json(target_points)};
const challengePoints={_safe_json(challenge_points)};
const hover=document.getElementById('hover');
function axes(ctx,w,h){{ctx.clearRect(0,0,w,h);ctx.strokeStyle='#bbb';ctx.beginPath();ctx.moveTo(55,15);ctx.lineTo(55,h-35);ctx.lineTo(w-20,h-35);ctx.stroke();}}
function xy(canvas,p){{return [55+(canvas.width-75)*(p.x??0),15+(canvas.height-50)*(1-(p.y??0))];}}
const q=document.getElementById('quadrant'),qc=q.getContext('2d');axes(qc,q.width,q.height);
function drawPoints(ctx,canvas,points,r){{for(const p of points){{const [x,y]=xy(canvas,p);ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);ctx.stroke();}}}}
drawPoints(qc,q,targetPoints,5);drawPoints(qc,q,challengePoints,3);
const t=document.getElementById('track'),tc=t.getContext('2d');axes(tc,t.width,t.height);
const allStarts=targetPoints.map(p=>p.start);const maxRef=Math.max(1,...targetPoints.map(p=>p.end));
function trackPoint(canvas,p,value){{return [55+(canvas.width-75)*(p.start/maxRef),15+(canvas.height-50)*(1-value)];}}
for(const p of targetPoints){{const [x,y]=trackPoint(t,p,p.y);tc.strokeRect(x-3,y-3,6,6);}}
for(const p of challengePoints){{const [x,y]=trackPoint(t,p,p.x);tc.strokeRect(x-2,y-2,4,4);}}
q.addEventListener('mousemove',e=>{{
  const r=q.getBoundingClientRect(),mx=(e.clientX-r.left)*q.width/r.width,my=(e.clientY-r.top)*q.height/r.height;
  let best=null,bestD=144;
  for(const p of targetPoints){{const [x,y]=xy(q,p),d=(x-mx)**2+(y-my)**2;if(d<bestD){{best={{series:'target',p}},bestD=d;}}}}
  for(const p of challengePoints){{const [x,y]=xy(q,p),d=(x-mx)**2+(y-my)**2;if(d<bestD){{best={{series:'challenge',p}},bestD=d;}}}}
  hover.textContent=best?`${{best.series}} · ${{best.p.start}}-${{best.p.end}} · similarity ${{best.p.x??'n/a'}} · target conservation ${{best.p.y??'n/a'}}`:'Move over a point for read-only evidence details.';
}});
</script>
</body></html>"""
