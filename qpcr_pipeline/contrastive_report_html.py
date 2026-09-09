"""Self-contained, read-only HTML report for contrastive conservation."""

from __future__ import annotations

import html
import json

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
    """Render deterministic offline HTML with read-only explorer controls."""
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
    candidate_regions = [
        {
            "id": item.region.region_id,
            "rank": item.region.rank,
            "start": item.region.reference_start,
            "end": item.region.reference_end,
            "peakStart": item.region.peak_start,
            "peakEnd": item.region.peak_end,
            "contributing": item.contributing_windows,
            "worstDataset": item.worst_dataset_name,
            "criticality": item.worst_dataset_criticality,
            "similarity": item.worst_similarity,
            "margin": item.contrast_margin,
            "meanConservation": item.region.mean_conservation,
            "minimumConservation": item.region.minimum_conservation,
        }
        for item in candidates
    ]
    dataset_names = sorted(
        {row.dataset_name for row in dataset_evidence}, key=str.casefold
    )

    candidate_rows = "".join(
        f'<tr class="candidate-row" data-region-id="{_cell(item.region.region_id)}" tabindex="0">'
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

    script = """<script>
const targetPoints=%s;
const challengePoints=%s;
const candidateRegions=%s;
const datasetNames=%s;
const visibleDatasets=new Set();
const hover=document.getElementById('hover');
const candidateDetail=document.getElementById('candidate-detail');
const filters=document.getElementById('dataset-filters');
const q=document.getElementById('quadrant'),qc=q.getContext('2d');
const t=document.getElementById('track'),tc=t.getContext('2d');
const plot={left:55,top:18,right:20,bottom:35};
const fullStart=targetPoints.length?Math.min(...targetPoints.map(p=>p.start)):0;
const fullEnd=targetPoints.length?Math.max(...targetPoints.map(p=>p.end)):1;
let viewStart=fullStart;
let viewEnd=Math.max(fullStart+1,fullEnd);
let selectedRegionId=null;
let overviewHits=[];
let trackHits=[];
let trackRegionHits=[];
let dragging=false;
let dragMoved=false;
let dragStartX=0;
let dragViewStart=0;
let dragViewEnd=0;

function axes(ctx,canvas,xLabel){
  ctx.clearRect(0,0,canvas.width,canvas.height);
  ctx.strokeStyle='#9aa5b2';
  ctx.lineWidth=1;
  ctx.beginPath();
  ctx.moveTo(plot.left,plot.top);
  ctx.lineTo(plot.left,canvas.height-plot.bottom);
  ctx.lineTo(canvas.width-plot.right,canvas.height-plot.bottom);
  ctx.stroke();
  ctx.fillStyle='#657181';
  ctx.font='12px system-ui,sans-serif';
  ctx.fillText('1.0',16,plot.top+5);
  ctx.fillText('0.0',16,canvas.height-plot.bottom+4);
  ctx.fillText(xLabel,plot.left,canvas.height-8);
}
function inView(p){return p.end>=viewStart&&p.start<=viewEnd;}
function overviewXY(canvas,p){
  return [
    plot.left+(canvas.width-plot.left-plot.right)*(p.x??0),
    plot.top+(canvas.height-plot.top-plot.bottom)*(1-(p.y??0))
  ];
}
function referenceX(canvas,position){
  const span=Math.max(1,viewEnd-viewStart);
  return plot.left+(canvas.width-plot.left-plot.right)*((position-viewStart)/span);
}
function metricY(canvas,value){
  return plot.top+(canvas.height-plot.top-plot.bottom)*(1-(value??0));
}
function drawPoint(ctx,x,y,r,fill,stroke){
  ctx.beginPath();ctx.arc(x,y,r,0,Math.PI*2);
  ctx.fillStyle=fill;ctx.fill();ctx.strokeStyle=stroke;ctx.stroke();
}
function drawOverview(){
  axes(qc,q,'Non-target similarity →');
  overviewHits=[];
  for(const p of targetPoints){
    if(!inView(p)||typeof p.x!=='number')continue;
    const [x,y]=overviewXY(q,p);
    drawPoint(qc,x,y,4,'#245e72','#173f4d');
    overviewHits.push({x,y,series:'worst',p});
  }
  for(const p of challengePoints){
    if(!inView(p)||!visibleDatasets.has(p.dataset))continue;
    const [x,y]=overviewXY(q,p);
    drawPoint(qc,x,y,2.5,'#fff','#8b5e34');
    overviewHits.push({x,y,series:'dataset',p});
  }
}
function drawLine(ctx,canvas,points,valueKey,stroke){
  const visible=points.filter(inView).filter(p=>typeof p[valueKey]==='number').sort((a,b)=>a.start-b.start);
  if(!visible.length)return;
  ctx.beginPath();
  visible.forEach((p,index)=>{
    const x=referenceX(canvas,(p.start+p.end)/2),y=metricY(canvas,p[valueKey]);
    if(index===0)ctx.moveTo(x,y);else ctx.lineTo(x,y);
  });
  ctx.strokeStyle=stroke;ctx.lineWidth=2;ctx.stroke();
}
function drawTrack(){
  axes(tc,t,'Reference position →');
  trackHits=[];
  trackRegionHits=[];
  const plotBottom=t.height-plot.bottom;
  for(const region of candidateRegions){
    if(region.end<viewStart||region.start>viewEnd)continue;
    const x1=referenceX(t,Math.max(viewStart,region.start));
    const x2=referenceX(t,Math.min(viewEnd,region.end));
    tc.fillStyle=region.id===selectedRegionId?'rgba(36,94,114,.20)':'rgba(36,94,114,.08)';
    tc.fillRect(x1,plot.top,Math.max(2,x2-x1),plotBottom-plot.top);
    tc.strokeStyle=region.id===selectedRegionId?'#245e72':'#8ca8b2';
    tc.lineWidth=region.id===selectedRegionId?2:1;
    tc.strokeRect(x1,plot.top,Math.max(2,x2-x1),plotBottom-plot.top);
    trackRegionHits.push({x1,x2,region});
  }
  drawLine(tc,t,targetPoints,'y','#245e72');
  drawLine(tc,t,targetPoints,'x','#9b4d2d');
  for(const name of visibleDatasets){
    drawLine(tc,t,challengePoints.filter(p=>p.dataset===name),'x','#8b5e34');
  }
  for(const p of targetPoints){
    if(!inView(p))continue;
    const mid=(p.start+p.end)/2;
    trackHits.push({x:referenceX(t,mid),y:metricY(t,p.y),series:'target',p});
    if(typeof p.x==='number')trackHits.push({x:referenceX(t,mid),y:metricY(t,p.x),series:'worst',p});
  }
  for(const p of challengePoints){
    if(!inView(p)||!visibleDatasets.has(p.dataset))continue;
    trackHits.push({x:referenceX(t,(p.start+p.end)/2),y:metricY(t,p.x),series:'dataset',p});
  }
  tc.fillStyle='#245e72';tc.fillText('Target conservation',70,32);
  tc.fillStyle='#9b4d2d';tc.fillText('Worst challenge similarity',190,32);
  if(visibleDatasets.size){tc.fillStyle='#8b5e34';tc.fillText('Selected dataset overlays',370,32);}
}
function redraw(){drawOverview();drawTrack();}
function nearestHit(hits,mx,my){
  let best=null,bestD=144;
  for(const hit of hits){
    const d=(hit.x-mx)**2+(hit.y-my)**2;
    if(d<bestD){best=hit;bestD=d;}
  }
  return best;
}
function hoverText(hit){
  if(!hit)return 'Move over a point for evidence details.';
  const p=hit.p;
  const dataset=hit.series==='target'?'target conservation':(p.dataset||'worst challenge');
  const bits=[dataset,`${p.start}-${p.end}`];
  if(typeof p.x==='number')bits.push(`similarity ${Number(p.x).toFixed(3)}`);
  if(typeof p.y==='number')bits.push(`target conservation ${Number(p.y).toFixed(3)}`);
  if(p.sequence)bits.push(`sequence ${p.sequence}`);
  if(p.orientation)bits.push(`orientation ${p.orientation}`);
  return bits.join(' · ');
}
function mousePoint(canvas,event){
  const r=canvas.getBoundingClientRect();
  return [(event.clientX-r.left)*canvas.width/r.width,(event.clientY-r.top)*canvas.height/r.height];
}
function referenceAtEvent(canvas,event){
  const [mx]=mousePoint(canvas,event);
  const ratio=Math.max(0,Math.min(1,(mx-plot.left)/(canvas.width-plot.left-plot.right)));
  return viewStart+ratio*(viewEnd-viewStart);
}
function zoomAt(referencePosition,factor){
  const fullWidth=Math.max(1,fullEnd-fullStart);
  const width=Math.max(1,viewEnd-viewStart);
  const next=Math.max(Math.min(100,fullWidth),Math.min(fullWidth,width*factor));
  const ratio=(referencePosition-viewStart)/width;
  viewStart=Math.max(fullStart,referencePosition-next*ratio);
  viewEnd=Math.min(fullEnd,viewStart+next);
  viewStart=Math.max(fullStart,viewEnd-next);
  redraw();
}
function selectRegion(regionId){
  selectedRegionId=regionId;
  document.querySelectorAll('.candidate-row').forEach(row=>{
    row.classList.toggle('selected',row.dataset.regionId===regionId);
  });
  const region=candidateRegions.find(item=>item.id===regionId);
  candidateDetail.textContent=region
    ? `${region.id} · rank ${region.rank} · region ${region.start}-${region.end} · anchor ${region.peakStart}-${region.peakEnd} · worst ${region.worstDataset??'n/a'} · similarity ${region.similarity??'n/a'} · margin ${region.margin??'n/a'} · mean conservation ${region.meanConservation} · minimum conservation ${region.minimumConservation}`
    : 'Select a candidate region to inspect its recorded evidence.';
  drawTrack();
}
function setAllDatasets(visible){
  visibleDatasets.clear();
  filters.querySelectorAll('input[type="checkbox"]').forEach(input=>{
    input.checked=visible;
    if(visible)visibleDatasets.add(input.value);
  });
  redraw();
}
for(const name of datasetNames){
  const label=document.createElement('label');
  const input=document.createElement('input');
  input.type='checkbox';input.value=name;
  input.addEventListener('change',()=>{
    if(input.checked)visibleDatasets.add(name);else visibleDatasets.delete(name);
    redraw();
  });
  label.append(input,document.createTextNode(` ${name}`));
  filters.appendChild(label);
}
document.getElementById('show-all').addEventListener('click',()=>setAllDatasets(true));
document.getElementById('hide-all').addEventListener('click',()=>setAllDatasets(false));
document.getElementById('zoom-in').addEventListener('click',()=>zoomAt((viewStart+viewEnd)/2,0.7));
document.getElementById('zoom-out').addEventListener('click',()=>zoomAt((viewStart+viewEnd)/2,1/0.7));
document.getElementById('reset-view').addEventListener('click',()=>{viewStart=fullStart;viewEnd=Math.max(fullStart+1,fullEnd);redraw();});
q.addEventListener('mousemove',event=>{
  const [mx,my]=mousePoint(q,event);hover.textContent=hoverText(nearestHit(overviewHits,mx,my));
});
q.addEventListener('mouseleave',()=>{hover.textContent='Move over a point for evidence details.';});
t.addEventListener('mousemove',event=>{
  if(dragging)return;
  const [mx,my]=mousePoint(t,event);hover.textContent=hoverText(nearestHit(trackHits,mx,my));
});
t.addEventListener('mouseleave',()=>{if(!dragging)hover.textContent='Move over a point for evidence details.';});
t.addEventListener('wheel',event=>{
  event.preventDefault();
  zoomAt(referenceAtEvent(t,event),event.deltaY<0?0.7:1/0.7);
},{passive:false});
t.addEventListener('pointerdown',event=>{
  dragging=true;dragMoved=false;dragStartX=event.clientX;dragViewStart=viewStart;dragViewEnd=viewEnd;
  t.style.cursor='grabbing';t.setPointerCapture(event.pointerId);
});
t.addEventListener('pointermove',event=>{
  if(!dragging)return;
  const rect=t.getBoundingClientRect();
  const delta=event.clientX-dragStartX;
  if(Math.abs(delta)>2)dragMoved=true;
  const width=dragViewEnd-dragViewStart;
  let start=dragViewStart-delta*(width/Math.max(1,rect.width));
  start=Math.max(fullStart,Math.min(fullEnd-width,start));
  viewStart=start;viewEnd=start+width;redraw();
});
t.addEventListener('pointerup',event=>{
  if(!dragging)return;
  dragging=false;t.style.cursor='grab';
  const [mx]=mousePoint(t,event);
  if(!dragMoved){
    const hit=trackRegionHits.find(item=>mx>=item.x1&&mx<=item.x2);
    if(hit)selectRegion(hit.region.id);
  }
  try{t.releasePointerCapture(event.pointerId);}catch(_error){}
});
t.addEventListener('pointercancel',()=>{dragging=false;t.style.cursor='grab';});
document.querySelectorAll('.candidate-row').forEach(row=>{
  const activate=()=>selectRegion(row.dataset.regionId);
  row.addEventListener('click',activate);
  row.addEventListener('keydown',event=>{
    if(event.key==='Enter'||event.key===' '){event.preventDefault();activate();}
  });
});
redraw();
</script>""" % (
        _safe_json(target_points),
        _safe_json(challenge_points),
        _safe_json(candidate_regions),
        _safe_json(dataset_names),
    )

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
.card{{border:1px solid #ddd;border-radius:12px;padding:14px}} canvas{{width:100%;height:320px;border:1px solid #ddd;border-radius:10px;touch-action:none}}
table{{width:100%;border-collapse:collapse;font-size:14px}} th,td{{text-align:left;border-bottom:1px solid #e7e7e7;padding:8px;vertical-align:top}}
#hover{{min-height:2em;padding:8px 0;color:#333}} .toolbar{{display:flex;gap:8px;flex-wrap:wrap;margin:12px 0}} button{{font:inherit;padding:7px 10px;border:1px solid #bbb;border-radius:8px;background:#fff;cursor:pointer}}
.filters{{display:flex;gap:10px 16px;flex-wrap:wrap;margin:8px 0 16px}} .filters label{{white-space:nowrap}} #track{{cursor:grab}}
.candidate-row{{cursor:pointer}} .candidate-row:hover,.candidate-row:focus{{background:#f4f7f8;outline:none}} .candidate-row.selected{{background:#e7f0f3}}
.legend{{display:flex;gap:18px;flex-wrap:wrap;font-size:13px;color:#555;margin:6px 0 16px}}
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
<div class="toolbar" aria-label="Explorer controls">
<button id="zoom-in" type="button">Zoom in</button>
<button id="zoom-out" type="button">Zoom out</button>
<button id="reset-view" type="button">Reset view</button>
<button id="show-all" type="button">Show all datasets</button>
<button id="hide-all" type="button">Hide all datasets</button>
</div>
<div id="dataset-filters" class="filters" aria-label="Challenge dataset overlays"></div>
<h2>Contrast overview</h2>
<p class="muted">Default view shows one point per target window using its worst recorded challenge similarity. Enable datasets above to overlay their individual evidence.</p>
<canvas id="quadrant" width="1100" height="360" role="img" aria-label="Target conservation versus non-target similarity"></canvas>
<h2>Reference track</h2>
<p class="muted">Use the mouse wheel or buttons to zoom, drag to pan, and click a candidate band to inspect it.</p>
<canvas id="track" width="1100" height="320" role="img" aria-label="Reference-position conservation and challenge similarity"></canvas>
<div class="legend"><span>Target conservation</span><span>Worst challenge similarity</span><span>Optional per-dataset overlays</span></div>
<div id="hover" aria-live="polite">Move over a point for evidence details.</div>
<div id="candidate-detail" class="card" aria-live="polite">Select a candidate region to inspect its recorded evidence.</div>
<h2>Candidate regions</h2>
<table><thead><tr><th>ID</th><th>Rank</th><th>Reference</th><th>Contributing windows</th><th>Worst challenge</th><th>Criticality</th><th>Similarity</th></tr></thead><tbody>{candidate_rows}</tbody></table>
<h2>Per-dataset evidence</h2>
<table><thead><tr><th>Window</th><th>Dataset</th><th>Criticality</th><th>Sequences</th><th>Best sequence</th><th>Orientation</th><th>Similarity</th></tr></thead><tbody>{dataset_rows}</tbody></table>
{script}
</body></html>"""
