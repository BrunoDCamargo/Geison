"""Deterministic, self-contained HTML for conservation results."""

from __future__ import annotations

import json

from qpcr_pipeline.config import ConservationConfig


def _window_values(item: object) -> list[int | float]:
    return [
        item.reference_start,
        item.reference_end,
        item.position_count,
        item.mean_conservation,
        item.minimum_conservation,
        item.mean_coverage,
        item.mean_gap_frequency,
        item.mean_entropy_bits,
    ]


def _safe_json(payload: object) -> str:
    """Serialize JSON so data cannot terminate its HTML script element."""
    return (
        json.dumps(payload, ensure_ascii=True, separators=(",", ":"), sort_keys=True)
        .replace("<", "\\u003c")
        .replace(">", "\\u003e")
        .replace("&", "\\u0026")
        .replace("\u2028", "\\u2028")
        .replace("\u2029", "\\u2029")
    )


def render_conservation_html(
    *,
    target_name: str,
    reference_id: str | None,
    sequence_count: int,
    config: ConservationConfig,
    windows: tuple[object, ...],
    annotations: tuple[object, ...],
) -> str:
    """Render an offline Canvas report with deterministic ranking and controls."""
    window_rows = [_window_values(item) for item in windows]
    top_windows = sorted(
        windows,
        key=lambda item: (
            -item.mean_conservation,
            -item.minimum_conservation,
            -item.mean_coverage,
            item.mean_entropy_bits,
            item.reference_start,
        ),
    )[:10]
    payload = {
        "identity": {
            "targetName": target_name,
            "referenceId": reference_id,
            "sequenceCount": sequence_count,
            "windowCount": len(windows),
            "windowSize": config.window_size,
            "stepSize": config.step_size,
        },
        "windows": window_rows,
        "topWindows": [_window_values(item) for item in top_windows],
        "annotations": [
            [item.feature_type, item.start, item.end, item.strand, item.label]
            for item in annotations
        ],
    }
    empty_hidden = "" if not windows else " hidden"
    return (
        _REPORT_TEMPLATE.replace("__REPORT_DATA__", _safe_json(payload)).replace(
            "__EMPTY_HIDDEN__", empty_hidden
        )
    )


_REPORT_TEMPLATE = """<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Geison conservation report</title>
<style>
:root{color-scheme:light;--ink:#172033;--muted:#596579;--line:#d8deea;--panel:#f5f7fb;--conservation:#1769aa;--coverage:#d65f00;--annotation:#805ad5;--peak:#fff0b3}
*{box-sizing:border-box}
body{font-family:system-ui,-apple-system,"Segoe UI",sans-serif;max-width:76rem;margin:0 auto;padding:1.5rem;color:var(--ink);background:#fff}
h1,h2{line-height:1.2}h1{margin:.2rem 0 1rem}h2{font-size:1.15rem;margin:0 0 .8rem}
.summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(11rem,1fr));gap:.75rem;margin:0 0 1.25rem}
.summary div,.panel{border:1px solid var(--line);border-radius:.6rem;background:var(--panel);padding:.8rem}
dt{font-size:.78rem;text-transform:uppercase;letter-spacing:.04em;color:var(--muted);font-weight:700}dd{margin:.25rem 0 0;font-weight:650;overflow-wrap:anywhere}
.toolbar{display:flex;align-items:center;flex-wrap:wrap;gap:.75rem;margin:.6rem 0}.toolbar button{border:1px solid #9aa7ba;border-radius:.4rem;background:#fff;padding:.45rem .7rem;cursor:pointer}.toolbar button:hover,.toolbar button:focus-visible{background:#e9eef7}
.legend{display:flex;gap:1rem;flex-wrap:wrap;font-size:.9rem}.key::before{content:"";display:inline-block;width:1.4rem;border-top:3px solid;margin-right:.35rem;vertical-align:middle}.key.conservation::before{border-color:var(--conservation)}.key.coverage::before{border-color:var(--coverage)}
.canvas-wrap{position:relative;width:100%;height:29rem;border:1px solid var(--line);border-radius:.5rem;background:#fff;overflow:hidden}canvas{display:block;width:100%;height:100%;touch-action:none;cursor:crosshair}
.instructions{color:var(--muted);font-size:.88rem;margin:.55rem 0 0}.workspace{display:grid;grid-template-columns:minmax(0,2fr) minmax(18rem,1fr);gap:1rem;margin-top:1rem}.stack{display:grid;gap:1rem;align-content:start}
#hover-details dl{display:grid;grid-template-columns:max-content 1fr;gap:.3rem .75rem;margin:0}#hover-details dt{text-transform:none;letter-spacing:0}#hover-details dd{margin:0;font-weight:500}
table{border-collapse:collapse;width:100%;font-variant-numeric:tabular-nums;font-size:.86rem}th,td{text-align:right;padding:.38rem;border-bottom:1px solid var(--line)}th:first-child,td:first-child{text-align:left}tbody tr{cursor:pointer}tbody tr:hover,tbody tr:focus{background:#fff}tbody tr:focus{outline:2px solid #1769aa;outline-offset:-2px}
#empty-state{border:1px solid #e0b000;background:#fff8d6;border-radius:.5rem;padding:1rem;margin:1rem 0}
@media(max-width:52rem){body{padding:1rem}.workspace{grid-template-columns:1fr}.canvas-wrap{height:24rem}}
</style>
</head>
<body>
<main>
<h1>Genome conservation report</h1>
<dl class="summary">
<div><dt>Target</dt><dd id="target"></dd></div>
<div><dt>Reference</dt><dd id="reference"></dd></div>
<div><dt>Discovery sequences</dt><dd id="sequence-count"></dd></div>
<div><dt>Windows</dt><dd id="window-summary"></dd></div>
</dl>
<p id="empty-state" role="status" aria-live="polite"__EMPTY_HIDDEN__>No conservation windows are available.</p>
<section class="panel" aria-labelledby="plot-heading">
<h2 id="plot-heading">Conservation across the reference genome</h2>
<div class="toolbar">
<button id="reset-zoom" type="button">Reset zoom</button>
<span id="view-interval" role="status" aria-live="polite"></span>
<span class="legend" aria-label="Plot legend"><span class="key conservation">Conservation</span><span class="key coverage">Coverage</span></span>
</div>
<div class="canvas-wrap"><canvas id="conservation-canvas" aria-label="Genome conservation and coverage plot"></canvas></div>
<p class="instructions"><strong>Zoom:</strong> use the mouse wheel over the plot. <strong>Pan:</strong> drag horizontally. <strong>Reset:</strong> use Reset zoom. Select a top window to focus it.</p>
</section>
<div class="workspace">
<section id="hover-details" class="panel" aria-labelledby="details-heading">
<h2 id="details-heading">Window details</h2>
<dl>
<dt>Interval</dt><dd id="detail-interval">Hover over the plot</dd>
<dt>Mean conservation</dt><dd id="detail-mean">Not selected</dd>
<dt>Minimum conservation</dt><dd id="detail-minimum">Not selected</dd>
<dt>Mean coverage</dt><dd id="detail-coverage">Not selected</dd>
<dt>Mean gaps</dt><dd id="detail-gaps">Not selected</dd>
<dt>Mean entropy</dt><dd id="detail-entropy">Not selected</dd>
<dt>Overlapping annotations</dt><dd id="detail-annotations">None</dd>
</dl>
</section>
<div class="stack">
<section class="panel" aria-labelledby="annotations-heading"><h2 id="annotations-heading">Reference annotations</h2><p>Feature spans are drawn in the annotation lane below the genomic traces. Hover a window to list overlapping labels.</p></section>
<section class="panel" aria-labelledby="top-heading">
<h2 id="top-heading">Top conservation windows</h2>
<div style="overflow-x:auto"><table id="top-windows"><thead><tr><th>Interval</th><th>Mean</th><th>Minimum</th><th>Coverage</th><th>Entropy</th></tr></thead><tbody></tbody></table></div>
</section>
</div>
</div>
<p class="instructions">Per-position metrics and reference-coordinate consensus sequences are available in the conservation artifact directory.</p>
</main>
<script id="geison-report-data" type="application/json">__REPORT_DATA__</script>
<script>
"use strict";
const reportData=JSON.parse(document.getElementById("geison-report-data").textContent);
const canvas=document.getElementById("conservation-canvas");
const context=canvas.getContext("2d");
const resetButton=document.getElementById("reset-zoom");
const windows=reportData.windows;
const annotations=reportData.annotations;
const margins={left:58,right:20,top:22,bottom:70};
const genomeStart=windows.length?windows[0][0]:1;
const genomeEnd=windows.length?Math.max(...windows.map(item=>item[1])):1;
let viewStart=genomeStart;
let viewEnd=genomeEnd;
let cssWidth=1;
let cssHeight=1;
let dragging=false;
let dragX=0;
let dragStart=1;
let dragEnd=1;

document.getElementById("target").textContent=reportData.identity.targetName;
document.getElementById("reference").textContent=reportData.identity.referenceId??"None";
document.getElementById("sequence-count").textContent=String(reportData.identity.sequenceCount);
document.getElementById("window-summary").textContent=`${reportData.identity.windowCount} windows; ${reportData.identity.windowSize} bases, step ${reportData.identity.stepSize}`;

function formatFraction(value){return Number(value).toFixed(3)}
function formatCoordinate(value){return Math.round(value).toLocaleString("en-US")}
function plotLeft(){return margins.left}
function plotRight(){return Math.max(margins.left+1,cssWidth-margins.right)}
function plotTop(){return margins.top}
function annotationTop(){return Math.max(margins.top+60,cssHeight-margins.bottom+22)}
function plotBottom(){return annotationTop()-26}
function xFor(referencePosition){return plotLeft()+((referencePosition-viewStart)/Math.max(1,viewEnd-viewStart))*(plotRight()-plotLeft())}
function yFor(value){return plotBottom()-Number(value)*(plotBottom()-plotTop())}

function clampView(start,end){
  const fullWidth=Math.max(0,genomeEnd-genomeStart);
  const width=Math.min(Math.max(0,end-start),fullWidth);
  let nextStart=start;
  if(nextStart<genomeStart)nextStart=genomeStart;
  if(nextStart+width>genomeEnd)nextStart=genomeEnd-width;
  viewStart=nextStart;
  viewEnd=nextStart+width;
}

function visibleWindows(){return windows.filter(item=>item[1]>=viewStart&&item[0]<=viewEnd)}

function drawTrace(items,index,color){
  if(!items.length)return;
  context.beginPath();
  items.forEach((item,offset)=>{
    const center=(item[0]+item[1])/2;
    const x=xFor(center);
    const y=yFor(item[index]);
    if(offset===0)context.moveTo(x,y);else context.lineTo(x,y);
  });
  context.strokeStyle=color;
  context.lineWidth=2;
  context.lineJoin="round";
  context.stroke();
}

function draw(){
  context.clearRect(0,0,cssWidth,cssHeight);
  context.fillStyle="#ffffff";
  context.fillRect(0,0,cssWidth,cssHeight);
  document.getElementById("view-interval").textContent=windows.length?`Reference ${formatCoordinate(viewStart)}–${formatCoordinate(viewEnd)}`:"No reference interval";
  if(!windows.length)return;
  const items=visibleWindows();
  const topKeys=new Set(reportData.topWindows.map(item=>`${item[0]}:${item[1]}`));
  for(const item of items){
    if(topKeys.has(`${item[0]}:${item[1]}`)){
      context.fillStyle="#fff0b3";
      context.fillRect(xFor(item[0]),plotTop(),Math.max(1,xFor(item[1])-xFor(item[0])),plotBottom()-plotTop());
    }
  }
  context.strokeStyle="#d8deea";
  context.fillStyle="#596579";
  context.font="12px system-ui,sans-serif";
  context.textAlign="right";
  for(let tick=0;tick<=4;tick+=1){
    const value=tick/4;
    const y=yFor(value);
    context.beginPath();context.moveTo(plotLeft(),y);context.lineTo(plotRight(),y);context.stroke();
    context.fillText(value.toFixed(2),plotLeft()-7,y+4);
  }
  drawTrace(items,3,"#1769aa");
  drawTrace(items,5,"#d65f00");
  context.textAlign="left";
  context.fillStyle="#596579";
  context.fillText(formatCoordinate(viewStart),plotLeft(),plotBottom()+18);
  context.textAlign="right";
  context.fillText(formatCoordinate(viewEnd),plotRight(),plotBottom()+18);
  context.textAlign="left";
  context.fillStyle="#596579";
  context.fillText("Reference annotations",plotLeft(),annotationTop()-7);
  for(const annotation of annotations){
    if(annotation[2]<viewStart||annotation[1]>viewEnd)continue;
    const start=Math.max(annotation[1],viewStart);
    const end=Math.min(annotation[2],viewEnd);
    const x=xFor(start);
    const width=Math.max(2,xFor(end)-x);
    context.fillStyle="#805ad5";
    context.fillRect(x,annotationTop(),width,13);
    if(width>32){
      context.save();context.beginPath();context.rect(x,annotationTop(),width,13);context.clip();
      context.fillStyle="#ffffff";context.font="10px system-ui,sans-serif";context.textAlign="left";
      context.fillText(annotation[4],x+3,annotationTop()+10);context.restore();
    }
  }
}

function resizeCanvas(){
  const bounds=canvas.getBoundingClientRect();
  cssWidth=Math.max(1,bounds.width);
  cssHeight=Math.max(1,bounds.height);
  const ratio=Math.max(1,window.devicePixelRatio||1);
  canvas.width=Math.round(cssWidth*ratio);
  canvas.height=Math.round(cssHeight*ratio);
  context.setTransform(ratio,0,0,ratio,0,0);
  draw();
}

function resetZoom(){viewStart=genomeStart;viewEnd=genomeEnd;draw()}

function onWheel(event){
  if(!windows.length)return;
  event.preventDefault();
  const bounds=canvas.getBoundingClientRect();
  const fraction=Math.min(1,Math.max(0,(event.clientX-bounds.left-plotLeft())/Math.max(1,plotRight()-plotLeft())));
  const width=Math.max(1,viewEnd-viewStart);
  const minimumWidth=Math.min(Math.max(1,reportData.identity.windowSize),Math.max(1,genomeEnd-genomeStart));
  const nextWidth=Math.min(Math.max(minimumWidth,width*(event.deltaY>0?1.25:0.8)),Math.max(1,genomeEnd-genomeStart));
  const anchor=viewStart+fraction*width;
  clampView(anchor-fraction*nextWidth,anchor+(1-fraction)*nextWidth);
  draw();
}

function onPointerDown(event){
  if(!windows.length)return;
  dragging=true;dragX=event.clientX;dragStart=viewStart;dragEnd=viewEnd;
  canvas.setPointerCapture(event.pointerId);
}

function nearestWindow(referencePosition){
  let low=0;let high=windows.length-1;
  while(low<high){const middle=Math.floor((low+high)/2);if(windows[middle][0]<referencePosition)low=middle+1;else high=middle}
  const candidates=[windows[low],windows[Math.max(0,low-1)]].filter(Boolean).filter(item=>item[1]>=viewStart&&item[0]<=viewEnd);
  if(!candidates.length)return null;
  return candidates.reduce((best,item)=>Math.abs((item[0]+item[1])/2-referencePosition)<Math.abs((best[0]+best[1])/2-referencePosition)?item:best);
}

function showWindow(item){
  if(!item)return;
  document.getElementById("detail-interval").textContent=`${formatCoordinate(item[0])}–${formatCoordinate(item[1])}`;
  document.getElementById("detail-mean").textContent=formatFraction(item[3]);
  document.getElementById("detail-minimum").textContent=formatFraction(item[4]);
  document.getElementById("detail-coverage").textContent=formatFraction(item[5]);
  document.getElementById("detail-gaps").textContent=formatFraction(item[6]);
  document.getElementById("detail-entropy").textContent=`${formatFraction(item[7])} bits`;
  const overlapping=annotations.filter(annotation=>annotation[2]>=item[0]&&annotation[1]<=item[1]).map(annotation=>`${annotation[4]} (${annotation[0]}, ${annotation[1]}–${annotation[2]})`);
  document.getElementById("detail-annotations").textContent=overlapping.length?overlapping.join("; "):"None";
}

function onPointerMove(event){
  if(!windows.length)return;
  if(dragging){
    const referencePerPixel=(dragEnd-dragStart)/Math.max(1,plotRight()-plotLeft());
    const shift=(dragX-event.clientX)*referencePerPixel;
    clampView(dragStart+shift,dragEnd+shift);draw();return;
  }
  const bounds=canvas.getBoundingClientRect();
  const x=Math.min(plotRight(),Math.max(plotLeft(),event.clientX-bounds.left));
  const referencePosition=viewStart+((x-plotLeft())/Math.max(1,plotRight()-plotLeft()))*(viewEnd-viewStart);
  showWindow(nearestWindow(referencePosition));
}

function onPointerUp(event){
  if(!dragging)return;
  dragging=false;
  if(canvas.hasPointerCapture(event.pointerId))canvas.releasePointerCapture(event.pointerId);
}

function focusWindow(item){
  const padding=Math.max(reportData.identity.stepSize,(item[1]-item[0]+1)*0.35);
  clampView(item[0]-padding,item[1]+padding);draw();showWindow(item);
  canvas.focus();
}

function populateTopWindows(){
  const body=document.querySelector("#top-windows tbody");
  for(const item of reportData.topWindows){
    const row=document.createElement("tr");row.tabIndex=0;
    const values=[`${formatCoordinate(item[0])}–${formatCoordinate(item[1])}`,formatFraction(item[3]),formatFraction(item[4]),formatFraction(item[5]),formatFraction(item[7])];
    for(const value of values){const cell=document.createElement("td");cell.textContent=value;row.appendChild(cell)}
    row.addEventListener("click",()=>focusWindow(item));
    row.addEventListener("keydown",event=>{if(event.key==="Enter"||event.key===" "){event.preventDefault();focusWindow(item)}});
    body.appendChild(row);
  }
}

canvas.tabIndex=0;
canvas.addEventListener("wheel",onWheel,{passive:false});
canvas.addEventListener("pointerdown",onPointerDown);
canvas.addEventListener("pointermove",onPointerMove);
canvas.addEventListener("pointerup",onPointerUp);
canvas.addEventListener("pointercancel",onPointerUp);
resetButton.addEventListener("click",resetZoom);
populateTopWindows();
if("ResizeObserver" in window)new ResizeObserver(resizeCanvas).observe(canvas);else window.addEventListener("resize",resizeCanvas);
resizeCanvas();
</script>
</body>
</html>
"""
