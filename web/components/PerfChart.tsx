"use client";
import { useEffect, useRef } from "react";

export default function PerfChart() {
  const canvasRef = useRef<HTMLCanvasElement>(null);

  useEffect(() => {
    const c = canvasRef.current;
    if (!c) return;
    const ctx = c.getContext("2d");
    if (!ctx) return;

    const W = c.width = (c.parentElement?.clientWidth || 720);
    const H = c.height = 420;
    const pad = { top: 30, right: 30, bottom: 60, left: 70 };
    const pw = W - pad.left - pad.right;
    const ph = H - pad.top - pad.bottom;

    const maxMb = 100;
    const xLabels = [0, 1, 5, 10, 25, 50, 100];
    const xToPx = (mb: number) => pad.left + (mb / maxMb) * pw;

    const maxTimeS = 360000;
    const yLogMax = Math.log10(maxTimeS);
    const yToPx = (t: number) => pad.top + ph - ((Math.log10(Math.max(t, 1)) / yLogMax) * ph);

    // Baseline: O(L) linear, hits ~5Mb memory wall
    const baselineTimes: [number, number][] = [];
    for (let mb = 0.1; mb <= 100; mb *= 1.15) {
      const t = mb <= 5 ? mb * 100 : Infinity;
      if (t <= 500000) baselineTimes.push([mb, t]);
    }
    baselineTimes.push([5.5, 550]);
    baselineTimes.push([6, 600]);

    // Stitched: O(L/k) wall time, k chunks of ~400Kb
    const stitchedTimes: [number, number][] = [];
    const chunkSizeKb = 400;
    const timePerKb = 0.25;
    const stitchOverhead = 5;
    for (let mb = 0.1; mb <= 100; mb *= 1.15) {
      const kb = mb * 1000;
      const numChunks = Math.max(1, Math.ceil(kb / chunkSizeKb));
      const perChunkTime = timePerKb * chunkSizeKb;
      const wallTime = perChunkTime + stitchOverhead + (numChunks > 1 ? 10 : 0);
      stitchedTimes.push([mb, wallTime]);
    }

    // Canvas background
    ctx.fillStyle = "#0f172a";
    ctx.fillRect(0, 0, W, H);

    // Grid
    ctx.strokeStyle = "#1e293b";
    ctx.lineWidth = 0.5;
    for (let pow = 0; pow <= 5; pow++) {
      const t = Math.pow(10, pow);
      const y = yToPx(t);
      ctx.beginPath();
      ctx.moveTo(pad.left, y);
      ctx.lineTo(W - pad.right, y);
      ctx.stroke();

      ctx.fillStyle = "#64748b";
      ctx.font = "11px Inter, sans-serif";
      ctx.textAlign = "right";
      const label = t >= 3600 ? `${(t / 3600).toFixed(0)}h`
        : t >= 60 ? `${(t / 60).toFixed(0)}min`
        : `${t.toFixed(0)}s`;
      ctx.fillText(label, pad.left - 8, y + 4);
    }

    // Y-axis label
    ctx.save();
    ctx.translate(16, pad.top + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.fillStyle = "#94a3b8";
    ctx.font = "13px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Wall time (log scale)", 0, 0);
    ctx.restore();

    // X-axis labels
    ctx.fillStyle = "#64748b";
    ctx.font = "11px Inter, sans-serif";
    ctx.textAlign = "center";
    for (const mb of xLabels) {
      const x = xToPx(mb);
      ctx.fillText(mb === 0 ? "0" : `${mb}Mb`, x, H - pad.bottom + 20);
      ctx.beginPath();
      ctx.moveTo(x, pad.top + ph);
      ctx.lineTo(x, pad.top + ph + 5);
      ctx.strokeStyle = "#334155";
      ctx.stroke();
    }
    ctx.fillText("Region size", pad.left + pw / 2, H - 6);

    // Draw lines helper
    const drawLine = (points: [number, number][], color: string, dash: number[], width: number) => {
      if (points.length < 2) return;
      ctx.save();
      ctx.strokeStyle = color;
      ctx.lineWidth = width;
      ctx.setLineDash(dash);
      ctx.lineJoin = "round";
      ctx.beginPath();
      const first = points[0];
      ctx.moveTo(xToPx(first[0]), yToPx(first[1]));
      let cliff = false;
      for (let i = 1; i < points.length; i++) {
        const pt = points[i];
        if (pt[1] >= 500000) {
          cliff = true;
          ctx.stroke();
          ctx.setLineDash([6, 8]);
          ctx.strokeStyle = color + "80";
          ctx.beginPath();
          const prev = points[i - 1];
          ctx.moveTo(xToPx(prev[0]), yToPx(prev[1]));
          ctx.lineTo(xToPx(pt[0]), pad.top + 2);
          ctx.stroke();
          break;
        }
        ctx.lineTo(xToPx(pt[0]), yToPx(pt[1]));
      }
      if (!cliff) ctx.stroke();
      ctx.restore();
    };

    // Baseline line (red, slopes up then hits wall)
    drawLine(baselineTimes, "#ef4444", [], 2.5);

    // Stitched line (green, near-flat)
    drawLine(stitchedTimes, "#22c55e", [], 2.5);

    // OOM cliff vertical marker at ~5Mb
    ctx.strokeStyle = "#ef444460";
    ctx.lineWidth = 1;
    ctx.setLineDash([3, 5]);
    const wallX = xToPx(5);
    ctx.beginPath();
    ctx.moveTo(wallX, pad.top);
    ctx.lineTo(wallX, pad.top + ph);
    ctx.stroke();
    ctx.setLineDash([]);
    ctx.fillStyle = "#ef4444";
    ctx.font = "11px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("single-machine", wallX, pad.top + ph - 2);
    ctx.fillText("memory limit ~5 Mb", wallX, pad.top + ph + 13);

    // Speedup callout
    const callX = xToPx(50);
    const callY = yToPx(stitchedTimes[stitchedTimes.length - 1][1]) - 60;
    ctx.fillStyle = "#22c55e20";
    ctx.beginPath();
    ctx.roundRect(callX - 60, callY - 10, 230, 34, 6);
    ctx.fill();
    ctx.fillStyle = "#22c55e";
    ctx.font = "bold 13px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Stitched: near-constant wall time", callX + 55, callY + 6);
    ctx.fillText("regardless of region size →", callX + 55, callY + 24);

    // Title
    ctx.fillStyle = "#f1f5f9";
    ctx.font = "bold 15px Inter, sans-serif";
    ctx.textAlign = "center";
    ctx.fillText("Time Complexity: Baseline vs Parallel Stitched", W / 2, pad.top - 8);

    ctx.fillStyle = "#94a3b8";
    ctx.font = "11px Inter, sans-serif";
    ctx.fillText(
      "5 haplotypes · PGGB v0.6.0 · All-pairs alignment O(n\u00b2\u00d7L) \u2192 O(L) for n=5 · Chunks: 400 Kb each",
      W / 2, H - 28
    );
  }, []);

  return (
    <canvas ref={canvasRef} style={{ width: "100%", maxWidth: 740, borderRadius: 10, display: "block" }} />
  );
}