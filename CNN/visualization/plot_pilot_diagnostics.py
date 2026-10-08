#!/usr/bin/env python3
"""Generate detailed diagnostic plots from the no-Poisson pilot checkpoint."""

from __future__ import annotations

import argparse
import base64
import io
import json
import os
from pathlib import Path
from typing import Any

os.environ.setdefault("MPLCONFIGDIR", "/tmp/cnn21d-matplotlib")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np
import torch
from PIL import Image

from training.cnn_dataset import CNNDataConfig
from training.cnn21d.metrics import slopes_to_angles
from training.train_cnn21d import build_model, load_config, make_loader, sample_types_to_labels


def collect_split(
    split: str,
    config: dict[str, Any],
    model: torch.nn.Module,
    device: torch.device,
    mode: str = "pilot",
) -> list[dict[str, Any]]:
    seed = int(config["seed"])
    task = config["task"]
    positive_types = tuple(int(x) for x in task["positive_sample_types"])
    include_types = tuple(int(x) for x in task["include_sample_types"])
    theta_values = task.get("signal_theta_range_mrad")
    signal_theta_range_mrad = None if theta_values is None else (float(theta_values[0]), float(theta_values[1]))
    data_config = CNNDataConfig(**config["data"]["dataset"], seed=seed)
    loader = make_loader(
        Path(config["data"]["hdf5_path"]).resolve(),
        split,
        data_config,
        batch_size=int(config["training"]["batch_size"]),
        workers=0,
        seed=seed + {"train": 0, "validation": 1, "test": 2}[split],
        subset_size=config["subsets"][mode][split],
        training=False,
        shuffle=False,
        include_sample_types=include_types,
        positive_sample_types=positive_types,
        signal_theta_range_mrad=signal_theta_range_mrad,
    )
    records: list[dict[str, Any]] = []
    model.eval()
    with torch.no_grad():
        for batch in loader:
            volume = batch["volume"].to(device)
            output = model(volume)
            probability = torch.sigmoid(output["presence_logit"]).cpu().numpy()
            slope_prediction = output["slope_xy"].cpu().numpy()
            slope_truth = batch["slope_xy"].numpy()
            sample_type = batch["sample_type"].numpy()
            labels = sample_types_to_labels(batch["sample_type"], positive_types).numpy()
            activation = volume.clamp_min(0).sum(dim=(1, 2)).cpu()
            rim = activation.clone(); rim[:, 3:-3, 3:-3] = 0
            border_fraction = (rim.sum((1, 2)) / activation.sum((1, 2)).clamp_min(1e-12)).numpy()
            hdf5_index = batch["metadata"]["hdf5_index"].numpy()
            dataset_index = batch["metadata"]["dataset_index"].numpy()
            volumes = volume.cpu().numpy()
            theta_true, _ = slopes_to_angles(slope_truth)
            theta_pred, _ = slopes_to_angles(slope_prediction)
            for index in range(len(probability)):
                records.append({
                    "split": split,
                    "dataset_index": int(dataset_index[index]),
                    "hdf5_index": int(hdf5_index[index]),
                    "sample_type": int(sample_type[index]),
                    "label": int(labels[index]),
                    "probability": float(probability[index]),
                    "predicted_label": int(probability[index] >= config["evaluation"]["threshold"]),
                    "slope_truth": slope_truth[index].astype(float),
                    "slope_prediction": slope_prediction[index].astype(float),
                    "theta_true_mrad": float(theta_true[index]),
                    "theta_pred_mrad": float(theta_pred[index]),
                    "theta_abs_error_mrad": float(abs(theta_pred[index] - theta_true[index])),
                    "border_fraction": float(border_fraction[index]),
                    "volume": volumes[index, 0],
                })
    return records


def save_angle_scatter(records: list[dict[str, Any]], path: Path) -> None:
    figure, axes = plt.subplots(1, 2, figsize=(11, 4.8), constrained_layout=True)
    colors = {0: "tab:blue", 1: "tab:orange"}
    names = {0: "pred hard", 1: "pred signal"}
    signal = [record for record in records if record["sample_type"] == 2]
    maximum = max(max(record["theta_true_mrad"], record["theta_pred_mrad"]) for record in signal)
    for axis, split in zip(axes, ("validation", "test"), strict=True):
        selected = [record for record in signal if record["split"] == split]
        for predicted in (0, 1):
            group = [record for record in selected if record["predicted_label"] == predicted]
            axis.scatter(
                [record["theta_true_mrad"] for record in group],
                [record["theta_pred_mrad"] for record in group],
                s=30, alpha=0.8, color=colors[predicted], label=f"{names[predicted]} (n={len(group)})",
            )
        axis.plot([0, maximum], [0, maximum], "--", color="0.4", linewidth=1, label="pred=true")
        axis.set(xlabel="theta true [mrad]", ylabel="theta predetto [mrad]", title=split.capitalize())
        axis.grid(alpha=0.2); axis.legend(fontsize=8)
    figure.suptitle("Signal del pilot: ricostruzione angolare e decisione hard/signal")
    figure.savefig(path, dpi=170); plt.close(figure)


def save_layer_sheet(record: dict[str, Any], path: Path, *, show_slope_arrows: bool = False) -> None:
    volume = record["volume"]
    vmin, vmax = np.quantile(volume, [0.02, 0.995])
    figure, axes = plt.subplots(8, 8, figsize=(11, 11), constrained_layout=True)
    for z, axis in enumerate(axes.flat):
        axis.axis("off")
        if z < 57:
            axis.imshow(volume[z], origin="lower", cmap="magma", vmin=vmin, vmax=vmax, interpolation="nearest")
            if show_slope_arrows:
                center = (9.5, 9.5); scale = 4.0
                truth = np.asarray(record["slope_truth"]) * scale
                prediction = np.asarray(record["slope_prediction"]) * scale
                axis.arrow(*center, truth[0], truth[1], color="cyan", width=0.10, head_width=0.65, length_includes_head=True)
                axis.arrow(*center, prediction[0], prediction[1], color="lime", width=0.10, head_width=0.65, length_includes_head=True)
            axis.set_title(f"z={z:02d}", fontsize=7)
    if show_slope_arrows:
        title = (
            f"{record['split']} signal | HDF5 {record['hdf5_index']} | "
            f"theta true={record['theta_true_mrad']:.1f}, pred={record['theta_pred_mrad']:.1f}, "
            f"|err|={record['theta_abs_error_mrad']:.1f} mrad | cyan=vera, verde=predetta"
        )
    else:
        title = f"{record['split']} hard→signal | HDF5 {record['hdf5_index']} | P(signal)={record['probability']:.3f}"
    figure.suptitle(title, fontsize=12)
    figure.savefig(path, dpi=130); plt.close(figure)


def draw_projection(axis: plt.Axes, record: dict[str, Any], *, show_rim: bool = False) -> None:
    projection = np.maximum(record["volume"], 0).sum(axis=0)
    vmax = np.quantile(projection, 0.995)
    axis.imshow(projection, origin="lower", cmap="magma", vmin=0, vmax=max(vmax, 1e-6))
    center = np.asarray([9.5, 9.5])
    scale = 10.0
    truth = np.asarray(record["slope_truth"]) * scale
    prediction = np.asarray(record["slope_prediction"]) * scale
    axis.arrow(*center, truth[0], truth[1], color="cyan", width=0.08, head_width=0.55, length_includes_head=True)
    axis.arrow(*center, prediction[0], prediction[1], color="lime", width=0.08, head_width=0.55, length_includes_head=True)
    if show_rim:
        axis.add_patch(Rectangle((2.5, 2.5), 14, 14, fill=False, edgecolor="white", linewidth=1.2, linestyle="--"))
    axis.set_xticks([]); axis.set_yticks([])


def save_large_errors(records: list[dict[str, Any]], path: Path) -> list[dict[str, Any]]:
    chosen: list[dict[str, Any]] = []
    for split in ("validation", "test"):
        signal = sorted(
            (record for record in records if record["split"] == split and record["sample_type"] == 2),
            key=lambda record: record["theta_abs_error_mrad"], reverse=True,
        )
        chosen.extend(signal[:3])
    figure, axes = plt.subplots(2, 3, figsize=(10, 7), constrained_layout=True)
    for axis, record in zip(axes.flat, chosen, strict=True):
        draw_projection(axis, record)
        axis.set_title(
            f"{record['split']} HDF5 {record['hdf5_index']}\n"
            f"true={record['theta_true_mrad']:.1f}, pred={record['theta_pred_mrad']:.1f}, |err|={record['theta_abs_error_mrad']:.1f} mrad\n"
            f"P(signal)={record['probability']:.3f}", fontsize=8,
        )
    figure.suptitle("Signal con maggiore errore angolare — cyan: slope vera, verde: predetta")
    figure.savefig(path, dpi=170); plt.close(figure)
    return chosen


def save_central_edge(records: list[dict[str, Any]], path: Path) -> list[dict[str, Any]]:
    validation_signal = [record for record in records if record["split"] == "validation" and record["sample_type"] == 2]
    central = min(validation_signal, key=lambda record: record["border_fraction"])
    edge = max(validation_signal, key=lambda record: record["border_fraction"])
    figure, axes = plt.subplots(1, 2, figsize=(8, 4.5), constrained_layout=True)
    for axis, record, name in zip(axes, (central, edge), ("Più centrale", "Più al bordo"), strict=True):
        draw_projection(axis, record, show_rim=True)
        axis.set_title(
            f"{name} — HDF5 {record['hdf5_index']}\n"
            f"frazione bordo={record['border_fraction']:.3f}, errore theta={record['theta_abs_error_mrad']:.1f} mrad\n"
            f"P(signal)={record['probability']:.3f}", fontsize=9,
        )
    figure.suptitle("Definizione diagnostica central/edge — linea tratteggiata: limite della cornice di 3 pixel")
    figure.savefig(path, dpi=180); plt.close(figure)
    return [central, edge]


def image_data_url(path: Path) -> str:
    # Keep the in-conversation visualization below 1 MB while retaining the
    # lossless PNG files on disk for detailed inspection.
    with Image.open(path) as source:
        image = source.convert("RGB")
        is_layer_sheet = "57_layers" in str(path)
        maximum = 500 if is_layer_sheet else 820
        image.thumbnail((maximum, maximum), Image.Resampling.LANCZOS)
        buffer = io.BytesIO()
        image.save(buffer, format="JPEG", quality=28 if is_layer_sheet else 45, optimize=True)
    return "data:image/jpeg;base64," + base64.b64encode(buffer.getvalue()).decode("ascii")


def save_visualization_html(
    path: Path,
    records: list[dict[str, Any]],
    false_positive_paths: list[tuple[dict[str, Any], Path]],
    large_error_layer_paths: list[tuple[dict[str, Any], Path]],
    large_errors_path: Path,
    central_edge_path: Path,
) -> None:
    angle_data = [
        {
            "split": record["split"], "true": round(record["theta_true_mrad"], 3),
            "pred": round(record["theta_pred_mrad"], 3), "cls": record["predicted_label"],
            "p": round(record["probability"], 4), "id": record["hdf5_index"],
        }
        for record in records if record["sample_type"] == 2
    ]
    sheet_data = [
        {"label": f"{record['split']} · HDF5 {record['hdf5_index']} · P={record['probability']:.3f}", "src": image_data_url(image_path)}
        for record, image_path in false_positive_paths
    ]
    error_sheet_data = [
        {"label": f"{record['split']} · HDF5 {record['hdf5_index']} · |err|={record['theta_abs_error_mrad']:.1f} mrad", "src": image_data_url(image_path)}
        for record, image_path in large_error_layer_paths
    ]
    hard_control = (
        '<select id="fp-select" aria-label="Seleziona falso positivo"></select><img id="fp-image" alt="Contact sheet dei 57 layer dell hard classificato come signal">'
        if sheet_data else
        '<div>Nessun hard è stato classificato come signal alla soglia 0.5.</div><select id="fp-select" aria-label="Nessun falso positivo" style="display:none"></select><img id="fp-image" alt="Nessun falso positivo" style="display:none">'
    )
    fragment = f'''<div id="pilot-diag-root">
  <style>
    #pilot-diag-root {{ color: var(--foreground); font-family: var(--font-sans); display:grid; gap:28px; }}
    #pilot-diag-root h2 {{ margin:0 0 8px; font-size:18px; }}
    #pilot-diag-root .panels {{ display:grid; grid-template-columns:repeat(2,minmax(0,1fr)); gap:18px; }}
    #pilot-diag-root svg {{ width:100%; min-height:330px; overflow:visible; }}
    #pilot-diag-root .frame {{ fill:transparent; stroke:var(--border); }}
    #pilot-diag-root text {{ fill:var(--foreground); font-size:12px; }}
    #pilot-diag-root .legend {{ display:flex; gap:16px; flex-wrap:wrap; margin-bottom:5px; }}
    #pilot-diag-root .legend button {{ border:0; background:transparent; color:var(--foreground); padding:2px 0; }}
    #pilot-diag-root .sw {{ display:inline-block; width:10px; height:10px; margin-right:5px; }}
    #pilot-diag-root select {{ color:var(--foreground); background:var(--background); border:1px solid var(--border); padding:5px 7px; max-width:100%; }}
    #pilot-diag-root img {{ display:block; max-width:100%; max-height:740px; margin:8px auto 0; object-fit:contain; }}
    #pilot-diag-root .tooltip {{ position:absolute; pointer-events:none; display:none; background:var(--popover); color:var(--popover-foreground); border:1px solid var(--border); padding:6px 8px; font-size:12px; z-index:5; }}
    @media(max-width:600px) {{ #pilot-diag-root .panels {{ grid-template-columns:1fr; }} }}
  </style>
  <section><h2>Signal: theta vero vs predetto</h2><div class="legend"><button type="button" aria-pressed="true"><span class="sw" style="background:var(--viz-series-1)"></span>pred hard</button><button type="button" aria-pressed="true"><span class="sw" style="background:var(--viz-series-2)"></span>pred signal</button></div><div class="panels"><svg data-split="validation"></svg><svg data-split="test"></svg></div></section>
  <section><h2>Hard classificati come signal: tutti i 57 layer</h2>{hard_control}</section>
  <section><h2>Signal con maggiore errore sulla regressione angolare</h2><img src="{image_data_url(large_errors_path)}" alt="Sei signal con maggiore errore angolare, con slope vera e predetta"></section>
  <section><h2>Signal con errore maggiore: tutti i 57 layer e slope</h2><select id="error-select" aria-label="Seleziona signal con errore angolare grande"></select><img id="error-image" alt="Contact sheet dei 57 layer con slope vera e predetta"></section>
  <section><h2>Esempio central vs edge</h2><img src="{image_data_url(central_edge_path)}" alt="Confronto tra un signal centrale e uno al bordo"></section>
  <div class="tooltip" role="tooltip"></div>
  <script src="https://cdn.jsdelivr.net/npm/d3@7.9.0/dist/d3.min.js"></script>
  <script>
  (() => {{
    const root=document.getElementById('pilot-diag-root');
    const data={json.dumps(angle_data, separators=(',', ':'))};
    const sheets={json.dumps(sheet_data, separators=(',', ':'))};
    const errorSheets={json.dumps(error_sheet_data, separators=(',', ':'))};
    const select=root.querySelector('#fp-select'), image=root.querySelector('#fp-image');
    sheets.forEach((d,i)=>{{ const o=document.createElement('option'); o.value=i; o.textContent=d.label; select.appendChild(o); }});
    function showSheet() {{ const d=sheets[+select.value || 0]; if(d) {{ image.src=d.src; image.alt=d.label; }} }}
    select.addEventListener('change',showSheet); showSheet();
    const errorSelect=root.querySelector('#error-select'), errorImage=root.querySelector('#error-image');
    errorSheets.forEach((d,i)=>{{ const o=document.createElement('option'); o.value=i; o.textContent=d.label; errorSelect.appendChild(o); }});
    function showErrorSheet() {{ const d=errorSheets[+errorSelect.value || 0]; if(d) {{ errorImage.src=d.src; errorImage.alt=d.label; }} }}
    errorSelect.addEventListener('change',showErrorSheet); showErrorSheet();
    const tooltip=root.querySelector('.tooltip');
    function draw(svgNode) {{
      const split=svgNode.dataset.split, values=data.filter(d=>d.split===split);
      const width=Math.max(360,svgNode.parentElement.clientWidth), height=350;
      const m={{top:28,right:18,bottom:52,left:66}}, innerW=width-m.left-m.right, innerH=height-m.top-m.bottom;
      const svg=d3.select(svgNode).attr('viewBox',`0 0 ${{width}} ${{height}}`); svg.selectAll('*').remove();
      const max=d3.max(values,d=>Math.max(d.true,d.pred))*1.05;
      const x=d3.scaleLinear().domain([0,max]).nice().range([m.left,m.left+innerW]);
      const y=d3.scaleLinear().domain([0,max]).nice().range([m.top+innerH,m.top]);
      svg.append('rect').attr('class','frame').attr('data-chart-frame','').attr('x',m.left).attr('y',m.top).attr('width',innerW).attr('height',innerH);
      svg.append('g').attr('transform',`translate(0,${{m.top+innerH}})`).call(d3.axisBottom(x).ticks(width<500?4:6));
      svg.append('g').attr('transform',`translate(${{m.left}},0)`).call(d3.axisLeft(y).ticks(5));
      svg.append('line').attr('x1',x(0)).attr('y1',y(0)).attr('x2',x(max)).attr('y2',y(max)).attr('stroke','var(--border)').attr('stroke-dasharray','5 4');
      svg.append('text').attr('x',width/2).attr('y',18).attr('text-anchor','middle').attr('font-weight',600).text(split[0].toUpperCase()+split.slice(1));
      svg.append('text').attr('class','axis-title').attr('data-axis','x').attr('x',m.left+innerW/2).attr('y',height-8).attr('text-anchor','middle').text('theta vero [mrad]');
      svg.append('text').attr('class','axis-title').attr('data-axis','y').attr('transform',`translate(16,${{m.top+innerH/2}}) rotate(-90)`).attr('text-anchor','middle').text('theta predetto [mrad]');
      svg.selectAll('circle').data(values).join('circle').attr('cx',d=>x(d.true)).attr('cy',d=>y(d.pred)).attr('r',4).attr('fill',d=>d.cls?'var(--viz-series-2)':'var(--viz-series-1)').attr('opacity',.82)
        .on('pointerenter',(event,d)=>{{ tooltip.style.display='block'; tooltip.textContent=`HDF5 ${{d.id}} · true ${{d.true}} · pred ${{d.pred}} mrad · P(signal) ${{d.p}}`; }})
        .on('pointermove',event=>{{ const b=root.getBoundingClientRect(); tooltip.style.left=`${{event.clientX-b.left+10}}px`; tooltip.style.top=`${{event.clientY-b.top+10}}px`; }})
        .on('pointerleave',()=>tooltip.style.display='none');
    }}
    const redraw=()=>root.querySelectorAll('svg[data-split]').forEach(draw); redraw(); new ResizeObserver(redraw).observe(root);
  }})();
  </script>
</div>'''
    path.write_text(fragment, encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, default=Path("training/configs/cnn21d.yaml"))
    parser.add_argument("--checkpoint", type=Path, default=Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/best_model.pt"))
    parser.add_argument("--output-dir", type=Path, default=Path("runs/cnn21d_sample_type_no_poisson_theta10_50/pilot/diagnostics"))
    parser.add_argument("--visualization-html", type=Path)
    args = parser.parse_args()
    config = load_config(args.config); device = torch.device("cpu")
    model = build_model(config).to(device)
    checkpoint = torch.load(args.checkpoint, map_location=device, weights_only=False)
    model.load_state_dict(checkpoint["model_state"])
    records = collect_split("validation", config, model, device) + collect_split("test", config, model, device)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    save_angle_scatter(records, args.output_dir / "signal_angles_validation_test.png")
    false_positives = [record for record in records if record["sample_type"] == 1 and record["predicted_label"] == 1]
    false_positive_paths: list[tuple[dict[str, Any], Path]] = []
    layer_dir = args.output_dir / "hard_pred_signal_57_layers"; layer_dir.mkdir(exist_ok=True)
    for record in false_positives:
        path = layer_dir / f"{record['split']}_hdf5_{record['hdf5_index']}.png"
        save_layer_sheet(record, path); false_positive_paths.append((record, path))
    large = save_large_errors(records, args.output_dir / "signal_large_angle_errors.png")
    error_layer_dir = args.output_dir / "signal_large_error_57_layers"; error_layer_dir.mkdir(exist_ok=True)
    large_error_layer_paths: list[tuple[dict[str, Any], Path]] = []
    for record in large:
        path = error_layer_dir / f"{record['split']}_hdf5_{record['hdf5_index']}.png"
        save_layer_sheet(record, path, show_slope_arrows=True)
        large_error_layer_paths.append((record, path))
    central_edge = save_central_edge(records, args.output_dir / "central_vs_edge.png")
    manifest = {
        "false_positive_hard_count": len(false_positives),
        "false_positive_hard": [{key: value for key, value in record.items() if key != "volume"} for record in false_positives],
        "large_angle_errors": [{key: value for key, value in record.items() if key != "volume"} for record in large],
        "central_edge": [{key: value for key, value in record.items() if key != "volume"} for record in central_edge],
    }
    with (args.output_dir / "manifest.json").open("w", encoding="utf-8") as stream:
        json.dump(manifest, stream, indent=2, default=lambda value: value.tolist() if isinstance(value, np.ndarray) else value)
    if args.visualization_html:
        args.visualization_html.parent.mkdir(parents=True, exist_ok=True)
        save_visualization_html(
            args.visualization_html, records, false_positive_paths,
            large_error_layer_paths, args.output_dir / "signal_large_angle_errors.png",
            args.output_dir / "central_vs_edge.png",
        )
    print(json.dumps({"records": len(records), "false_positive_hard": len(false_positives), "output_dir": str(args.output_dir.resolve())}, indent=2))


if __name__ == "__main__":
    main()
