const SEGMENT_NAMES = [
  "N-term",
  "CX",
  "N-loop",
  "B1",
  "30s-loop",
  "B2",
  "40s-loop",
  "B3",
  "50s-loop",
  "Helix",
  "C-term"
];

document.addEventListener("DOMContentLoaded", () => {
  const segmenter = document.getElementById("segmenter");
  if (!segmenter) return;
  const sequence = segmenter.dataset.seq;
  const chars = sequence.split("");

  const container = document.createElement("div");
  container.id = "sequence-container";
  container.style.position = "relative";
  container.style.fontFamily = "monospace";
  container.style.whiteSpace = "pre";
  container.style.userSelect = "none";

  container.style.marginBottom = "20px";

  segmenter.appendChild(container);

  chars.forEach(ch => {
    const span = document.createElement("span");
    span.textContent = ch;
    container.appendChild(span);
  });

  const charWidth = container.firstChild.getBoundingClientRect().width;
  container.style.width = charWidth * chars.length + "px";

  const separators = [];
  const labels = [];
  const totalSeparators = SEGMENT_NAMES.length - 1;
  for (let i = 1; i <= totalSeparators; i++) {
    const sep = document.createElement("div");
    sep.className = "separator";
    sep.style.position = "absolute";
    sep.style.top = 0;
    sep.style.bottom = 0;
    sep.style.width = "2px";
    sep.style.background = "red";
    sep.style.cursor = "ew-resize";
    const pos = Math.round((i * chars.length) / SEGMENT_NAMES.length);
    sep.style.left = pos * charWidth + "px";
    sep.dataset.index = i - 1;
    container.appendChild(sep);

    const handle = document.createElement("div");
    handle.className = "separator-handle";
    handle.style.position = "absolute";
    handle.style.bottom = "-14px";
    handle.style.left = "-8px";
    handle.style.width = "0";
    handle.style.height = "0";
    handle.style.borderLeft = "8px solid transparent";
    handle.style.borderRight = "8px solid transparent";
    handle.style.borderTop = "14px solid red";
    handle.style.cursor = "ew-resize";
    sep.appendChild(handle);

    separators.push(sep);
  }

  // create labels for each segment
  SEGMENT_NAMES.forEach(name => {
    const label = document.createElement("div");
    label.className = "segment-label";
    label.textContent = name;
    label.style.position = "absolute";
    label.style.top = "-1.5em";
    label.style.fontSize = "0.8em";
    label.style.pointerEvents = "none";
    container.appendChild(label);
    labels.push(label);
  });

  function updateLabels() {
    const positions = [0, ...separators.map(sep => parseFloat(sep.style.left)), container.clientWidth];
    labels.forEach((label, i) => {
      const start = positions[i];
      const end = positions[i + 1];
      const center = (start + end) / 2 - label.offsetWidth / 2;
      label.style.left = center + "px";
    });
  }

  updateLabels();

  let active = null;
  container.addEventListener("mousedown", e => {
    if (e.target.classList.contains("separator")) {
      active = e.target;
      e.preventDefault();
    } else if (e.target.parentElement && e.target.parentElement.classList.contains("separator")) {
      active = e.target.parentElement;
      e.preventDefault();
    }
  });
  document.addEventListener("mousemove", e => {
    if (!active) return;
    const rect = container.getBoundingClientRect();
    let x = e.clientX - rect.left;
    const index = parseInt(active.dataset.index);
    const min = index === 0 ? 0 : parseFloat(separators[index - 1].style.left) + charWidth;
    const max = index === separators.length - 1 ? rect.width : parseFloat(separators[index + 1].style.left) - charWidth;
    if (x < min) x = min;
    if (x > max) x = max;
    x = Math.round(x / charWidth) * charWidth;
    if (x < min) x = min; // ensure snapping doesn't violate bounds
    if (x > max) x = max;
    active.style.left = x + "px";
    updateLabels();
  });
  document.addEventListener("mouseup", () => {
    active = null;
  });

  const form = document.getElementById("segment-form");
  form.addEventListener("submit", () => {
    const segmentsInput = document.getElementById("id_segments");
    const positions = separators
      .map(sep => parseFloat(sep.style.left) / charWidth)
      .map(pos => Math.round(pos));
    let prev = 1;
    const ranges = [];
    positions.forEach((p, idx) => {
      const end = p;
      ranges.push(`${prev}-${end}:${SEGMENT_NAMES[idx]}`);
      prev = end + 1;
    });
    ranges.push(`${prev}-${chars.length}:${SEGMENT_NAMES[SEGMENT_NAMES.length - 1]}`);
    segmentsInput.value = ranges.join(",");
  });
  
});