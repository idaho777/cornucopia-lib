(function (window) {
  // path to your runs.json (adjust if needed)
  const runsPath = "../runs/runs.json";

  window.fetchRuns = function (onSuccess) {
    fetch(runsPath)
      .then((res) => res.json())
      .then(onSuccess)
      .catch((err) => console.error("Failed to load runs.json:", err));
  };

  window.loadRunIndex = function (runName, onSuccess) {
    fetch(`../runs/${runName}/index.json`)
      .then((res) => res.json())
      .then(onSuccess)
      .catch((err) => console.error("Failed to load index:", err));
  };
})(window);
