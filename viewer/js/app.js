(function () {
  const runSelector = document.getElementById("runSelector");
  const sortSelector = document.getElementById("sortSelector");

  let currentData = [];
  let currentRun = "";

  // 1) get list of runs
  window.fetchRuns((runs) => {
    runs.forEach((rn) => {
      const opt = document.createElement("option");
      opt.value = opt.text = rn;
      runSelector.appendChild(opt);
    });
    // default‐select last
    if (runSelector.options.length > 1) {
      runSelector.selectedIndex = runSelector.options.length - 1;
      currentRun = runSelector.value;
      sortSelector.disabled = false;
      sortSelector.value = "default";
      loadAndRender(runSelector.value);
    }
  });

  // 2) on run change
  runSelector.addEventListener("change", () => {
    currentRun = runSelector.value;
    sortSelector.value = "default";
    sortSelector.disabled = !currentRun;
    loadAndRender(currentRun);
  });

  // 3) on sort change
  sortSelector.addEventListener("change", () => {
    if (currentData.length) {
      window.renderCurrent(currentData, sortSelector.value, currentRun);
    }
  });

  function loadAndRender(runName) {
    if (!runName) {
      document.getElementById("container").innerHTML = "";
      return;
    }
    window.loadRunIndex(runName, (data) => {
      currentData = data;
      window.renderCurrent(currentData, sortSelector.value, runName);
    });
  }
})();
