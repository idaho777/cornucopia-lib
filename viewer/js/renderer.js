// viewer/js/renderer.js
(function (window) {
  const container = document.getElementById("container");

  window.renderCurrent = function (data, sortValue, runName) {
    container.innerHTML = "";
    let sorted = [...data];
    switch (sortValue) {
      case "s":
        sorted.sort((a, b) => a.s[0] - b.s[0] || a.s[1] - b.s[1]);
        break;
      case "T":
        sorted.sort((a, b) => a.T[0] - b.T[0] || a.T[1] - b.T[1]);
        break;
      case "a":
        sorted.sort((a, b) => a.a - b.a);
        break;
      case "noise":
        sorted.sort((a, b) => a.noise - b.noise);
        break;
      // default: leave in original order
    }
    sorted.forEach((group) => window.renderGroup(group, runName));
  };

  window.renderGroup = function (group, runName) {
    const grp = document.createElement("div");
    grp.className = "group";

    // group label
    const lbl = document.createElement("div");
    lbl.className = "label";
    lbl.textContent =
      `Arc Length Range s=[${group.s[0]}, ${group.s[1]}],` +
      `Translation T=(${group.T[0]}, ${group.T[1]}),   ` +
      `Anticlothoid Circle Radius a=${group.a},   ` +
      `Vertex noise=${group.noise}`;
    grp.appendChild(lbl);

    // build image URLs
    const urls = group.images.map(
      (img) => `/runs/${runName}/img/${img.filename}`,
    );

    // build one caption per image, using its own rot
    const captions = group.images.map(
      (img) =>
        `Arc Length: [${group.s[0]} - ${group.s[1]}]<br>` +
        `Translation: (${group.T[0]}, ${group.T[1]})<br>` +
        `Rotation: ${img.rot}r<br>` +
        `Circle Radius: ${group.a}<br>` +
        `Noise: ${group.noise}`,
    );

    const row = document.createElement("div");
    row.className = "row";

    group.images.forEach((img, idx) => {
      const wrapper = document.createElement("div");
      wrapper.className = "img-container";

      const imgEl = document.createElement("img");
      imgEl.src = urls[idx];
      imgEl.title = `Rot = ${img.rot}`;
      imgEl.dataset.idx = idx;
      imgEl.addEventListener("click", (e) => {
        e.stopPropagation();
        window.showOverlay(urls, captions, +e.target.dataset.idx);
      });

      const lbl2 = document.createElement("div");
      lbl2.className = "img-label";
      lbl2.textContent = `Rot = ${img.rot}`;

      wrapper.append(imgEl, lbl2);
      row.appendChild(wrapper);
    });

    grp.appendChild(row);
    container.appendChild(grp);
  };
})(window);
