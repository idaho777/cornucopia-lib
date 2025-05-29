(function (window) {
  const overlay = document.getElementById("overlay");
  const overlayImg = document.getElementById("overlay-img");
  const overlayCaption = document.getElementById("overlay-caption");
  const btnClose = document.getElementById("overlay-close");
  const btnPrev = document.getElementById("overlay-prev");
  const btnNext = document.getElementById("overlay-next");

  let urls = [];
  let captions = [];
  let idx = 0;

  window.showOverlay = function (u, c, i) {
    urls = u;
    captions = c;
    idx = i;
    overlayImg.src = urls[idx];
    overlayCaption.innerHTML = captions[idx];
    overlay.style.display = "flex";
    console.log(captions);
  };

  function hideOverlay() {
    overlay.style.display = "none";
  }
  function showNext() {
    idx = (idx + 1) % urls.length;
    overlayImg.src = urls[idx];
    overlayCaption.innerHTML = captions[idx];
  }
  function showPrev() {
    idx = (idx - 1 + urls.length) % urls.length;
    overlayImg.src = urls[idx];
    overlayCaption.innerHTML = captions[idx];
  }

  btnClose.addEventListener("click", hideOverlay);
  btnNext.addEventListener("click", (e) => {
    e.stopPropagation();
    showNext();
  });
  btnPrev.addEventListener("click", (e) => {
    e.stopPropagation();
    showPrev();
  });
  overlay.addEventListener("click", (e) => {
    if (e.target === overlay) hideOverlay();
  });
  document.addEventListener("keydown", (e) => {
    if (overlay.style.display !== "flex") return;
    if (e.key === "ArrowRight") showNext();
    if (e.key === "ArrowLeft") showPrev();
    if (e.key === "Escape") hideOverlay();
  });
})(window);
